use std::{
    collections::HashSet,
    ops::{Deref, DerefMut},
    thread,
    time::Instant,
};

use clap::{self, Parser};
use crossbeam::channel::{Receiver, Sender};
use gskits::{
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{get_spin_pb, DEFAULT_INTERVAL},
    utils::ScopedTimer,
};
use num_cpus;
use rsabpoa::{
    abpoa::{abpoa_consensus_dna, abpoa_consensus_dna_core, msa_with_adaptive_fwd_rev, AbpoaParam},
    abpoa_sys::{abpoa_free, abpoa_free_para, abpoa_init, abpoa_para_t, abpoa_t},
};
use rust_htslib::bam::{self, header::HeaderRecord, Header, Read, Reader, Record, Records, Writer};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(short = 'i')]
    pub sbr_bam: String,
}

#[derive(Debug, Clone, Default)]
pub struct Subread {
    pub seq: String,
    pub name: String,
    pub channel: u32,
    pub cx: u8,
}

impl Subread {
    pub fn new(seq: String, name: String, channel: u32, cx: u8) -> Self {
        Self {
            seq,
            name,
            channel,
            cx,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ChannelSubreads(pub Vec<Subread>);
impl Deref for ChannelSubreads {
    type Target = Vec<Subread>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for ChannelSubreads {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Subread>> for ChannelSubreads {
    fn from(value: Vec<Subread>) -> Self {
        Self(value)
    }
}

struct ConsensusResult {
    channel_id: u32,
    consensus_str: String,
}

impl ConsensusResult {
    pub fn new(channel_id: u32, consensus_str: String) -> Self {
        Self {
            channel_id,
            consensus_str,
        }
    }
}

pub struct ChannelSbrIter<'a> {
    records: Records<'a, bam::Reader>,
    channel_subreads: Option<ChannelSubreads>,
}

impl<'a> ChannelSbrIter<'a> {
    pub fn new(records: Records<'a, bam::Reader>) -> Self {
        Self {
            records,
            channel_subreads: Some(ChannelSubreads::default()),
        }
    }
}

impl<'a> Iterator for ChannelSbrIter<'a> {
    type Item = ChannelSubreads;
    fn next(&mut self) -> Option<Self::Item> {
        let res = loop {
            let record = self.records.next();
            if record.is_none() {
                break self.channel_subreads.take();
            }
            let record = record.unwrap().unwrap();
            let record_ext = BamRecordExt::new(&record);
            let seq = record_ext.get_seq();
            let name = record_ext.get_qname();
            let cx = record_ext.get_cx().unwrap();
            let channel_idx = record_ext.get_ch().unwrap();
            let subread = Subread::new(seq, name, channel_idx, cx);

            if let Some(channel_subreads) = &self.channel_subreads {
                if !channel_subreads.is_empty()
                    && channel_subreads.first().unwrap().channel != channel_idx
                {
                    let res = self.channel_subreads.take();
                    self.channel_subreads = Some(ChannelSubreads(vec![subread]));
                    break res;
                } else {
                    self.channel_subreads.as_mut().unwrap().push(subread);
                }
            } else {
                unreachable!();
            }
        };

        let distinct_channels = res
            .as_ref()
            .unwrap_or(&ChannelSubreads::default())
            .iter()
            .map(|sbr| sbr.channel)
            .collect::<HashSet<_>>();
        assert!(distinct_channels.len() <= 1);
        return res;
    }
}

fn bam_reader_worker(bam_file: &str, sender: Sender<ChannelSubreads>) {
    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(2).unwrap();

    let mut pre_channel = None;
    let mut channel_subreads = ChannelSubreads::default();

    let mut scoped_timer = ScopedTimer::new();

    let mut _timer = scoped_timer.perform_timing();

    let mut instant = Instant::now();
    // let records = reader.records();

    for record in reader.records() {
        if instant.elapsed().as_secs() > 60 {
            println!(
                "bam_reader_speed:{:.4}iter/sec",
                _timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }
        let record = record.unwrap();
        let record_ext = BamRecordExt::new(&record);
        let seq = record_ext.get_seq();
        let name = record_ext.get_qname();
        let cx = record_ext.get_cx().unwrap();
        let channel_idx = record_ext.get_ch().unwrap();

        if pre_channel.is_some() && *pre_channel.as_ref().unwrap() != channel_idx {
            _timer.done_with_cnt(1);

            match sender.send(channel_subreads) {
                Ok(_) => {}
                Err(err) => {
                    eprintln!("{}", err);
                    return;
                }
            }

            _timer = scoped_timer.perform_timing();
            channel_subreads = ChannelSubreads::default();
        }

        let subread = Subread::new(seq, name, channel_idx, cx);
        channel_subreads.push(subread);
        pre_channel = Some(channel_idx);
    }

    if channel_subreads.len() > 0 {
        match sender.send(channel_subreads) {
            Ok(_) => {}
            Err(err) => {
                eprintln!("{}", err);
            }
        }
    }
}

fn consensus_worker(
    abpoa_param: &AbpoaParam,
    recv: Receiver<ChannelSubreads>,
    send: Sender<ConsensusResult>,
    idx: usize,
) {
    let mut scoped_timer = ScopedTimer::new();
    let mut instant = Instant::now();
    let mut abpt = abpoa_param.to_abpoa_para_t();
    let abpt: *mut abpoa_para_t = &mut abpt as *mut abpoa_para_t;
    let ab = unsafe { abpoa_init() };

    for channel_subreads in recv {
        if instant.elapsed().as_secs() > 60 && idx == 0 {
            println!(
                "consensus_worker:{:.4}iter/sec",
                scoped_timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }
        let mut _timer = scoped_timer.perform_timing();

        let cns_res = consensus_core(Some(ab), abpt, channel_subreads);
        _timer.done_with_cnt(1);
        if let Some(consensus_res) = cns_res {
            match send.send(consensus_res) {
                Ok(_) => {}
                Err(err) => {
                    eprintln!("{}", err);
                    break;
                }
            }
        }
    }
}

fn consensus_core(
    ab: Option<*mut abpoa_t>,
    abpt: *mut abpoa_para_t,
    mut channel_subreads: ChannelSubreads,
) -> Option<ConsensusResult> {
    channel_subreads = filter_and_sort_subreads(channel_subreads);
    let mut seqs = channel_subreads
        .iter_mut()
        .map(|sbr| unsafe { sbr.seq.as_bytes_mut() })
        .collect::<Vec<_>>();
    let clean_up_ab = ab.is_none();
    let ab = if let Some(ab_) = ab {
        ab_
    } else {
        unsafe { abpoa_init() }
    };

    let msa_result = abpoa_consensus_dna_core(ab, abpt, &mut seqs);
    if clean_up_ab {
        unsafe {
            abpoa_free(ab);
        }
    }

    if msa_result.is_none() {
        return None;
    }
    let consensus_seq = msa_result.as_ref().unwrap().max_num_of_seq_consensus_seq();
    if consensus_seq.is_none() {
        return None;
    }

    let consensus_res = ConsensusResult::new(
        channel_subreads[0].channel,
        consensus_seq.unwrap().to_string(),
    );
    return Some(consensus_res);
}

fn consensus_writer(recv: Receiver<ConsensusResult>, fname: &str) {
    // 创建 BAM 头部
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");

    header.push_record(&hd);

    // 创建 BAM 写入器
    let mut writer =
        Writer::from_path(fname, &header, bam::Format::Bam).expect("Failed to create BAM writer");

    writer.set_threads(2).unwrap();
    let pbar = get_spin_pb("writing...".to_string(), DEFAULT_INTERVAL);
    let mut scoped_timer = ScopedTimer::new();
    let mut instant = Instant::now();
    for consensus_res in recv {
        if instant.elapsed().as_secs() > 60 {
            println!(
                "consensus_writer:{:.4}iter/sec",
                scoped_timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }

        let _timer = scoped_timer.perform_timing();
        let mut record = Record::new();
        record.set_qname(format!("channel_{}", consensus_res.channel_id).as_bytes());
        record.set(
            format!("channel/{}/poa", consensus_res.channel_id).as_bytes(),
            None,
            consensus_res.consensus_str.as_bytes(),
            &vec![255; consensus_res.consensus_str.len()],
        );

        pbar.inc(1);

        if let Err(e) = writer.write(&record) {
            eprintln!("Error writing to BAM file: {}", e);
            break;
        }

        _timer.done_with_cnt(1);
    }
    pbar.finish();
}

fn filter_and_sort_subreads(channel_subreads: ChannelSubreads) -> ChannelSubreads {
    let mut lengths = channel_subreads
        .iter()
        .filter(|sbr| sbr.cx == 3)
        .map(|sbr| sbr.seq.len())
        .collect::<Vec<_>>();

    lengths.sort();

    let max_len = match lengths.len() {
        0 => 1000000,
        1 => lengths[0] * 2,
        n => lengths[n / 2] + lengths[n / 2 - 1],
    };

    let sbureads = channel_subreads.0;
    let mut filtered_subreads = sbureads
        .into_iter()
        .filter(|sbr| sbr.seq.len() < max_len)
        .collect::<Vec<_>>();
    filtered_subreads.sort_by(|left, right| {
        let median = max_len as f32 / 2.0;
        let left_len = left.seq.len() as f32;
        let left_score = (left_len / median).min(median / left_len);

        let right_len = right.seq.len() as f32;
        let right_score = (right_len / median).min(median / right_len);

        right_score.partial_cmp(&left_score).unwrap()
    });

    filtered_subreads.into()
}

fn single_thread(bam_file: &str) {
    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(2).unwrap();
    let channel_sbr_iter = ChannelSbrIter::new(reader.records());

    let param = AbpoaParam::channel_draft_param();
    let mut abpt = param.to_abpoa_para_t();
    let mut succ = 0;
    let mut fail = 0;
    let pbar = get_spin_pb("processing...".to_string(), DEFAULT_INTERVAL);
    let ab = unsafe { abpoa_init() };

    for channel_sbr in channel_sbr_iter {
        if let Some(v) = consensus_core(None, &mut abpt, channel_sbr) {
            succ += 1;
        } else {
            fail += 1;
        }
        pbar.inc(1);
    }
    pbar.finish();
    println!("succ:{succ}, fail:{fail}");

}

fn main() {
    let cli = Cli::parse();
    let inp_bam = &cli.sbr_bam;
    let oup_bam = format!("{}.poa.bam", inp_bam.rsplit_once(".").unwrap().0);
    let oup_bam = &oup_bam;

    single_thread(inp_bam);

    // let threads = cli.threads.unwrap_or(num_cpus::get_physical() / 2);
    // let threads = 1;
    // thread::scope(|thread_scope| {
    //     let (reader_sender, reader_recv) = crossbeam::channel::bounded(1000);
    //     thread_scope.spawn(move || {
    //         bam_reader_worker(&inp_bam, reader_sender);
    //     });

    //     // let mut cnt = 0;
    //     // let instant = Instant::now();
    //     // for v in reader_recv {
    //     //     cnt += 1;
    //     // }
    //     // println!("cnt:{cnt}, {}", instant.elapsed().as_secs());

    //     let (consensus_result_sender, consensus_result_recv) = crossbeam::channel::bounded(1000);
    //     for thread_idx in 0..threads {
    //         let reader_recv_ = reader_recv.clone();
    //         let consensus_result_sender_ = consensus_result_sender.clone();
    //         thread_scope.spawn(move || {
    //             let mut align_param = AbpoaParam::default();
    //             align_param.match_score = 2;
    //             align_param.mismatch_score = 5;
    //             align_param.gap_open1 = 2;
    //             align_param.gap_open2 = 24;
    //             align_param.gap_ext1 = 1;
    //             align_param.gap_ext2 = 0;
    //             align_param.out_msa = false;
    //             align_param.out_consensus = true;
    //             align_param.set_align_mode(rsabpoa::abpoa::AlignMode::LOCAL);
    //             consensus_worker(
    //                 &align_param,
    //                 reader_recv_,
    //                 consensus_result_sender_,
    //                 thread_idx,
    //             );
    //         });
    //     }

    //     consensus_writer(consensus_result_recv, oup_bam);
    // });
}

#[cfg(test)]
mod test {
    use crate::{filter_and_sort_subreads, Subread};

    #[test]
    fn test_filter_and_sort_subreads() {
        let subreads = vec![
            Subread::new("ACTTTG".to_string(), "name".to_string(), 0, 3),
            Subread::new("ACTTTTTTG".to_string(), "name".to_string(), 0, 3),
            Subread::new(
                "ACTTTGGGGGACTTTGGGGGACTTTGGGGGACTTTGGGGG".to_string(),
                "name".to_string(),
                0,
                3,
            ),
            Subread::new("CCCCAAAGT".to_string(), "name".to_string(), 0, 3),
            Subread::new("CA".to_string(), "name".to_string(), 0, 3),
        ];

        let res = filter_and_sort_subreads(subreads.into());
        println!("{:?}", res);
    }
}
