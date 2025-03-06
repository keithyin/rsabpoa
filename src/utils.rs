pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    // DNA 碱基的互补映射表
    let complement = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N', // 未知碱基默认用 'N'
    };

    // 反向遍历并互补
    dna.iter().rev().map(|&b| complement(b)).collect()
}
