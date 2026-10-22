// Rust implementations of umi.

#[derive(Clone, Debug, Default)]
pub struct UmiDigest {
    pub records: i64,
    pub hash: u64,
}

#[derive(Clone, Debug, Default)]
pub struct UmiKeep {
    pub bitmap: Vec<u64>,
    pub fragments: i64,
    pub molecules: i64,
    pub digest: UmiDigest,
}

fn umi_check(pattern: &str) {
    if pattern.is_empty() || !pattern.bytes().all(|c| c == b'N' || c == b'X') || !pattern.contains('N') {
        panic!("extractUmi: pattern '{}' must be N and X only, with at least one N", pattern);
    }
}

// The UMI and the read without it.
fn umi_cut(pattern: &[u8], e: &BioFastqEntry) -> (String, BioFastqEntry) {
    let (b, q) = (&e.bases.bases, &e.bases.quals);
    if b.len() < pattern.len() {
        panic!("extractUmi: read '{}' is shorter than the pattern", e.header);
    }
    let mut umi = String::with_capacity(pattern.len());
    let (mut kb, mut kq) = (Vec::with_capacity(b.len()), Vec::with_capacity(b.len()));
    for (i, &p) in pattern.iter().enumerate() {
        if p == b'N' {
            umi.push(b[i] as char);
        } else {
            kb.push(b[i]);
            kq.push(q[i]);
        }
    }
    kb.extend_from_slice(&b[pattern.len()..]);
    kq.extend_from_slice(&q[pattern.len()..]);
    (umi, BioSequence { header: e.header.clone(), bases: BioBaseQual { bases: kb, quals: kq } })
}

fn umi_tag(header: &str, umi: &str) -> String {
    match header.find(|c: char| c == ' ' || c == '\t') {
        Some(i) => format!("{}_{}{}", &header[..i], umi, &header[i..]),
        None => format!("{}_{}", header, umi),
    }
}

pub fn umi_extract_pairs(pattern: &String, xs: &Vec<(BioFastqEntry, BioFastqEntry)>) -> Vec<(BioFastqEntry, BioFastqEntry)> {
    umi_check(pattern);
    xs.iter()
        .map(|(r1, r2)| {
            let (umi, mut a) = umi_cut(pattern.as_bytes(), r1);
            a.header = umi_tag(&a.header, &umi);
            let mut b = r2.clone();
            b.header = umi_tag(&b.header, &umi);
            (a, b)
        })
        .collect()
}

pub fn umi_extract_reads(pattern: &String, xs: &Vec<BioFastqEntry>) -> Vec<BioFastqEntry> {
    umi_check(pattern);
    xs.iter()
        .map(|r| {
            let (umi, mut a) = umi_cut(pattern.as_bytes(), r);
            a.header = umi_tag(&a.header, &umi);
            a
        })
        .collect()
}

pub fn umi_of_header(h: &String) -> String {
    let first = h.split(|c: char| c == ' ' || c == '\t').next().unwrap_or("");
    match first.rfind('_') {
        Some(i) => first[i + 1..].to_string(),
        None => String::new(),
    }
}

fn umi_one_apart(a: &str, b: &str) -> bool {
    a.len() == b.len() && a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count() == 1
}

fn umi_clusters_idx(us: &[(&str, i64)]) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..us.len()).collect();
    order.sort_by(|&i, &j| us[j].1.cmp(&us[i].1).then(us[i].0.cmp(us[j].0)));
    let adj: Vec<Vec<usize>> = (0..us.len())
        .map(|i| (0..us.len()).filter(|&j| j != i && umi_one_apart(us[i].0, us[j].0) && us[i].1 >= 2 * us[j].1 - 1).collect())
        .collect();
    let mut found = vec![false; us.len()];
    let mut out = Vec::new();
    for &n in &order {
        if found[n] {
            continue;
        }
        let mut comp = vec![n];
        let mut seen = vec![false; us.len()];
        seen[n] = true;
        let mut k = 0;
        while k < comp.len() {
            for &m in &adj[comp[k]] {
                if !seen[m] {
                    seen[m] = true;
                    comp.push(m);
                }
            }
            k += 1;
        }
        for &m in &comp {
            found[m] = true;
        }
        out.push(comp);
    }
    out
}

pub fn umi_clusters(xs: &Vec<(String, i64)>) -> Vec<Vec<String>> {
    let us: Vec<(&str, i64)> = xs.iter().map(|(u, c)| (u.as_str(), *c)).collect();
    umi_clusters_idx(&us).into_iter().map(|c| c.into_iter().map(|i| xs[i].0.clone()).collect()).collect()
}

fn umi_pack(u: &str) -> u64 {
    if u.len() > 21 {
        panic!("UMI '{}' is longer than 21 bases", u);
    }
    let mut x = 0u64;
    for c in u.bytes() {
        let code = match c {
            b'A' => 1,
            b'C' => 2,
            b'G' => 3,
            b'T' => 4,
            b'N' => 5,
            _ => panic!("UMI '{}' holds '{}', not A, C, G, T or N", u, c as char),
        };
        x = x << 3 | code;
    }
    x
}

fn umi_unpack(mut x: u64) -> String {
    let mut out = Vec::new();
    while x != 0 {
        out.push(b"?ACGTN??"[(x & 7) as usize]);
        x >>= 3;
    }
    out.reverse();
    String::from_utf8(out).unwrap()
}

pub fn umi_records(_anchor: i64, keys: &DupKeys, umis: &Vec<String>, ords: &Vec<u64>) -> Vec<u64> {
    if keys.keyed.len() != umis.len() || umis.len() != ords.len() {
        panic!("umiRecords: {} keys, {} UMIs and {} indices", keys.keyed.len(), umis.len(), ords.len());
    }
    let mut out = Vec::with_capacity(4 * keys.words.len() / 2);
    let mut w = 0;
    for (i, &k) in keys.keyed.iter().enumerate() {
        if k {
            out.extend_from_slice(&[keys.words[w], keys.words[w + 1], umi_pack(&umis[i]), ords[i]]);
            w += 2;
        }
    }
    out
}

pub fn umi_stream_indices(at: u64, index: &Vec<u32>) -> Vec<u64> {
    index.iter().map(|&i| at + i as u64).collect()
}

fn umi_mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

fn umi_record_hash(r: &[u64]) -> u64 {
    r.iter().fold(0x9e3779b97f4a7c15u64, |h, &x| umi_mix(h ^ x))
}

pub fn umi_digest(xs: &Vec<u64>) -> UmiDigest {
    if xs.len() % 4 != 0 {
        panic!("umiDigest: {} words are not whole records", xs.len());
    }
    let mut d = UmiDigest::default();
    for r in xs.chunks_exact(4) {
        d.records += 1;
        d.hash = d.hash.wrapping_add(umi_record_hash(r));
    }
    d
}

pub fn umi_merge_digest(a: &UmiDigest, b: &UmiDigest) -> UmiDigest {
    UmiDigest { records: a.records + b.records, hash: a.hash.wrapping_add(b.hash) }
}

pub fn umi_same_digest(a: &UmiDigest, b: &UmiDigest) -> bool {
    a.records == b.records && a.hash == b.hash
}

// One position's UMIs: (packed UMI, fragments, where its fragments' stream
// indices start in `ords`, ascending). Marks every fragment of a molecule
// but the first in the stream as a duplicate.
fn umi_close_group(group: &mut Vec<(u64, i64, usize)>, ords: &mut Vec<u64>, keep: &mut UmiKeep) {
    if group.is_empty() {
        return;
    }
    let names: Vec<String> = group.iter().map(|g| umi_unpack(g.0)).collect();
    let us: Vec<(&str, i64)> = names.iter().map(|n| n.as_str()).zip(group.iter().map(|g| g.1)).collect();
    let mut assigned = vec![false; group.len()];
    for comp in umi_clusters_idx(&us) {
        let own: Vec<usize> = comp.into_iter().filter(|&k| !assigned[k]).collect();
        let first = own.iter().map(|&k| ords[group[k].2]).min().unwrap();
        for &k in &own {
            assigned[k] = true;
            let (_, n, at) = group[k];
            for &o in &ords[at..at + n as usize] {
                let w = (o / 64) as usize;
                if w >= keep.bitmap.len() {
                    panic!("umiKeepSorted: stream index {} is beyond the stream's fragments", o);
                }
                if o != first {
                    keep.bitmap[w] |= 1 << (o % 64);
                }
            }
        }
        keep.molecules += 1;
    }
    group.clear();
    ords.clear();
}

pub fn umi_keep_sorted<P>(n: i64, pull: P) -> UmiKeep
where
    P: Fn() -> (bool, String, Vec<u64>),
{
    let n = n.max(0) as usize;
    let mut keep = UmiKeep { bitmap: vec![0; (n + 63) / 64], ..Default::default() };
    let mut key: Option<(u64, u64)> = None;
    let mut group: Vec<(u64, i64, usize)> = Vec::new();
    let mut ords: Vec<u64> = Vec::new();
    let mut prev: Option<[u64; 4]> = None;
    loop {
        let (ok, msg, xs) = pull();
        if !ok {
            panic!("umiKeepSorted: {}", msg);
        }
        if xs.is_empty() {
            break;
        }
        if xs.len() % 4 != 0 {
            panic!("umiKeepSorted: {} words are not whole records", xs.len());
        }
        for r in xs.chunks_exact(4) {
            let r4 = [r[0], r[1], r[2], r[3]];
            if prev.map_or(false, |p| p > r4) {
                panic!("umiKeepSorted: records are not ascending");
            }
            prev = Some(r4);
            keep.fragments += 1;
            keep.digest.records += 1;
            keep.digest.hash = keep.digest.hash.wrapping_add(umi_record_hash(r));
            if key != Some((r[0], r[1])) {
                umi_close_group(&mut group, &mut ords, &mut keep);
                key = Some((r[0], r[1]));
            }
            match group.last_mut() {
                Some(g) if g.0 == r[2] => g.1 += 1,
                _ => group.push((r[2], 1, ords.len())),
            }
            ords.push(r[3]);
        }
    }
    umi_close_group(&mut group, &mut ords, &mut keep);
    keep
}

// A fragment is kept unless it duplicates an earlier fragment of its
// molecule; fragments without a key have no molecule and are all kept.
pub fn umi_keep_fragments(t: &UmiKeep, ords: &Vec<u64>) -> Vec<bool> {
    ords.iter()
        .map(|&i| match t.bitmap.get((i / 64) as usize) {
            Some(w) => w >> (i % 64) & 1 == 0,
            None => panic!("keepFragments: stream index {} is beyond the stream's fragments", i),
        })
        .collect()
}

pub fn umi_dedup_text(sample: &String, k: &UmiKeep) -> String {
    let (n, m) = (k.fragments, k.molecules);
    let pct = if n > 0 { 100.0 * (n - m) as f64 / n as f64 } else { 0.0 };
    format!("#id: 'umi_dedup'\n#section_name: 'UMI deduplication'\n#description: 'Fragments with a genome position, molecules after directional UMI deduplication, and the percent removed.'\n#plot_type: 'table'\nSample\tFragments\tMolecules\tPercent duplicates\n{}\t{}\t{}\t{}\n", sample, n, m, pct)
}
