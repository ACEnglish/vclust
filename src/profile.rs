use itertools::Itertools;
use rust_htslib::bam::{self, IndexedReader, Record};

pub type CigarOp = rust_htslib::bam::record::Cigar;

pub type Region<'a> = (&'a str, i64, i64);

#[derive(Debug, PartialEq, Clone)]
pub struct Cigar {
    pub ref_pos: i64,
    pub ops: Vec<CigarOp>,
}

#[derive(Debug)]
pub struct Prof {
    pub alts: Vec<f64>,
    pub depth: f64,
}

pub fn get_profile(bam: &mut IndexedReader, region: Region) -> Result<(Prof, bool), String> {
    let prof_len = (region.2 - region.1) as usize;
    let mut covs = vec![0; prof_len];
    let mut alts = vec![0; prof_len];
    bam.fetch(region).map_err(|e| e.to_string())?;
    let mut any_alt = 0;
    for (index, rec) in bam::Read::records(bam).enumerate() {
        let rec = rec.map_err(|e| e.to_string())?;

        if rec.is_secondary() || rec.is_supplementary() || rec.mapq() < 50 {
            continue;
        }
        any_alt += update_profs(rec, &mut covs, &mut alts, region) as usize;

        // Absolute max depth
        if index >= 200 {
            return Err("High depth".to_string());
        }
    }

    let depth = get_mean(&covs);

    let alts = alts
        .iter()
        .map(|v| *v as f64 / depth.max(1.0))
        .collect_vec();

    Ok((Prof { alts, depth }, any_alt >= 3))
}

pub fn update_profs(rec: Record, covs: &mut [u32], alts: &mut [u32], region: Region) -> bool {
    assert_eq!(covs.len() as i64, region.2 - region.1);
    assert_eq!(covs.len(), alts.len());

    let mut ref_pos = rec.pos();
    let region_start = region.1;
    let region_end = region.2;
    let mut any_alt = false;
    for op in rec.cigar().iter() {
        let op_len = get_ref_len(op);

        // Skip operations entirely before the region
        if ref_pos + op_len <= region_start {
            ref_pos += op_len;
            continue;
        }

        // Stop if we’re entirely past the region
        if ref_pos >= region_end {
            break;
        }

        // Determine overlap with the region
        let clipped_start = ref_pos.max(region_start);
        let clipped_end = (ref_pos + op_len).min(region_end);
        let clipped_len = (clipped_end - clipped_start) as usize;

        let index = (clipped_start - region_start) as usize;

        match op {
            CigarOp::Match(_) | CigarOp::Equal(_) => {
                let slice = &mut covs[index..index + clipped_len];
                for cov in slice.iter_mut() {
                    *cov += 1;
                }
            }
            CigarOp::Diff(_) | CigarOp::Del(_) => {
                let cov_slice = &mut covs[index..index + clipped_len];
                let alt_slice = &mut alts[index..index + clipped_len];
                for (cov, alt) in cov_slice.iter_mut().zip(alt_slice.iter_mut()) {
                    *cov += 1;
                    *alt += 1;
                }
                any_alt |= clipped_len >= 5;
            }
            CigarOp::Ins(len) => {
                // Insertions don't consume reference but we can still bump alt at insertion site
                if ref_pos >= region_start && ref_pos < region_end {
                    let idx = (ref_pos - region_start) as usize;
                    alts[idx] += *len;
                }
                any_alt |= clipped_len >= 5;
            }
            CigarOp::SoftClip(_) => {
                if ref_pos >= region_start && ref_pos < region_end {
                    let idx = (ref_pos - region_start) as usize;
                    alts[idx] += 1;
                }
            }
            CigarOp::HardClip(_) | CigarOp::Pad(_) | CigarOp::RefSkip(_) => {
                panic!("Unexpected operation {:?}", op);
            }
        }

        // Advance the reference position if this op consumes ref
        ref_pos += match op {
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
            _ => op_len,
        };
    }
    //Return
    any_alt
}

fn get_mean(vals: &[u32]) -> f64 {
    vals.iter().sum::<u32>() as f64 / vals.len() as f64
}

fn get_ref_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::RefSkip(len)
        | CigarOp::Del(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len) => *len as i64,
        CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}
