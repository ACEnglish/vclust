use crate::locus::Locus;
use crate::models::{MODEL_REF, MODEL_VC, PRIOR_REF, PRIOR_VC, RADIUS};
use crate::profile::{get_profile, Prof};
use itertools::Itertools;
use logaddexp::LogAddExp;
use rust_htslib::bam::IndexedReader;

pub fn get_extension_offsets(
    locus: &Locus,
    bams: &mut Vec<IndexedReader>,
) -> Option<(i64, i64, i64)> {
    let region = extend_region(locus).ok()?;

    let mut ns = 0;
    // let alt_minimum = 0.35;
    let alt_depth = 5.0;
    // Add as you go
    let mut sum_alts: Option<Vec<f64>> = None;
    let mut sum_depth: f64 = 0.0;
    let mut count: usize = 0;

    for bam in bams {
        let (prof, any_alt) = get_profile(bam, region).ok()?;

        if let Some(ref mut alts) = sum_alts {
            for (sum, alt) in alts.iter_mut().zip(prof.alts.iter()) {
                *sum += alt;
                //any_alt |= *alt >= alt_minimum;
            }
        } else {
            //any_alt |= prof.alts.iter().any(|&alt| alt >= alt_minimum);
            sum_alts = Some(prof.alts.clone());
        }
        sum_depth += prof.depth;
        count += 1;
        if any_alt & (prof.depth >= alt_depth) {
            ns += 1;
        }
    }

    let prof = if let Some(sum_alts) = sum_alts {
        let alts = sum_alts.into_iter().map(|sum| sum / count as f64).collect();
        let depth = sum_depth / count as f64;
        Prof { alts, depth }
    } else {
        return None;
    };

    if prof.depth < 5.0 || prof.depth > 150.0 {
        return None;
    }

    let alts = discretize(&prof.alts);

    let span = (RADIUS, RADIUS + locus.end - locus.start);
    let span = extend_to_ref_flanks(&alts, span, 150)?;
    let span = extend_to_ref_flanks(&alts, span, 50)?;
    let span = extend_to_ref_flanks(&alts, span, 25)?;
    let span = extend_to_ref_flanks(&alts, span, 10)?;

    let lf_offset = RADIUS - span.0;
    let rf_offset = span.1 - (RADIUS + locus.end - locus.start);

    Some((lf_offset, rf_offset, ns))
}

fn extend_region(locus: &Locus) -> Result<(&str, i64, i64), String> {
    if locus.start < RADIUS {
        Err("Locus too close to chromosome start".to_string())
    } else {
        Ok((&locus.chrom[..], locus.start - RADIUS, locus.end + RADIUS))
    }
}

fn discretize(vals: &[f64]) -> Vec<u8> {
    vals.iter()
        .map(|val| {
            if *val <= 0.10 {
                0
            } else if *val < 0.25 {
                1
            } else if *val < 0.75 {
                2
            } else if *val < 1.50 {
                3
            } else if *val < 5.00 {
                4
            } else {
                5
            }
        })
        .collect()
}

fn extend_to_ref_flanks(alts: &[u8], span: (i64, i64), window_len: i64) -> Option<(i64, i64)> {
    let mut lf_pos = span.0 - window_len;
    while lf_pos >= 0 {
        let window = &alts[lf_pos as usize..(lf_pos + window_len) as usize];
        let window = window.iter().rev().copied().collect_vec();
        let prob_ref = assess_window(&window[..]);
        if prob_ref >= 0.5 {
            break;
        }
        lf_pos -= 1;
    }

    if lf_pos == 0 {
        return None;
    }

    let mut rf_pos = span.1;
    while rf_pos <= alts.len() as i64 - window_len {
        let window = &alts[rf_pos as usize..(rf_pos + window_len) as usize];
        let prob_ref = assess_window(window);

        if prob_ref >= 0.5 {
            break;
        }
        rf_pos += 1;
    }

    if alts.len() as i64 - window_len < rf_pos {
        return None;
    }

    Some((lf_pos + window_len, rf_pos))
}

fn assess_window(vals: &[u8]) -> f64 {
    let ll_norm = get_loglik(vals, &MODEL_REF) + PRIOR_REF.ln();
    let ll_poly = get_loglik(vals, &MODEL_VC) + PRIOR_VC.ln();
    let ll_sum = ll_norm.ln_add_exp(ll_poly);

    (ll_norm - ll_sum).exp()
}

fn get_loglik(prof: &[u8], model: &[f64; 1500]) -> f64 {
    let mut ll = 0.0;
    for (pos, val) in prof.iter().enumerate() {
        ll += model[pos * 6 + *val as usize].ln();
    }
    ll
}
