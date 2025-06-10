use chrono::Datelike;
use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use locus::{load_loci, Locus};
use rust_htslib::bam::IndexedReader;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::thread::{self, JoinHandle};
use workflow::run_workflow;

mod extend;
mod locus;
mod models;
mod profile;
mod workflow;

#[derive(Parser)]
#[command(name="HIFI-VCLUST",
          about="HiFi Variation Cluster Analysis Tool", 
          long_about = None,
          after_help = format!("Copyright (C) 2004-{} Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
#[command(arg_required_else_help(true))]
pub struct CliParams {
    #[clap(required = true)]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(long = "reads")]
    #[clap(help = "BAM file with aligned HiFi reads")]
    #[clap(value_name = "READS")]
    #[clap(num_args = 1)]
    #[arg(value_parser = check_file_exists)]
    pub reads_paths: PathBuf,

    #[clap(required = true)]
    #[clap(long = "regions")]
    #[clap(help = "BED file with region coordinates")]
    #[clap(value_name = "REGIONS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,

    #[clap(long = "threads")]
    #[clap(help = "Number of threads to use")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value_t = 1)]
    pub threads: usize,
}

type InputType = Option<Locus>;
type OutputType = Option<String>;

// Return some kind of Result/Status or something.
fn task_thread(
    reads_paths: Vec<PathBuf>,
    task_receiver: Receiver<InputType>,
    result_sender: Sender<OutputType>,
) -> Result<(), String> {
    let mut bams = Vec::new();
    for path in reads_paths {
        let bam = IndexedReader::from_path(&path).map_err(|e| e.to_string())?;
        bams.push(bam);
    }
    loop {
        match task_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(locus)) => match run_workflow(&mut bams, &locus) {
                Err(message) => {
                    log::warn!("{message}");
                }
                Ok(result) => {
                    result_sender.send(Some(result)).unwrap();
                }
            },
        }
    }

    result_sender.send(None).unwrap();

    Ok(())
}

fn read_bam_paths(file_path: PathBuf) -> std::io::Result<Vec<PathBuf>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut paths = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if !line.trim().is_empty() {
            paths.push(PathBuf::from(line));
        }
    }

    Ok(paths)
}

fn main() -> Result<(), String> {
    let args = CliParams::parse();

    let paths = read_bam_paths(args.reads_paths).map_err(|e| e.to_string())?;
    // Create channels for communication between threads
    let (task_sender, task_receiver): (Sender<InputType>, Receiver<InputType>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    let task_handles: Vec<JoinHandle<Result<(), String>>> = (0..args.threads)
        .map(|_| {
            let m_reads = paths.clone();
            let m_receiver = task_receiver.clone();
            let m_result_sender = result_sender.clone();

            thread::spawn(move || task_thread(m_reads, m_receiver, m_result_sender))
        })
        .collect();

    // Push each of the loci to the channel
    let loci = load_loci(args.repeats_path)?;
    for locus in loci {
        task_sender.send(Some(locus)).unwrap();
    }

    // Signal worker threads to exit
    for _ in 0..args.threads {
        task_sender.send(None).unwrap();
    }

    // Collect results
    let mut n_done = 0;
    while n_done < args.threads {
        match result_receiver.recv() {
            Ok(None) | Err(_) => {
                n_done += 1;
            }
            Ok(Some(result)) => {
                println!("{result}");
            }
        }
    }

    // Close up
    for handle in task_handles {
        let _ = handle.join().unwrap();
    }

    // For now, we'll just have the task_handles hold the lines
    Ok(())
}

fn check_file_exists(path: &str) -> Result<PathBuf, String> {
    let path = Path::new(path);
    if path.exists() {
        Ok(path.to_path_buf())
    } else {
        Err(format!("File does not exist: {}", path.display()))
    }
}
