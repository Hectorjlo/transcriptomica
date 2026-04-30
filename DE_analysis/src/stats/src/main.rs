use std::path::Path;
use std::{env, process};
use stats::{collect_hisat2, collect_star, Entry, Config};

// Creates a normalized "fasta-like" stats for "HISAT2" and "STAR" 
// Usage:
//      ./stats <results_dir>
//      ./stats results/

fn main() {
    // Build the parser
    let config: Config = Config::build(env::args()).unwrap_or_else(|err| {
        eprint!("[ERROR] While parsing arguments {err}, missing directory");
        process::exit(1);
    });
    // Path of the results where hisat2 and star subfolder exists
    let results_dir = Path::new(&config.results_path);

    // Check if the results directory exists, else stop the execution
    if !results_dir.exists() {
        eprintln!("[ERROR] Directory not found");
        process::exit(1);
    }

    // Create the vector of all Entry structs
    let mut all: Vec<Entry> = Vec::new();

    // Call to a function to generate the vector of Entry of hisat2
    all.extend(collect_hisat2(results_dir));

    // Call to a function to generate the vector of Entry of star
    all.extend(collect_star(results_dir));

    // Print the Entry structs results
    for e in &all {
        println!(">ASSEMBLER: {} READ: {} ({})", e.assembler, e.srr, e.mode);
        println!("  Overall alignment: {}%", e.overall);
        println!("\tNon aligned: {}%", e.non_aligned);
        println!("\tAligned only once: {}%", e.once);
        println!("\tAligned >1 times: {}%", e.multi);
    }
}