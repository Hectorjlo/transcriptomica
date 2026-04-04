// Use module "filesystem", for manipulation of files
// Use struct "Path" of the module path, for inspecting system paths
use std::fs;
use std::path::Path;


// It also uses the method "new" of the struct "Regex" in the module "regex"

// Returns four owned String types
// Takes one argument: A reference to a string slice
// Parses the summary of hisat2
fn parse_hisat2_summary(content: &str) -> (String, String, String, String) {
    // Searchs pattern for any number until the percentaje before "overall alignment rate"
    let regex_overall_alignment_rate = regex::Regex::new(r"([\d.]+)% overall alignment rate").unwrap();
    
    // For paired-end summary files
    // The regex pattern will look for any number until a "%" before <String> that matches 
    let regex_aligned_zero_pe = regex::Regex::new(r"([\d.]+)%\) aligned concordantly 0 times").unwrap();
    let regex_aligned_once_pe = regex::Regex::new(r"([\d.]+)%\) aligned concordantly exactly 1 time").unwrap();
    let regex_multi_aligned_pe = regex::Regex::new(r"([\d.]+)%\) aligned concordantly >1 times").unwrap();

    // For single-end
    // Same principle as single-end
    let regex_align_zero = regex::Regex::new(r"([\d.]+)%\) aligned 0 times").unwrap();
    let regex_align_once = regex::Regex::new(r"([\d.]+)%\) aligned exactly 1 time").unwrap();
    let regex_align_multi = regex::Regex::new(r"([\d.]+)%\) aligned >1 times").unwrap();

    // Creating the own types with the number gathered from previous expressions

    // For the overal aligment
    let overall: String = regex_overall_alignment_rate // Pass the pattern 
        .captures(content) // Search the pattern in the content
        .map(|c| c[1].to_string()) // Take the first group "([\d.]+)" of the match
        .unwrap_or("N/A".into()); // Get the number or if not found return N/A as String

    let zero: String = regex_aligned_zero_pe
        .captures(content)
        .or_else(|| regex_align_zero.captures(content))
        .map(|c| c[1].to_string())
        .unwrap_or("N/A".into());

    let once: String = regex_aligned_once_pe
        .captures(content)
        .or_else(|| regex_align_once.captures(content))
        .map(|c| c[1].to_string())
        .unwrap_or("N/A".into());

    let multi: String = regex_multi_aligned_pe
        .captures(content)
        .or_else(|| regex_align_multi.captures(content))
        .map(|c| c[1].to_string())
        .unwrap_or("N/A".into());

    // Return the 4 String
    (overall, zero, once, multi)
}

// Returns four owned String types
// Takes one argument: A reference to a string slice
// Parses the summary of star
fn parse_star_summary(content: &str) -> (String, String, String, String) {
    // Regex pattern that will search any digits after <String> and group with parenthesis
    let regex_unique_mapped = regex::Regex::new(r"Uniquely mapped reads %\s*\|\s*([\d.]+)%").unwrap();
    let regex_multi_mapped = regex::Regex::new(r"% of reads mapped to multiple loci\s*\|\s*([\d.]+)%").unwrap();
    let regex_toomany_mapped = regex::Regex::new(r"% of reads mapped to too many loci\s*\|\s*([\d.]+)%").unwrap();

    // Gathering the numbers of the stats

    let unique: f64 = regex_unique_mapped // Use this pattern
        .captures(content) // Captures the match in the "content"
        .map(|c| c[1].parse().unwrap_or(0.0)) // Map to an Option<f64> 
        .unwrap_or(0.0); // Return the number if found, or return 0.0 if not found

    let multi: f64 = regex_multi_mapped
        .captures(content)
        .map(|c| c[1].parse().unwrap_or(0.0))
        .unwrap_or(0.0);

    let toomany: f64 = regex_toomany_mapped
        .captures(content)
        .map(|c| c[1].parse().unwrap_or(0.0))
        .unwrap_or(0.0);

    // Add up all the aligns
    let overall = unique + multi + toomany;
    // Calculate the non aligned 
    let non_aligned = 100.0 - overall;

    // Convert to String type with 2 decimals using the macro format!
    (
        format!("{:.2}", overall),
        format!("{:.2}", non_aligned),
        format!("{:.2}", unique),
        format!("{:.2}", multi + toomany),
    )
}

//Creating the struct that manages the output of the "fasta-like" 
pub struct Entry {
    pub assembler: String,
    pub srr: String,
    pub mode: String,
    pub overall: String,
    pub non_aligned: String,
    pub once: String,
    pub multi: String,
}

// Collect the files of hisat2 summaries
// Returns: A vector of Entry structs
pub fn collect_hisat2(results_dir: &Path) -> Vec<Entry> {
    // Arguments:
    //      (results_dir: &Path): A reference to a Path struct

    // Create the Regex pattern, to search for the files
    let regex_pattern_file = regex::Regex::new(r"summary_of_(SRR\d+)_(pe|1)\.txt").unwrap();
    // Initialize the Vector of Entry
    let mut entries: Vec<Entry> = Vec::new();

    // Search for the subdirectories in hisat2 
    for subdir in &["paired_end", "single_end"] {
        let dir = results_dir.join("hisat2").join(subdir);
        // If the directory does not exists search for the next one
        if !dir.exists() {
            continue;
        }

        // Create the files variable, that contains a Vector of PathBuf
        let mut files: Vec<_> = fs::read_dir(&dir) // Reads the directory
            .unwrap() // Unwraps the result
            .filter_map(|e| e.ok()) // Checks if the Result enum is ok and discards the None variants
            .map(|e| e.path()) // Gathers the Path of the remaining files
            .filter(|p| { // Filter to Keep only Paths 
                p.file_name() // Gathers the file name
                    .and_then(|n| n.to_str()) // Calls .to_str() if Some(T) case
                    .map(|n| n.starts_with("summary_of_")) // Maps the Option<&str> to Option<bool> 
                    .unwrap_or(false) // Keep those that passed the filter
            })
            .collect(); // Create the Vec<PathBuf>
        files.sort(); // Sort the files

        // For each file path in the Vector of File paths
        for fpath in files {
            // Gathers the file name given the file path
            let fname = fpath.file_name()
                .unwrap()
                .to_str()
                .unwrap();

            // If the regex capture is the Some(caps) variant
            if let Some(caps) = regex_pattern_file.captures(fname) {
                // Selects the first group of the regex_pattern_file -> "SRR11111"
                let srr = caps[1].to_string();
                // Selects the second group of the regex_pattern_file -> "pe|1"
                // and then uses "PE" or "SE"
                let mode = if &caps[2] == "pe" { "PE" } else { "SE" };
                // Reads the file to a string
                let content = fs::read_to_string(&fpath).unwrap();
                // Call to the function
                let (overall, non_aligned, once, multi) = parse_hisat2_summary(&content);
                // Pushing the Entry structs for the HISAT2
                entries.push(Entry {
                    assembler: "HISAT2".into(), // Assembler name
                    srr, // Pass the srr match
                    mode: mode.into(), // Pass the mode into a String
                    overall, // Pass the stats gather from the function
                    non_aligned,
                    once,
                    multi,
                });
            }
        }
    }
    entries
}

// Same priciple, but with the expression of Star summary files
pub fn collect_star(results_dir: &Path) -> Vec<Entry> {
    // Pattern to search for the file names
    let re_file = regex::Regex::new(r"output_(SRR\d+)_(pe|1)_Log\.final\.out").unwrap();
    // Initialize the Vector of Entry
    let mut entries = Vec::new();

    // Search in the subdirectory
    for subdir in &["paired_end", "single_end"] {
        // Creates the sudirectory path, given the results directory
        let dir = results_dir.join("star").join(subdir);
        
        // If those subdirectories are not found pass to next
        if !dir.exists() {
            continue;
        }
        // Creates the variable files
        // It will contain a Vector of PathBuf
        let mut files: Vec<_> = fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .map(|n| n.ends_with("_Log.final.out"))
                    .unwrap_or(false)
            })
            .collect();
        files.sort();

        // Process each file to gather the stats
        for fpath in files {
            let fname = fpath.file_name().unwrap().to_str().unwrap();
            if let Some(caps) = re_file.captures(fname) {
                let srr = caps[1].to_string();
                let mode = if &caps[2] == "pe" { "PE" } else { "SE" };
                let content = fs::read_to_string(&fpath).unwrap();
                let (overall, non_aligned, once, multi) = parse_star_summary(&content);
                // Push the Entry struct
                entries.push(Entry {
                    assembler: "STAR".into(),
                    srr,
                    mode: mode.into(),
                    overall,
                    non_aligned,
                    once,
                    multi,
                });
            }
        }
    }
    // Return the entries vec
    entries
}

// Struct to manage parse of arguments
pub struct Config {
    pub results_path: String,
}

// Implementation block for the Config struct
impl Config {
    // Function to build a Config struct
    pub fn build(
        // Take mutability of the argument, acepts an data type that implements Iterator with elements of String
        mut args: impl Iterator<Item = String>,
        // Returns a Result enum, Ok(Config) or Err(&'static str)
    ) -> Result<Config, &'static str> {
        // Pass the first element that is path of the executable
        args.next();

        let results_path = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a directory path")
        };

        Ok(Config { results_path })
    }
}