use dynamic_weighted_index::DynamicWeightedIndex;
use rand::Rng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use pcg_rand::Pcg64;
use rand_distr::{StandardNormal, Distribution};
use std::{
    hint::black_box,
    time::Instant,
};
use std::path::Path;

fn setup(rng: &mut Pcg64, size: usize) -> DynamicWeightedIndex<f64> {
    let mut dw_index = DynamicWeightedIndex::new(size);
    for i in 0..size {
        let weight: f64 = StandardNormal.sample(rng);
        let abs_weight = weight.abs();
        dw_index.set_weight(i, abs_weight);
    }
    dw_index
}

fn sample_fixed(
    rng: &mut Pcg64,
    dw_index: &DynamicWeightedIndex<f64>,
    n: usize,
) -> Vec<usize> {
    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        samples.push(dw_index.sample(rng).unwrap());
    }
    samples
}

fn sample_variable(
    rng: &mut Pcg64,
    dw_index: &mut DynamicWeightedIndex<f64>,
    n: usize,
) -> Vec<usize> {
    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        samples.push(dw_index.sample(rng).unwrap());
        let index_to_modify = rng.gen_range(0..n);
        dw_index.remove_weight(index_to_modify);
        let new_weight: f64 = StandardNormal.sample(rng);
        let abs_new_weight = new_weight.abs();
        dw_index.set_weight(index_to_modify, abs_new_weight);
    }
    samples
}

fn sample_decreasing(
    rng: &mut Pcg64,
    dw_index: &mut DynamicWeightedIndex<f64>,
    removal_order: &[usize],
) -> Vec<usize> {
    let steps = removal_order.len() - removal_order.len() / 10;
    let mut samples = Vec::with_capacity(steps);
    for &index_to_remove in removal_order.iter().take(steps) {
        samples.push(dw_index.sample(rng).unwrap());
        dw_index.remove_weight(index_to_remove);
    }
    samples
}


fn numerical_check(rng: &mut Pcg64) -> f64 {
    // Initialize a DynamicWeightedIndex with 3 slots
    let mut dw_index = DynamicWeightedIndex::new(3);

    // Set the initial weights: [0.1, 0.9, 100000.0]
    dw_index.set_weight(0, 0.1);
    dw_index.set_weight(1, 0.9);
    dw_index.set_weight(2, 9.0001*1.0e15);

    // “Update” index 2 to weight = 0.0
    dw_index.set_weight(2, 0.0);

    // Draw 1_000 samples and count how often we get index 0
    let mut count_zero = 0;
    for _ in 0..1_000_00 {
            if dw_index.sample(&mut rand::thread_rng()) == Some(0) {
                count_zero += 1;
            }
    }
     println!("{:}", count_zero);

    // Return the empirical probability of drawing 0
    (count_zero as f64) / 1_000_00.0
}

use std::fs::File;
use std::io::{self, BufWriter, Write};

fn empirical_probabilities(
    rng: &mut Pcg64,
    n: usize,
    t: usize,
) -> Vec<f64> {
    let mut dw_index = DynamicWeightedIndex::new(n);

    let mut weights: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let i1 = (i + 1) as f64;
        let base_i = 2.0 + (1.0 / ((100 * n) as f64)) * i1;
        let w0 = base_i.powf(1000.0);
        weights.push(w0);
    }

    for i in 0..n {
        dw_index.set_weight(i, weights[i]);
    }

    for _ in 0..t {
        for i in 0..n {
            let i1 = (i + 1) as f64;
            let base_i = 2.0 + (1.0 / ((100 * n) as f64)) * i1;
            weights[i] /= base_i;
            dw_index.set_weight(i, weights[i]);
        }
    }

    let num_samples = 1_000_000;
    let mut counts: Vec<usize> = vec![0; n];
    for _ in 0..num_samples {
        if let Some(idx) = dw_index.sample(rng) {
            counts[idx] += 1;
        }
    }

    counts.into_iter()
        .map(|c| (c as f64) / (num_samples as f64))
        .collect()
}

use std::error::Error;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};

/// Read `path` into memory, append `column_name` with one entry per row from `values`,
/// and atomically rewrite the original CSV.
fn append_column_to_csv(
    path: &Path,
    column_name: &str,
    values: &[f64],
) -> Result<(), Box<dyn Error>> {
    // 1) Read existing CSV
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;
    let headers = rdr.headers()?.clone();
    let mut records: Vec<StringRecord> = rdr
        .records()
        .enumerate()
        .map(|(i, r)| {
            let mut rec = r?;
            if i >= values.len() {
                Err("CSV has more rows than values provided")?
            }
            Ok(rec)
        })
        .collect::<Result<_, Box<dyn Error>>>()?;

    // 2) Prepare a new temporary writer
    let tmp_path = path.with_extension("tmp");
    let mut wtr = WriterBuilder::new()
        .has_headers(true)
        .from_path(&tmp_path)?;

    // 3) Write out headers + new column name
    let mut new_headers = headers.clone();
    new_headers.push_field(column_name);
    wtr.write_record(&new_headers)?;

    // 4) Write each record with the appended value
    for (i, mut rec) in records.into_iter().enumerate() {
        let val = values
            .get(i)
            .ok_or("Not enough values to fill all rows")?
            .to_string();
        rec.push_field(&val);
        wtr.write_record(&rec)?;
    }
    wtr.flush()?;

    // 5) Atomically replace the original file
    std::fs::rename(tmp_path, path)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = Pcg64::seed_from_u64(42);

    let s: usize = 100;

    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));

    let numerical = manifest_dir
        .parent().unwrap()
        .join("data")
        .join("forest_of_trees_numerical.csv");

    let file = File::create(numerical)
        .expect("Failed to create probabilities_mc.csv");
    let mut writer = BufWriter::new(file);

    for t in 1..=52 {
        // Compute empirical probabilities for this t
        let probs: Vec<f64> = empirical_probabilities(&mut rng, s, t);

        // Convert each f64 to string and write as one CSV row
        let row: Vec<String> = probs.iter()
            .map(|p| p.to_string())
            .collect();
        writeln!(writer, "{}", row.join(","))
            .expect("Failed to write CSV row");
        writer.flush();
        print!("{}", t);
        io::stdout().flush().unwrap();
    }

    writer.flush().expect("Failed to flush writer");

    let sample_sizes: Vec<usize> = (3..=7).map(|i| 10usize.pow(i)).collect();
    let mut repetitions = 50;

    // --- 1) static (fixed) sampling ---
    let mut static_medians = Vec::with_capacity(sample_sizes.len());
    for &size in &sample_sizes {
        if size > 100000 {repetitions = 5;}
        let dw_index = setup(&mut rng, size);
        let mut times: Vec<u128> = Vec::with_capacity(repetitions);
        for _ in 0..repetitions {
            let idx = black_box(&dw_index);
            let n   = black_box(size);
            let start = Instant::now();
            sample_fixed(&mut rng, idx, n);
            times.push(start.elapsed().as_nanos());
        }
        times.sort();
        let median_ns = if repetitions % 2 == 0 {
            (times[repetitions / 2 - 1] + times[repetitions / 2]) / 2
        } else {
            times[repetitions / 2]
        };
        let median_per_sample = median_ns as f64 / size as f64;
        println!(
            "Size: {:>8}, Fixed ({} reps): {:.2} ns/sample (median)",
            size, repetitions, median_per_sample
        );
        static_medians.push(median_per_sample);
    }

    // write static_medians out
    let static_csv = manifest_dir
        .parent().unwrap()
        .join("data")
        .join("static.csv");
    append_column_to_csv(&static_csv, "FOREST_OF_TREES", &static_medians)
        .expect("Failed to append FOREST_OF_TREES to static.csv");

    repetitions = 50;

    // --- 2) increasing-range sampling ---
    let mut variable_medians = Vec::with_capacity(sample_sizes.len());
    for &size in &sample_sizes {

        if size > 100000 {repetitions = 5;}
        let mut times: Vec<u128> = Vec::with_capacity(repetitions);
        for _ in 0..repetitions {
            let mut dw_index = setup(&mut rng, size);
            let idx = black_box(&mut dw_index);
            let n   = black_box(size);
            let start = Instant::now();
            sample_variable(&mut rng, idx, n);
            times.push(start.elapsed().as_nanos());
        }
        times.sort();
        let median_ns = if repetitions % 2 == 0 {
            (times[repetitions / 2 - 1] + times[repetitions / 2]) / 2
        } else {
            times[repetitions / 2]
        };
        let median_per_sample = median_ns as f64 / size as f64;
        println!(
            "Size: {:>8}, Variable ({} reps): {:.2} ns/sample (median)",
            size, repetitions, median_per_sample
        );
        variable_medians.push(median_per_sample);
    }

    // write variable_medians out
    let variable_csv = manifest_dir
        .parent().unwrap()
        .join("data")
        .join("dynamic_fixed.csv");
    append_column_to_csv(&variable_csv, "FOREST_OF_TREES", &variable_medians)
        .expect("Failed to append FOREST_OF_TREES to dynamic_fixed.csv");

    repetitions = 50;

    // --- 3) decreasing-range sampling ---
    let mut decreasing_medians = Vec::with_capacity(sample_sizes.len());
    for &size in &sample_sizes {
        if size > 100000 {repetitions = 5;}
        let mut removal_order: Vec<usize> = (0..size).collect();
        removal_order.shuffle(&mut rng);
        let steps = size - size / 10;
        let mut times: Vec<u128> = Vec::with_capacity(repetitions);
        for _ in 0..repetitions {
            let mut dw_index = setup(&mut rng, size);
            let idx = black_box(&mut dw_index);
            let order = black_box(&removal_order);
            let start = Instant::now();
            sample_decreasing(&mut rng, idx, order);
            times.push(start.elapsed().as_nanos());
        }
        times.sort();
        let median_ns = if repetitions % 2 == 0 {
            (times[repetitions / 2 - 1] + times[repetitions / 2]) / 2
        } else {
            times[repetitions / 2]
        };
        let median_per_sample = median_ns as f64 / steps as f64;
        println!(
            "Size: {:>8}, Decreasing ({} reps): {:.2} ns/sample (median)",
            size, repetitions, median_per_sample
        );
        decreasing_medians.push(median_per_sample);
    }

    let decreasing_csv = manifest_dir
        .parent().unwrap()
        .join("data")
        .join("dynamic_decreasing.csv");
    append_column_to_csv(&decreasing_csv, "FOREST_OF_TREES", &decreasing_medians)
        .expect("Failed to append FOREST_OF_TREES to dynamic_decreasing.csv");

    Ok(())
}

