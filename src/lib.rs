use anyhow::{bail, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};

use crate::utils::distance_functions;

pub fn kmer_similarity(
    base_seq: &[char],
    mod_sequence: &[char],
    distance_function: &str,
    k: usize,
    step: usize
) -> Result<f64> {
    let base_kmers = generate_kmers_with_step(base_seq, k, step)?;
    let mod_kmers = generate_kmers_with_step(mod_sequence, k, step)?;

    distance_functions::match_distance_function(
        base_kmers,
        mod_kmers,
        distance_function
    )
}

// This function simply takes a string and returns the set of k-mers.
pub fn generate_kmers(sequence: &[char], k: usize) -> Result<Vec<&[char]>> {
    generate_kmers_with_step(sequence, k, 1)
}

pub fn generate_kmers_with_step(
    sequence: &[char],
    k: usize,
    step: usize
) -> Result<Vec<&[char]>> {
    let mut kmers = vec![];

    if k > sequence.len() {
        bail!(
            "K ({}) is less than string length {}",
            k,
            sequence.len()
        );
    }

    for i in (0..=(sequence.len() - k)).step_by(step) {
        kmers.push(&sequence[i..i + k]);
    }
    Ok(kmers)
}

/* This function calculates the euclidean distance between kmer-vector representations
* of minimizer representations of strings. */
pub fn minimizer_similarity(
    base_seq: &[char],
    mod_seq: &[char],
    distance_function: &str,
    k: usize,
    w: usize,
    step: usize
) -> Result<f64> {
    let base_minimizers = generate_minimizers(base_seq, k, w, step)?;
    let mod_minimizers = generate_minimizers(mod_seq, k, w, step)?;

    distance_functions::match_distance_function(
        base_minimizers,
        mod_minimizers,
        distance_function
    )
}

/* This function generates a set of minimizers from a string. It uses my_hash_function
 * with a set seed (69). It only uses ONE hash function.
 */
pub fn generate_minimizers(
    seq: &[char],
    k: usize,
    w: usize,
    step: usize
) -> Result<Vec<&[char]>> {
    if k > w {bail!("k ({}) > w ({})", k, w);}
    if w > seq.len() {bail!("w ({}) > sequence length ({})", w, seq.len());}

    let mut minimizers: Vec<&[char]> = Vec::new();
    for idx in (0..seq.len() - w).step_by(step) {
        let window = &seq[idx..idx + w];
        let minimizer = generate_single_minimizer(window, k)?;
        minimizers.push(minimizer);
    }
    Ok(minimizers)
}

pub fn generate_single_minimizer(window: &[char], k: usize) -> Result<&[char]> {
    let seed: u64 = 69;
    let kmers = generate_kmers(window, k)?;
    let minimizer = kmers
            .iter()
            .min_by_key(|item| my_hash_function(item, seed))
            .unwrap();
    Ok(minimizer)
}

/* This is a simple hash function that maps items (kmers) to u64 integers.
* It uses a seed, so this is technically a series of hash functions.
*/
fn my_hash_function<T: Hash> (item: &T, seed: u64) -> u64 {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish()
}

/* This struct manages ownership of string slices. Strobemers are constructed from strobes,
 * each of which is its own borrowed string slice. These slices need to be concatenated
 * into a bigger string slice, so this struct takes ownership so the concatenated string
 * slice stays alive while it's being used for calculations.
 */
struct Strobemer {
    chars: Vec<char>
}

impl Strobemer {
    pub fn new(strobes: Vec<&[char]>) -> Self {
        let strobemer_length: usize = strobes.iter()
            .map(|strobe| strobe.len()).sum();
        let mut strobemer = Vec::with_capacity(strobemer_length);
        for strobe in strobes {
            strobemer.extend_from_slice(strobe);
        }
        Strobemer {chars: strobemer}
    }

    pub fn data_ref(&self) -> &[char] {
        &self.chars
    }
}

/*  This function generates a string's set of strobemers. Strobemers are quirky,
so here's an intuitive breakdown of what strobemers are and how they're generated:

* Stromebers are the concatenation of strobes. There are *order* strobes in each
strobemer.
* Each strobemer is drawn from a span of length strobemer_span_length =
            strobe_length + (strobe_window_gap + strobe_window_length) * order,
representing the first strobe and all of the remaining strobe windows.
* In any sequence, there are
            (sequence.len() - strobemer_span_length) / step
strobemers. 
* The n'th strobemer is drawn from
            sequence[n * step..n * step + order * l].
* In each strobemer there are *order* strobes.
* In each strobemer, the first strobe's index is alwoys at the very start of the
strobemer's span. It is NOT determined by the argmin of a hash function.
* The n'th strobe (not including the first strobe, of course) is the l-mer minimizer
of a window within the strobemer's span. Specifically, that window is:
        window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n
        window_end = window_start - strobe_window_length
        strobemer_span[window_start..window_end].
*/
pub fn strobemer_similarity(
    base_seq: &[char],
    mod_seq: &[char],
    distance_function: &str,        // distance function to be used (e.g. cosine similarity).
    order: usize,                   // the number of concatenated strobes in each strobemer.
    strobe_length: usize,           // the length of each strobe.
    strobe_window_gap: usize,       // the gap between strobe selection windows.
    strobe_window_length: usize,    // the size of each strobe's selection window.
    step: usize                     // the space between strobemer selection windows.
)-> Result<f64> {
    if strobe_window_length < strobe_length {
        bail!("Strobe window length ({}) is smaller than the length of the strobe ({})",
            strobe_window_length, strobe_length);
    }

    let base_strobemers = generate_strobemers(
        base_seq,
        order,
        strobe_length,
        strobe_window_gap,
        strobe_window_length,
        step
    )?;
    let base_strobemers_chars: Vec<&[char]> = base_strobemers
        .iter().map(|strobemer| strobemer.data_ref()).collect();
    let mod_strobemers = generate_strobemers(
        mod_seq,
        order,
        strobe_length,
        strobe_window_gap,
        strobe_window_length,
        step
    )?;
    let mod_strobemers_chars: Vec<&[char]> = mod_strobemers
        .iter().map(|strobemer| strobemer.data_ref()).collect();

    distance_functions::match_distance_function(base_strobemers_chars, mod_strobemers_chars, distance_function)
}

/* HELPER FUNCTION FOR strobemer_euclidean_distance()
 * This function generates a string's set of strobemers.
 */
fn generate_strobemers(
    seq: &[char],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize,
    step: usize
) -> Result<Vec<Strobemer>> {
    let mut strobemers: Vec<Strobemer> = Vec::new(); // custom struct to manage slice ownership.
    let strobemer_span = strobe_length +
        order * (strobe_window_length + strobe_window_gap);
    let last_strobemer_start_index = seq.len() - strobemer_span; // try + 1?

    for idx in (0..last_strobemer_start_index).step_by(step) {
        let strobemer_window = &seq[idx..idx + strobemer_span];
        let strobemer = generate_single_strobemer(
            strobemer_window,
            order,
            strobe_length,
            strobe_window_gap,
            strobe_window_length
        )?;
        strobemers.push(strobemer);
    }
    Ok(strobemers)
}

/* This function generates the first strobemer from a designated window. This means the first
 * strobe is always seq[0..=l]. The subsequent strobes are generated by generate_single_minimizer()
 * from the window designated by the strobemer paper.
 *
 * To generate a string's entire set of strobemers, use this function as a HELPER FUNCTION.
 * Pass it a reference to the substring you want to generate a strobemer from, with the first
 * strobemer being seq[0..=l].
 */
fn generate_single_strobemer(
    strobemer_window: &[char],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize
) -> Result<Strobemer> {
    let mut strobes: Vec<&[char]> = Vec::new();
    let first_strobe = &strobemer_window[0..strobe_length];
    strobes.push(first_strobe);
    for n in 1..order {
        let window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n;
        let window_start = window_end - strobe_window_length;
        let strobe_window = &strobemer_window[window_start..window_end];
        let strobe = generate_single_minimizer(strobe_window, strobe_length)?;
        strobes.push(strobe);
    }
    Ok(Strobemer::new(strobes))
}

/* This is the final goal of my project! */
pub fn tensor_slide_sketch(base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    w: usize
) -> Result<f64> {
    unimplemented!();
}

#[cfg(test)]
mod similarity_method_tests {
    use pretty_assertions::assert_eq;

    use crate::utils::methods::strobemer_similarity;

    #[test]
    fn strobemer_euclidean_distance_test_1() {
        let distance_function = "cosine_similarity";
        let order = 2;
        let l = 3;
        let w_min = 6;
        let w_max  = 12;
        let step = (w_max - w_min) / 2;

        let seq1: Vec<char> = "ACGTTGCA".chars().collect();
        let seq2 = seq1.clone();
        let expected = 0.0;
        let res = strobemer_similarity(&seq1, &seq2, distance_function, order, l, w_min, w_max, step);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
}