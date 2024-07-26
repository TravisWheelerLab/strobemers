use anyhow::{bail, ensure, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use itertools::Itertools;
use rand::seq::SliceRandom;
use rand::thread_rng;


// This function simply takes a string and returns the set of k-mers.
pub fn generate_kmers(sequence: &[u8], k: usize) -> Result<Vec<Vec<char>>> {
    generate_kmers_with_step(sequence, k, 1)
}

pub fn generate_kmers_with_step(
    sequence: &[u8],
    k: usize,
    step: usize
) -> Result<Vec<Vec<char>>> {
    let mut kmers = vec![];

    if k > sequence.len() {
        bail!(
            "K ({}) is less than string length {}",
            k,
            sequence.len()
        );
    }

    for i in (0..=(sequence.len() - k)).step_by(step) {
        let kmer: Vec<char> = std::str::from_utf8(&sequence[i..i + k])?.chars().collect();
        kmers.push(kmer);
    }
    Ok(kmers)
}


pub fn generate_n_masks(k: usize, spaces: usize, n: usize) -> Vec<Vec<char>> {
    let mut rng = thread_rng();
    let mut masks = Vec::new();
    let mut initial_mask = vec!['1'; k - 2];
    initial_mask.extend(vec!['0'; spaces]);
    for _ in 0..n {
        let mut perm = initial_mask.clone();
        perm.shuffle(&mut rng);
        let mask: Vec<char> = std::iter::once('1')
            .chain(perm.iter().copied())
            .chain(std::iter::once('1'))
            .collect();
        masks.push(mask);
    }
    masks
}
pub fn generate_g_spaced_kmers_with_step(
    sequence: &[u8],
    k: usize,
    spaces: usize,
    step: usize
) -> Result<HashMap<Vec<char>, Vec<Vec<char>>>> {
    ensure!(k + spaces <= sequence.len(), format!(
        "k + spaces ({:?}) is greater than string length ({:?})", k + spaces, sequence.len())
    );
    ensure!(k > 1, format!("k ({:?}) < 2", k));

    let mut spacemers: HashMap<Vec<char>, Vec<Vec<char>>> = Default::default();

    // create an initial mask to permute. Two "care" indices will always exist
    // at the start and the end of the mask, so we don't permute those.
    let mut initial_mask = vec!['1'; k - 2];
    initial_mask.extend(vec!['0'; spaces]);
    let permutations: HashSet<_> = initial_mask.iter().permutations(initial_mask.len()).collect();
    for permutation in permutations {
        let mask: Vec<char> = std::iter::once('1')
            .chain(permutation.into_iter().copied())
            .chain(std::iter::once('1'))
            .collect();

        for i in (0..(sequence.len() - k - spaces + 1)).step_by(step) { // start of window
            let window = &sequence[i..i+k+spaces];
            let spacemer = window
                .iter()
                .zip(&mask)
                .filter_map(|(&window_char, &mask_char)| if mask_char == '1' { Some(window_char as char) } else { None })
                .collect();
            spacemers.entry(mask.clone()).or_insert(Vec::new()).push(spacemer);
        }
    }
    Ok(spacemers)
}
pub fn generate_masked_seeds(
    sequence: &[char],
    mask: &Vec<char>,
    step: usize
) -> Result<Vec<Vec<char>>> {
    let mut masked_seeds = Vec::new();
    for i in (0..sequence.len() - mask.len() + 1).step_by(step) {
        masked_seeds.push(sequence[i..i+mask.len()].to_vec()
            .iter()
            .zip(mask)
            .filter_map(|(&window_char, &mask_char)| if mask_char == '1' { Some(window_char) } else { None })
            .collect()
        );
    }

    Ok(masked_seeds)
}

/* This function generates a set of minimizers from a string. It uses my_hash_function
 * with a set seed (69). It only uses ONE hash function.
 */
pub fn generate_minimizers(
    seq: &[u8],
    k: usize,
    w: usize,
    step: usize,
    hash_function: fn(&Vec<char>, u64) -> u64,
    seed: u64
) -> Result<Vec<Vec<char>>> {
    if k > w {bail!("k ({}) > w ({})", k, w);}
    if w > seq.len() {bail!("w ({}) > sequence length ({})", w, seq.len());}

    let mut minimizers: Vec<Vec<char>> = Vec::new();
    for idx in (0..seq.len() - w + 1).step_by(step) {
        let window = &seq[idx..idx + w];
        let minimizer = generate_single_minimizer(window, k, hash_function, seed)?;
        minimizers.push(minimizer);
    }
    Ok(minimizers)
}

pub fn generate_single_minimizer(
    window: &[u8],
    k: usize,
    hash_function: fn(&Vec<char>, u64) -> u64,
    seed: u64
) -> Result<Vec<char>> {
    let kmers = generate_kmers(window, k)?;
    let minimizer = kmers
            .iter()
            .min_by_key(|&item| hash_function(item, seed))
            .unwrap();
    Ok(minimizer.to_vec())
}

/* This is a simple hash function that maps items (kmers) to u64 integers.
* It uses a seed, so this is technically a series of hash functions.
*/
fn my_hash_function (item: &Vec<char>, seed: u64) -> u64 {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish()
}

/// This function generates a string's set of strobemers. Strobemers
/// are quirky, so here's an intuitive breakdown of what strobemers
/// are and how they're generated:
///
/// Stromebers are the concatenation of strobes. There are *order*
/// strobes in each strobemer.
/// Each strobemer is drawn from a span of length strobemer_span_length =
///    strobe_length + (strobe_window_gap + strobe_window_length) * order,
/// representing the first strobe and all of the remaining strobe windows.
/// In any sequence, there are
///    (sequence.len() - strobemer_span_length) / step
/// strobemers.
///
/// The n'th strobemer is drawn from
///    sequence[n * step..n * step + order * l].
/// In each strobemer there are *order* strobes.
/// In each strobemer, the first strobe's index is alwoys at the very
/// start of the strobemer's span. It is NOT determined by the argmin
/// of a hash function.
/// The n'th strobe (not including the first strobe, of
/// course) is the l-mer minimizer of a window within the strobemer's
/// span. Specifically, that window is:
///    window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n
///    window_end = window_start - strobe_window_length
///    strobemer_span[window_start..window_end].

/* HELPER FUNCTION FOR strobemer_euclidean_distance()
 * This function generates a string's set of strobemers.
 */
pub fn generate_strobemers(
    seq: &[u8],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize,
    step: usize,
    hash_function: Option<fn(&Vec<char>, u64) -> u64>
) -> Result<Vec<Vec<char>>> {
    ensure!(strobe_window_length > strobe_length, "Strobe length is equal or greater to the window it is selected from.");

    let mut strobemers: Vec<Vec<char>> = Vec::new(); // custom struct to manage slice ownership.
    let strobemer_span = strobe_length +
        (order - 1) * (strobe_window_length + strobe_window_gap);
    if seq.len() < strobemer_span { // the sequence is just too short...
        return Ok(strobemers)
    }


    let last_strobemer_start_index = seq.len() - strobemer_span; // try + 1?

    for idx in (0..=last_strobemer_start_index).step_by(step) {
        let strobemer_window = &seq[idx..idx + strobemer_span];
        let strobemer = generate_randstrobemer(
            strobemer_window,
            order,
            strobe_length,
            strobe_window_gap,
            strobe_window_length,
            hash_function.unwrap_or(my_hash_function)
        )?;
        strobemers.push(strobemer);
    }
    Ok(strobemers)
}

pub fn generate_minstrobemer(
    strobemer_window: &[u8],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize,
    hash_function: fn(&Vec<char>, u64) -> u64
) -> Result<Vec<char>> {
    let mut strobemer: Vec<char> = Vec::new();
    let first_strobe = std::str::from_utf8(&strobemer_window[0..strobe_length])?;
    let first_strobe: Vec<char> = first_strobe.chars().collect();
    strobemer.extend(first_strobe);
    for n in 1..order {
        let window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n;
        let window_start = window_end - strobe_window_length;
        let strobe_window = &strobemer_window[window_start..window_end];
        let strobe = generate_single_minimizer(strobe_window, strobe_length, hash_function, 69)?;
        strobemer.extend(strobe);
    }
    Ok(strobemer)
}

pub fn generate_randstrobemer(
    strobemer_window: &[u8],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize,
    hash_function: fn(&Vec<char>, u64) -> u64
) -> Result<Vec<char>> {
    let mut strobemer: Vec<char> = Vec::new();
    let first_strobe = std::str::from_utf8(&strobemer_window[0..strobe_length])?;
    let first_strobe: Vec<char> = first_strobe.chars().collect();
    strobemer.extend(first_strobe);
    for n in 1..order {
        let window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n;
        let window_start = window_end - strobe_window_length;
        let strobe_window = &strobemer_window[window_start..window_end];

        let lmers = generate_kmers_with_step(strobe_window, strobe_length, 1)?;

        strobemer = lmers
            .into_iter()
            .map(|lmer| vec![strobemer.clone(), lmer].concat())
            .min_by_key(|mish_mash| hash_function(mish_mash, 69))
            .unwrap();
    }
    Ok(strobemer)
}

/* This is the final goal of my project! */
pub fn tensor_slide_sketch(_base_seq: &[char],
    _mod_seq: &[char],
    _k: usize,
    _w: usize
) -> Result<f64> {
    unimplemented!();
}

// The following contains set and vector similarity methods
pub fn frequency_vectors<'a>(base_item_bag: &'a Vec<Vec<char>>, mod_item_bag: &'a Vec<Vec<char>>
) -> (HashMap<&'a Vec<char>, f64>, HashMap<&'a Vec<char>, f64>) {
    let mut base_item_counts: HashMap<&Vec<char>, f64> = HashMap::new();
    let mut mod_item_counts: HashMap<&Vec<char>, f64> = HashMap::new();
    for item in base_item_bag {
        base_item_counts.entry(item).and_modify(|count| *count += 1.0).or_insert(1.0);
        mod_item_counts.entry(item).or_insert(0.0);
    }
    for item in mod_item_bag {
        mod_item_counts.entry(item).and_modify(|count| *count += 1.0).or_insert(1.0);
        base_item_counts.entry(item).or_insert(0.0);
    }
    (base_item_counts, mod_item_counts)
}


pub fn euclidean_distance(base_item_bag: &Vec<Vec<char>>, mod_item_bag: &Vec<Vec<char>>,
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) = frequency_vectors(base_item_bag, mod_item_bag);

    let mut distance = 0.0;
    for (item, base_count) in base_item_counts.iter() {
        let mod_count = mod_item_counts.get(item).unwrap_or(&0.0);
        distance += (base_count.max(*mod_count) - base_count.min(*mod_count)).powf(2.0);
    }
    Ok(distance.sqrt())
}

/* This function calculates the cosine similarity between two strings' kmeer vectors.
 * cosine similarity = (v1 * v2) / (||v1|| * ||v2||)
 */
pub fn cosine_similarity(base_item_bag: &Vec<Vec<char>>, mod_item_bag: &Vec<Vec<char>>
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) = frequency_vectors(base_item_bag, mod_item_bag);
    
    let mut base_magnitude = 0.0;
    let mut mod_magnitude = 0.0;
    let mut dot_product = 0.0;
    for (item, base_count) in base_item_counts.iter() {
        let mod_count = mod_item_counts.get(item).unwrap();
        base_magnitude += base_count * base_count;
        mod_magnitude += mod_count * mod_count;
        dot_product += base_count * mod_count;
    }
    if base_magnitude == 0.0 || mod_magnitude == 0.0 {
        return Ok(0.0);
    }
    Ok(dot_product / (base_magnitude.sqrt() * mod_magnitude.sqrt()))
}

pub fn jaccard_similarity(base_item_bag: &Vec<Vec<char>>, mod_item_bag: &Vec<Vec<char>>
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) = frequency_vectors(base_item_bag, mod_item_bag);

    let mut intersection = 0.0;
    let mut union = 0.0;
    for (item, base_count) in base_item_counts.iter() {
        let mod_count = mod_item_counts.get(item).unwrap();
        union += base_count + mod_count;
        intersection += base_count.min(*mod_count) * 2.0;
    }
    Ok(intersection / union)
}

pub fn minhash_similarity(base_item_bag: &Vec<Vec<char>>, mod_item_bag: &Vec<Vec<char>>
) -> Result<f64> {
    fn hash_function<T: Hash + ?Sized> (item: &T, seed: u64) -> u64 {
        let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
        item.hash(&mut hasher);
        hasher.finish()
    }
    let mut agreement_hash = 0;
    let mut total_hash = 0;
    for i in 0..64 {
        let base_signature = base_item_bag.iter().map(|item| hash_function(item, i)).min();
        let mod_signature = mod_item_bag.iter().map(|item| hash_function(item, i)).min();
        total_hash += 1;
        if base_signature == mod_signature {agreement_hash += 1}
    }
    Ok(agreement_hash as f64 / total_hash as f64)
}