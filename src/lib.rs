use anyhow::{bail, ensure, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use std::os::unix::net;
use std::vec;
use itertools::Itertools;
use rand::seq::{self, SliceRandom};
use rand::thread_rng;
pub mod cli;


pub struct SeedObject {
    pub identifier: u64,
    pub ref_index: usize, // I have no idea what this is, it's just stored here for reference.
    pub local_index: usize,
    pub word_indices: Vec<usize>
}

static SEQ_NT4_TABLE: [u8; 256] = {
    let mut table = [4u8; 256]; // Default all values to 4 (error or invalid character)

    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;

    table
};

// This function takes a string represented as a slice of u8 and returns the set of k-mers.
pub fn generate_kmers(sequence: &[u8], k: usize, ref_index: usize) -> Result<Vec<SeedObject>> {
    let mut kmers = vec![];
    let mask: u64 = (1 << (2 * k)) - 1; // 2 * k 1's
    let mut kmer_bits = 0_u64;
    let mut count = 0;
    let mut l = 0; // how many characters we have "loaded"

    for i in 0..=sequence.len()-1 {
        let i_char = sequence.get(i).unwrap();
        if *i_char < 4 { // Valid nucleotide (not "N", which I have no idea what N is anyway)
            kmer_bits = (kmer_bits.wrapping_shl(2) | *i_char as u64) & mask;
            /* Step 1: kmer_bits.wrapping_shl(2)
            // * shift character-encoding bits 1 character (2 bits) up the u64 kmer representation
            // Step 2: | *i_char as u64
            // * Note: the 2 bits that wrap around will always be zero (see step 3)
            // * i_char is really just a u2. Remember that the first two bits are zero.
            // * Thus, the bitwise-OR copies the i_char information to the first two bits of the u64!
            // Step 3: & mask
            // * Cut off the character "left behind" (turn it into 00) */
            l += 1;
            if l >= k { // we have at least k characters loaded into kmer_bits
                kmers.push(SeedObject {
                    identifier: kmer_bits,
                    ref_index: ref_index,
                    local_index: i,
                    word_indices: vec![],
                });
                count += 1;
            }
        }
        else { // if there is an "N", restart
            l = 0;
            kmer_bits = 0; 
        }
    }
    Ok(kmers)
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
pub fn generate_2_order_randstrobemers(
    seq: &[u8],
    order: usize,
    strobe_length: usize,
    w_min: usize,
    w_max: usize,
    step: usize,
    hash_function: Option<fn(&Vec<char>, u64) -> u64>,
    ref_index: usize,
) -> Result<Vec<SeedObject>> {
    let seed_vector: Vec<SeedObject> = Vec::new();

    if seq.len() < w_max { // Sahlin's code -- it probably serves a purpose, IDK what tho
        return Ok(seed_vector);
    }

    // 
    let first_strobe_mask: u64 = (1 << (2 * strobe_length)) - 1;
    let x = 2_u64.pow(16) - 1;

    let lmer_hash_at_index: HashMap<usize, u64> = hash_string_to_kmers(&seq, strobe_length);

    for idx in 0..=seq.len() {
        let first_strobe_hash = lmer_hash_at_index[idx];

        // idx is the start of each strobemer.
        let strobe_start_range = {
            if idx + w_max <= seq.len() - 1 {
                idx + w_min..idx + w_max
            } else if idx + w_min + 1 < seq.len() && seq.len() <= idx + w_max {
            idx + w_min..seq.len() - 1
            } else {
                return Ok(seed_vector)
            }
        };

        let (next_strobe_index, next_strobe_hash) = strobe_start_range
            .map(|i| (i, lmer_hash_at_index.get(&i).unwrap()))
            .min_by_key(|&(_, hash)| hash ^ first_strobe_hash) // need to take first strobe into account
            .unwrap();

        seed_vector.push(SeedObject {
            identifier: first_strobe_hash/2 + next_strobe_hash/3,
            ref_index: ref_index,
            local_index: idx,
            word_indices: vec![next_strobe_index],
        });
    }
    Ok(seed_vector)
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