use anyhow::{anyhow, Result};
use cli::*;
use rand::seq;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::vec;
use std::str::FromStr;
use std::fmt;
pub mod cli;

// --------------------------------------------------
// Seed Object (how do we represent a seed?)

#[derive(Debug)]
pub struct SeedObject {
    pub identifier: u64, // seed representation, such as a hash
    pub ref_index: usize, // I have no idea what this is, it's just stored here for reference.
    pub local_index: usize, // index in the sequence it's from
    pub word_indices: Vec<usize> // if non-contiguous, the start positions of its sub-strings.
}
impl PartialEq for SeedObject {
    fn eq(&self, other: &Self) -> bool {
        self.identifier == other.identifier
    }
}
impl  Eq for SeedObject{}
impl core::hash::Hash for SeedObject {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.identifier.hash(state);
    }
}

// --------------------------------------------------
// Kmer code

// return u2-encoded nt
// it has to be a u8 because u2 type doesn't exist
pub fn nt_encoded_as_u2(nt_utf8: &u8) -> u8 {
    SEQ_NT4_TABLE[*nt_utf8 as usize]
}

pub fn pack_nt_onto_kmer(kmer_bits: &mut u64, nt: u8, mask: &u64) {
    *kmer_bits <<= 2; // make room for a new nt with bit-shift-left.
    *kmer_bits |= nt as u64; // turn the first 2 bits into the nt using bitwise OR
    *kmer_bits &= mask; // erase bits over the length of k with bitwise AND
}

// This turns a sequence of DNA (represented by a u8 slice) into kmers.
pub fn seq_to_kmers(sequence: &[u8], kmer_args: &KmerSpecificArgs) -> Result<Vec<SeedObject>> {
    if (kmer_args.k == 0) | (sequence.len() < kmer_args.k) {
        return Ok(Vec::<SeedObject>::new());
    }

    let mask: u64 = (0b1 << (2 * kmer_args.k)) - 1; // 1 repeated 2*k times
    let mut kmers = vec![];
    let mut kmer_bits = 0_u64;
    let mut l = 0; // how many characters we have "loaded" (must be >= k to be a valid kmer)

    for i in 0..sequence.len() { // i is the beginning of a k-mer
        let nt_utf8 = sequence.get(i).unwrap(); // i will never exceed seq.len()
        let nt = nt_encoded_as_u2(&nt_utf8);
        if nt >= 4 {l = 0;kmer_bits = 0;} else { // restart if nt == 'N'. Idk why, its from sahlin.
            pack_nt_onto_kmer(&mut kmer_bits, nt, &mask);
            if {l += 1; l >= kmer_args.k} { // kmer_bits represents a kmer!
                kmers.push(SeedObject {
                    identifier: kmer_bits.clone(),
                    ref_index: kmer_args.ref_index.clone(),
                    local_index: i + 1 - kmer_args.k,
                    word_indices: vec![],
                });
            }
        }
    }
    Ok(kmers)
}

// --------------------------------------------------
// Tests for kmer methods
#[cfg(test)]
mod kmer_tests {
    use anyhow::Result;
    use crate::{nt_encoded_as_u2, pack_nt_onto_kmer, seq_to_kmers, KmerSpecificArgs, SeedObject};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_nt_encoded_as_u2() {
        let seq = b"ACGTNx";
        assert_eq!(0u8 as u8, nt_encoded_as_u2(&seq[0]));
        assert_eq!(1u8 as u8, nt_encoded_as_u2(&seq[1]));
        assert_eq!(2u8 as u8, nt_encoded_as_u2(&seq[2]));
        assert_eq!(3u8 as u8, nt_encoded_as_u2(&seq[3]));
        assert_eq!(4u8 as u8, nt_encoded_as_u2(&seq[4]));
        assert_eq!(5u8 as u8, nt_encoded_as_u2(&seq[5]));
    }

    #[test]
    fn test_pack_nt_onto_kmer() {
        let mut kmer_bits = 0b0001; // AC
        pack_nt_onto_kmer(&mut kmer_bits, 0b10, &0b1111);

        let true_kmer_bits = 0b0110;
        assert_eq!(true_kmer_bits, kmer_bits);
    }
    
    #[test]
    fn test_seq_to_kmers_demonstration() -> Result<()>{
        let kmers_returned =seq_to_kmers(
            b"ACTG", 
            &KmerSpecificArgs{k:3, ref_index:0}
        )?;
        let mut kmers_real = Vec::new();
        kmers_real.push(SeedObject{
            identifier: 0b000111,
            ref_index: 0,
            local_index: 0,
            word_indices: vec![]
        });
        kmers_real.push(SeedObject{
            identifier: 0b011110,
            ref_index: 0,
            local_index: 1,
            word_indices: vec![]
        });
        assert_eq!(kmers_real, kmers_returned);
        Ok(())
    }

    #[test]
    fn test_seq_to_kmers_k0() -> Result<()>{
        let kmers_returned =seq_to_kmers(
            b"ACTG", 
            &KmerSpecificArgs{k:0, ref_index:0}
        )?;
        assert_eq!(Vec::<SeedObject>::new(), kmers_returned);
        Ok(())
    }

    #[test]
    fn test_seq_to_kmers_k_greater_than_sequence_length() -> Result<()>{
        let kmers_returned =seq_to_kmers(
            b"ACTG", 
            &KmerSpecificArgs{k:5, ref_index:0}
        )?;
        assert_eq!(Vec::<SeedObject>::new(), kmers_returned);
        Ok(())
    }
}



// --------------------------------------------------
// Hashing code -- used for fast index-to-kmer lookups in strobemers

// Takes a u64 (a bit-representation of a seed) which we want to hash.
// The mask restricts the output hash to 0..= 2*k. I don't really know why.
pub fn hash64(item: &u64, mask: &u64) -> u64 {
    let mut key = !item.saturating_add(item << 21) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = key.saturating_add(key << 3) + (key << 8) & mask; // key * 265
    key = key ^ key >> 14;
    key = key.saturating_add(key << 2) + (key << 4) & mask; // key * 21
    key = key ^ key >> 28;
    key = key.saturating_add(key << 31) & mask;
    key
}

// Lookup by index in the sequence
pub fn kmer_hash_table(seq: &[u8], k: &usize) -> HashMap<usize, u64> {
    let mut hash_at_idx = HashMap::new();
    let mask: u64 = (0b1 << (2 * k)) - 1; // 1 repeated 2*k times
    let mut kmer_bits = 0_u64;
    
    let mut l = 0;
    for idx in 0..seq.len() { // idx never greater than seq.len()!
        let nt_utf8 = seq.get(idx).unwrap();
        let nt = nt_encoded_as_u2(nt_utf8);
        if nt < 4 { // valid character (not 'N')
            pack_nt_onto_kmer(&mut kmer_bits, nt, &mask);
            if {l += 1; l >= *k} {
                let hash = hash64(&kmer_bits, &mask);
                hash_at_idx.insert(idx + 1 - k, hash); // it will never have an entry
            }
        }
    }
    hash_at_idx
}

// --------------------------------------------------
// Strobemers code

#[derive(Debug, Clone)]
pub enum Protocol {
    Rand,
    Min,
    Hybrid,
}

impl FromStr for Protocol {
    type Err = String;

    fn from_str(input: &str) -> Result<Protocol, Self::Err> {
        match input {
            "rand" => Ok(Protocol::Rand),
            "min" => Ok(Protocol::Min),
            "hybrid" => Ok(Protocol::Hybrid),
            _ => Err(format!("Invalid protocol: {}", input)),
        }
    }
}

impl fmt::Display for Protocol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let protocol_str = match self {
            Protocol::Rand => "rand",
            Protocol::Min => "min",
            Protocol::Hybrid => "hybrid",
        };
        write!(f, "{}", protocol_str)
    }
}


/// This function generates a string's set of strobemers.
/// Here's an intuitive breakdown of what strobemers
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
pub fn seq_to_strobemers(seq: &[u8], args: &StrobemerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut seed_vector: Vec<SeedObject> = Vec::new();

    if seq.len() < args.w_max {
        return Ok(seed_vector);
    }

    // pre-compute hashes for l-mers at every index
    let lmer_hash_at_index: HashMap<usize, u64> = kmer_hash_table(&seq, &args.strobe_length);

    for strobemer_start_idx in 0..=seq.len() { // One strobemer for each strobe start index
        let mut current_strobemer_hash = lmer_hash_at_index[&strobemer_start_idx];
        // ^ the first strobe is fixed. current_strobemer_hash will be mutated.
        let mut strobe_indices = vec![strobemer_start_idx];
        // ^ we will append indices as we select (order - 1) more strobes.
        for strobe_number in 0..args.order - 1 { // one strobe per iteration
            let end_of_last_window =
                strobemer_start_idx + &args.strobe_length + strobe_number * &args.w_max;
            let strobe_selection_range = {
                if end_of_last_window + &args.w_max <= seq.len() - &args.strobe_length - 1 {
                    end_of_last_window + &args.w_min..end_of_last_window + &args.w_max
                } else if end_of_last_window + &args.w_min + 1 < seq.len()
                  && seq.len() <= end_of_last_window + &args.w_max {
                    end_of_last_window + &args.w_min..seq.len() - 1
                } else {
                    return Ok(seed_vector)
                }
            };

            let (selected_strobe_index, selected_strobemer_hash) = match args.protocol {
                Protocol::Rand => select_randstrobe_index_and_hash(
                    &lmer_hash_at_index, strobe_selection_range, current_strobemer_hash)?,
                Protocol::Min => select_randstrobe_index_and_hash(
                    &lmer_hash_at_index, strobe_selection_range, 0)?,
                Protocol::Hybrid => unimplemented!(),
            };
            
            strobe_indices.push(selected_strobe_index);
            current_strobemer_hash = selected_strobemer_hash;
        }

        // construct the SeedObject vector
        seed_vector.push(SeedObject {
            identifier: current_strobemer_hash,
            ref_index: args.ref_index.clone(),
            local_index: strobemer_start_idx,
            word_indices: strobe_indices,
        });
    }
    Ok(seed_vector)
}


// Randstrobes are conditionally dependent on previous strobes' hashes. If we fix the previous
// strobes hashes, the conditional dependence disappears and it becomes equivalent to
// minstrobes.  
pub fn select_randstrobe_index_and_hash(
    lmer_hash_at_index: &HashMap<usize, u64>,
    strobe_selection_range: std::ops::Range<usize>,
    current_strobemer_hash: u64,
) -> Result<(usize, u64)> {

    let (next_strobe_index, next_strobe_hash) = strobe_selection_range
        .map(|i| (i, lmer_hash_at_index.get(&i).unwrap() ^ current_strobemer_hash))
        .min_by_key(|&(_, hash)| hash)
        .unwrap();

    Ok((next_strobe_index, next_strobe_hash))
}

// It's exactly what it sounds like.
pub fn jaccard_similarity(base_seed_bag: &Vec<SeedObject>, mod_seed_bag: &Vec<SeedObject>
) -> Result<f64> {

    let union = (base_seed_bag.len() + mod_seed_bag.len()) as f64;
    let mut intersection = 0.0;
    let mut base_hash_set = HashMap::new();

    for base_seed in base_seed_bag {
        base_hash_set.entry(base_seed)
            .and_modify(|value| *value += 1)
            .or_insert(1);
    }
    for mod_seed in mod_seed_bag {
        base_hash_set.entry(mod_seed).and_modify(|value| {
            *value -= 1;
            intersection += 2.0;
        });
    }

    Ok(intersection / union)
}

// --------------------------------------------------
// Match/NAMs related code.

#[allow(warnings)]
pub fn find_nams<'a>(
    query_seed_bag: &Vec<SeedObject>,
    ref_seed_bag: &Vec<SeedObject>,
    seed_hashmap: &mut HashMap<u64, (usize, &SeedObject)>,
    k: usize,
    //acc_to_idx: ,
    query_id: &str,
    //ref_id_map: HashMap<&str, idk>,
    filter_cutoff: usize,
) -> Result<Vec<Match<'a>>> {
    let mut hit_count_all = 0;
    let mut hit_count_reduced = 0;
    let mut total_seeds = 0;

    let mut hits: HashMap<u64, Vec<Match>> = HashMap::new(); // ref_id: 
    let mut nams = Vec::new();

    for query_seed in query_seed_bag {
        if let Entry::Occupied(mut mer) = seed_hashmap.entry(query_seed.identifier) {
            let (count, offset) = mer.get(); // get extracts the values from OccupiedEntry
            if *count <= filter_cutoff {
                //for j in offset..offset+count{ }

                
            }
        }

    }

    Ok(nams)
} 

pub struct Match<'a> {
    pub ref_id: &'a str,
    pub query_id: &'a str,
    pub ref_idx: usize,
    pub query_idx: usize
}

static SEQ_NT4_TABLE: [u8; 256] = {
    let mut table = [5u8; 256]; // Default all values to 5 (error or invalid character)

    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;

    table[b'N' as usize] = 4;

    table
};

/* Tests for kmer_hash_table and hash64. */
#[cfg(test)]
mod hashing_tests {
    use std::collections::HashMap;
    use crate::kmer_hash_table;
    use crate::hash64;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_kmer_hash_table_basic(){
        let k = 2;
        let mask = (1 << (2 * k)) - 1;
        let kmer_hash_at_idx = kmer_hash_table(b"ACGT", &k);
        let mut my_hash_hap = HashMap::new();
        my_hash_hap.insert(0, hash64(&0b0001, &mask));
        my_hash_hap.insert(1, hash64(&0b0110, &mask));
        my_hash_hap.insert(2, hash64(&0b1011, &mask));
        assert_eq!(kmer_hash_at_idx, my_hash_hap);
    }

    #[test]
    fn test_kmer_hash_table_16(){
        let k = 16;
        let mask = (1 << (2 * k)) - 1;
        let kmer_hash_at_idx = kmer_hash_table(b"AAAAAAAAAAAAAAAA", &k);
        let mut my_hash_hap = HashMap::new();
        my_hash_hap.insert(0, hash64(&0b0, &mask));
        assert_eq!(kmer_hash_at_idx, my_hash_hap);
    }

    #[test]
    fn test_kmer_hash_table_17(){
        let k = 17;
        let mask = (1 << (2 * k)) - 1;
        let kmer_hash_at_idx = kmer_hash_table(b"AAAAAAAAAAAAAAAAA", &k);
        let mut my_hash_hap = HashMap::new();
        my_hash_hap.insert(0, hash64(&0b0, &mask));
        assert_eq!(kmer_hash_at_idx, my_hash_hap);
    }

    #[test]
    fn test_hash64_basic() {
        let hash = hash64(&(0b1011 as u64), &(0b1111 as u64));
        assert_eq!(hash, 0b100);
    }
}
