use anyhow::Result;
use cli::*;
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
    pub identifier: u64,
    pub ref_index: usize, // I have no idea what this is, it's just stored here for reference.
    pub local_index: usize,
    pub word_indices: Vec<usize>
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

// This turns a sequence of DNA (represented by a u8 slice) into kmers.
pub fn seq_to_kmers(sequence: &[u8], kmer_args: &KmerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut kmers = vec![];
    let mask: u64 = (0b1 << (2 * &kmer_args.k)) - 1; // 2 * k 1's
    let mut kmer_bits = 0_u64;
    let mut _count = 0;
    let mut l = 0; // how many characters we have "loaded"

    for i in 0..=sequence.len()-1 {
        let char_utf8 = sequence.get(i).unwrap(); // this is utf-8 encoded (e.g. 0b65 is 'A')
        let char_u2_as_u8 = SEQ_NT4_TABLE[*char_utf8 as usize];
        if char_u2_as_u8 < 4 { // Valid nucleotide (not "N", which I have no idea what N is anyway)
            kmer_bits = (kmer_bits.wrapping_shl(2) | char_u2_as_u8 as u64) & mask;
            /* Step 1: kmer_bits.wrapping_shl(2)
            // * shift character-encoding bits 1 character (2 bits) up the u64 kmer representation
            // Step 2: | *i_char as u64
            // * Note: the 2 bits that wrap around will always be zero (see step 3)
            // * i_char is really just a u2. Remember that the first two bits are zero.
            // * Thus, the bitwise-OR copies the i_char information to the first two bits of the u64!
            // Step 3: & mask
            // * Cut off the character "left behind" (turn it into 00) */
            l += 1;
            if l >= kmer_args.k { // we have at least k characters loaded into kmer_bits
                kmers.push(SeedObject {
                    identifier: kmer_bits,
                    ref_index: kmer_args.ref_index.clone(),
                    local_index: i,
                    word_indices: vec![],
                });
                _count += 1;
            }
        }
        else { // if there is an "N", restart
            l = 0;
            kmer_bits = 0; 
        }
    }
    Ok(kmers)
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

    let mask = (1 << (2 * k)) - 1;
    let mut kmer_bits = 0_u64;
    
    let mut _count = 0;
    let mut l = 0;
    for idx in 0..seq.len() {
        let char_utf8 = seq.get(idx).unwrap();
        let char_u2_as_u8 = SEQ_NT4_TABLE[*char_utf8 as usize];

        if char_u2_as_u8 < 4 { // valid character (not 'N')
            kmer_bits = (kmer_bits << 2 | char_u2_as_u8 as u64) & mask;
            l += 1;
            if l >= *k {
                _count += 1;
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
    let mut table = [4u8; 256]; // Default all values to 4 (error or invalid character)

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
