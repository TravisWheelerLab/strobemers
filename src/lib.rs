use anyhow::{anyhow, Result};
use cli::*;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::vec;
use std::str::FromStr;
use std::fmt;
use std::cmp::{min, max};
use std::ops::Range;
pub mod cli;

// --------------------------------------------------
// Seed Object (how do we represent a seed?)

#[derive(Debug, Clone, Eq)]
pub struct SeedObject {
    pub identifier: u64, // seed representation, such as a hash
    pub ref_index: usize, // I have no idea what this is, it's just stored here for reference.
    pub local_index: usize, // index in the sequence it's from
    pub word_indices: Vec<usize> // if non-contiguous, the start positions of its sub-strings.
}
impl PartialEq for SeedObject {
    fn eq(&self, other: &Self) -> bool {
        (self.identifier == other.identifier) &
        (self.local_index == other.local_index) &
        (self.word_indices == other.word_indices)
    }
}
impl core::hash::Hash for SeedObject {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.identifier.hash(state);
    }
}

// --------------------------------------------------
// Kmer code

// A table to map utf-8 encoded characters to our custom u2 encoding
// for nucleotides.
//
// "A" and "a" map to 0 or 0b00
// "C" and "c" map to 1 or 0b01
// "G" and "g" map to 2 or 0b10
// "T" and "t" map to 3 or 0b11
// "N" maps to 4 or 0b100.
// ^ This is POSSIBLE because u2 doesn't exist and we HAVE to use u8 anyway.
// ^ This should probably change.
// Everything else maps to 5 or 0b101.
// ^ This should probably also change.
static SEQ_NT4_TABLE: [u8; 256] = {
    let mut table = [5; 256]; // Default 5
    table[b'N' as usize] = 4;

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

// return u2-encoded nt. See SEQ_NT4_TABLE directly above this function.
pub fn nt_encoded_as_u2(nt_utf8: &u8) -> u8 {
    SEQ_NT4_TABLE[*nt_utf8 as usize]
}

// modifies kmer_bits in place.
// Shift the bit-string 2 bits to the left, slap on the 2-bit encoded nucleotide,
// and cut off the extra 2 bits
pub fn pack_nt_onto_kmer(kmer_bits: &mut u64, nucleotide: u8, mask: &u64) {
    *kmer_bits = kmer_bits.wrapping_shl(2); // make room for a new nt with bit-shift-left.
    *kmer_bits |= nucleotide as u64; // turn the first 2 bits into the nt using bitwise OR
    *kmer_bits &= mask; // erase bits over the length of k with bitwise AND
}

// Turn a sequence of DNA (u8 slice of utf-8 encoded nucleotides) into kmers.
pub fn seq_to_kmers(sequence: &[u8], kmer_args: &KmerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut kmers = vec![];
    if (kmer_args.k == 0) | (sequence.len() < kmer_args.k) {
        return Ok(kmers);
    }

    let mask: u64 = (0b1 << (2 * kmer_args.k)) - 1; // 1 repeated 2*k times
    let mut kmer_bits = 0_u64;
    let mut l = 0; // how many characters we have "loaded" (must be >= k to be a valid kmer)

    for i in 0..sequence.len() { // i is the beginning of a k-mer
        let nt_utf8 = sequence.get(i).unwrap(); // i will never exceed seq.len()
        let nt = nt_encoded_as_u2(&nt_utf8);
        if nt >= 4 {l = 0;kmer_bits = 0;} else { // restart if nt == 'N'. Idk why, its from sahlin.
            pack_nt_onto_kmer(&mut kmer_bits, nt, &mask);
            if {l += 1; l >= kmer_args.k} { // kmer_bits IS the kmer identifier!
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
    // This test demonstrates that my kmer implementation works.
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
// Hashing code. This is used in strobemers for minimization.
// any k-mer has its own unique representation as a bitpacked integer.
// The issue is that these are ordered, so minimizing a window isn't
// a straightforward process.

// Takes a u64 (a bit-representation of a seed) which we want to hash.
// The mask restricts the output hash to 0..= 2*k. I don't really know why.
pub fn hash64(u2encoded_seed: &u64, mask: &u64) -> u64 {
    let mut key = u2encoded_seed.clone();
    key = (!key).wrapping_add(key.wrapping_shl(21));
    key &= mask;
    key ^= key.wrapping_shr(24);
    key = key.wrapping_add(key.wrapping_shl(3)).wrapping_add(key.wrapping_shl(8));
    key &= mask;
    key ^= key.wrapping_shr(14);
    key = key.wrapping_add(key.wrapping_shl(2)).wrapping_add(key.wrapping_shl(4));
    key &= mask;
    key ^= key.wrapping_shr(28);
    key = key.wrapping_add(key.wrapping_shl(31));
    key &= mask;
    key
}

// maps kmer indices to hash values NOT determined by their u2-encodings
pub fn kmer_hash_index(seq: &[u8], k: &usize) -> HashMap<usize, u64> {
    let mut hash_table = HashMap::new();
    let mask = match k {
        0 => {return hash_table;},
        1..=31 => 1_u64.wrapping_shl(*k as u32 * 2) - 1,
        32 => 0_u64.wrapping_sub(1), // max value for u64
        _ => {return hash_table;},
    };

    let mut kmer_bits = 0_u64;
    
    let mut l = 0;
    for idx in 0..seq.len() { // idx never greater than seq.len()!
        let nt_utf8 = seq.get(idx).unwrap();
        let nt = nt_encoded_as_u2(nt_utf8);
        if nt < 4 { // valid character (not 'N')
            pack_nt_onto_kmer(&mut kmer_bits, nt, &mask);
            if {l += 1; l >= *k} {
                let hash = hash64(&kmer_bits, &mask);
                hash_table.insert(idx + 1 - k, hash); // it will never have an entry
            }
        }
    }
    hash_table
}


// --------------------------------------------------
// Tests for hashing methods
#[cfg(test)]
mod hashing_tests {
    use std::collections::HashMap;
    use crate::kmer_hash_index;
    use crate::hash64;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_hash64_simple() {
        assert_eq!(0b00000011, hash64(&0b00, &0b00000011)); // A
        assert_eq!(0b00000011, hash64(&0b00, &0b00001111)); // AA
        assert_eq!(0b00000011, hash64(&0b00, &0b00111111)); // AAA
        assert_eq!(0b01000011, hash64(&0b00, &0b11111111)); // AAAA

        assert_eq!(0b00000010, hash64(&0b01, &0b00000011)); // C
        assert_eq!(0b00000110, hash64(&0b01, &0b00001111)); // CA
        assert_eq!(0b00000110, hash64(&0b01, &0b00111111)); // CAA
        assert_eq!(0b10000110, hash64(&0b01, &0b11111111)); // CAAA

        assert_eq!(0b00000001, hash64(&0b10, &0b00000011)); // G
        assert_eq!(0b00001001, hash64(&0b10, &0b00001111)); // GA
        assert_eq!(0b00001001, hash64(&0b10, &0b00111111)); // GAA
        assert_eq!(0b11001001, hash64(&0b10, &0b11111111)); // GAAA

        assert_eq!(0b00000000, hash64(&0b11, &0b00000011)); // T
        assert_eq!(0b00001100, hash64(&0b11, &0b00001111)); // TA
        assert_eq!(0b00001100, hash64(&0b11, &0b00111111)); // TAA
        assert_eq!(0b00001100, hash64(&0b11, &0b11111111)); // TAAA
    }

    #[test]
    fn test_hash64_complex() {
        assert_eq!(0b0100111, hash64(&0b01001100, &0b11111111)); // CATA
        assert_eq!(0b0010011, hash64(&0b11110000, &0b11111111)); // TTAA
    }

    #[test]
    fn test_kmer_hash_table_demonstration(){
        let hash_table = kmer_hash_index(b"ACGT", &2);
        let mut true_hash_table = HashMap::new();
        true_hash_table.insert(0, hash64(&0b0001, &0b1111)); // AC
        true_hash_table.insert(1, hash64(&0b0110, &0b1111)); // CG
        true_hash_table.insert(2, hash64(&0b1011, &0b1111)); // GT
        assert_eq!(hash_table, true_hash_table);
    }
    #[test]
    fn test_kmer_hash_table_seqlen33_k32(){
        let table = kmer_hash_index(b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA", &32);
        let mut true_table = HashMap::new();
        true_table.insert(0, hash64(
            &0b0000000001010101101010101111111100000000010101011010101011111111, 
            &0b1111111111111111111111111111111111111111111111111111111111111111)
        );
        true_table.insert(1, hash64(
            &0b0000000101010110101010111111110000000001010101101010101111111100, 
            &0b1111111111111111111111111111111111111111111111111111111111111111)
        );
        assert_eq!(table, true_table);
    }

    #[test]
    fn test_kmer_hash_table_k33(){
        let table = kmer_hash_index(b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA", &33);
        let true_table = HashMap::new();
        assert_eq!(table, true_table);
    }

    #[test]
    fn test_kmer_hash_table_k0(){
        let table = kmer_hash_index(b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA", &0);
        let true_table = HashMap::new();
        assert_eq!(table, true_table);
    }
}


// --------------------------------------------------
// Strobemers code

// We have Minstrobes, Randstrobes and Hybridstrobes.
// This makes our life easier for pattern matching.
#[derive(Debug, Clone)]
pub enum Protocol {
    Rand,
    Min,
    Hybrid,
}

impl FromStr for Protocol { // for pattern matching
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

impl fmt::Display for Protocol { // for formatting/displaying to a string.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let protocol_str = match self {
            Protocol::Rand => "rand",
            Protocol::Min => "min",
            Protocol::Hybrid => "hybrid",
        };
        write!(f, "{}", protocol_str)
    }
}


pub fn seq_to_strobemers(seq: &[u8], args: &StrobemerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut strobemers: Vec<SeedObject> = Vec::new();

    if seq.len() < args.w_max {
        return Ok(strobemers);
    }

    // pre-compute hashes for l-mers at every index
    let lmer_hash_at_index: HashMap<usize, u64> = kmer_hash_index(&seq, &args.strobe_length);

    for strobemer_start_idx in 0..=lmer_hash_at_index.len() { // One strobemer for each strobe start index
        let strobemer = create_strobemer(
            &args,
            &lmer_hash_at_index,
            strobemer_start_idx,
            &seq.len()
        )?;
        strobemers.push(strobemer);
    }
    Ok(strobemers)
}

pub fn create_strobemer(args: &StrobemerSpecificArgs, lmer_hash_at_index: &HashMap<usize, u64>, strobemer_start_idx: usize, seq_len: &usize) -> Result<SeedObject> {

    // Select first strobe
    let mut strobe_indices = vec![strobemer_start_idx];
    let mut strobemer_id = lmer_hash_at_index[&strobemer_start_idx];

    let w_u = min(args.w_max, (*seq_len-strobemer_start_idx)/(args.order - 1));
    let w_l = max(args.strobe_length, args.w_min - (args.w_max - w_u));
    
    for strobe_number in 2..=args.order { // one strobe per iteration
        let strobe_window = strobemer_start_idx + w_l + (strobe_number - 2)*w_u..strobemer_start_idx + (strobe_number - 1) * w_u;
        let strobe_idx = match args.protocol {
            Protocol::Min => argmin(&lmer_hash_at_index, strobe_window, &0)?,
            Protocol::Rand => argmin(&lmer_hash_at_index, strobe_window, &strobemer_id)?,
            Protocol::Hybrid => unimplemented!()
        };
        
        strobe_indices.push(strobe_idx);
        strobemer_id = strobemer_id/3 + lmer_hash_at_index[&strobe_idx]/2;
        // ^ this may be wrong
    }

    Ok(SeedObject {
        identifier: strobemer_id,
        ref_index: args.ref_index.clone(),
        local_index: strobemer_start_idx,
        word_indices: strobe_indices,
    })

}

pub fn argmin(lmer_hash_at_index: &HashMap<usize, u64>, strobe_window: Range<usize>, mask: &u64) -> Result<usize> {
    let mut min_key = None;
    let mut min_value = u64::MAX;

    for idx in strobe_window {
        match lmer_hash_at_index.get(&idx) {
            Some(&hash) => {
                if hash ^ mask < min_value { // if mask == 0, ^ has no effect
                    min_value = hash;
                    min_key = Some(idx);
                }
            }
            None => {
                return Err(anyhow!("Could not index lmer_hash_at_index at (idx={})", idx));
            }
        }
    }

    Ok(min_key.expect("min_key was never assigned: likely that strobe_window has length 0"))
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

// --------------------------------------------------
// Strobemer tests
// --------------------------------------------------

#[cfg(test)]
mod strobemer_tests {
    use crate::create_strobemer;
    use crate::SeedObject;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_jaccard_similarity_simple() {

    }
}


// --------------------------------------------------
// Similarity methods
// --------------------------------------------------


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
// Tests for similarity methods
// --------------------------------------------------

#[cfg(test)]
mod jaccard_similarity_tests {
    use crate::jaccard_similarity;
    use crate::SeedObject;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_jaccard_similarity_simple() {
        let seed_1 = SeedObject{
            identifier: 0b00,
            ref_index: 0,
            local_index: 0,
            word_indices: vec![]
        };
        let seed_1_clone = seed_1.clone();
        let mut seed_2 = seed_1.clone();
        seed_2.local_index = 1;
        let seed_2_clone = seed_2.clone();
        let mut seed_3 = seed_1.clone();
        seed_3.local_index = 2;
        let seed_3_clone = seed_3.clone();
        let mut seed_bag1 = Vec::<SeedObject>::new();
        seed_bag1.push(seed_1);
        let mut seed_bag2 = Vec::<SeedObject>::new();
        seed_bag2.push(seed_1_clone);
        assert_eq!(1.0, jaccard_similarity(&seed_bag1, &seed_bag2).unwrap());

        seed_bag1.push(seed_2);
        assert_eq!(2.0/3.0, jaccard_similarity(&seed_bag1, &seed_bag2).unwrap());
        seed_bag1.push(seed_3);
        assert_eq!(2.0/4.0, jaccard_similarity(&seed_bag1, &seed_bag2).unwrap());

        seed_bag2.push(seed_2_clone);
        seed_bag2.push(seed_3_clone);
        assert_eq!(1.0, jaccard_similarity(&seed_bag1, &seed_bag2).unwrap());

        
    }
}

// --------------------------------------------------
// Match/NAMs related code.
pub struct Match<'a> {
    pub ref_id: &'a str,
    pub query_id: &'a str,
    pub ref_idx: usize,
    pub query_idx: usize
}


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