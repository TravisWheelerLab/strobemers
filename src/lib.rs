//! This is a Rust implementation of Kristoffer Sahlin's [strobemers](https://github.com/ksahlin/strobemers):
//! Strobemers are alignment seeds that are robust to indels due to unfixed gaps between strobes.
//! 

use anyhow::{anyhow, Ok, Result};
use std::collections::{HashMap, hash_map::Entry};
use std::vec;
use std::{fmt, str::FromStr};
use std::cmp::{min, max};
use std::ops::Range;

pub mod cli;
use cli::*;

/// A struct for managing seed instances.
#[derive(Debug, Clone, Eq)]
pub struct SeedObject {
    /// How we identify the seed. Usually a hash or a bit-packed representation.
    pub identifier: u64, 
    /// I have no idea what this is supposed to do, I just set it to 1 for everything.
    pub ref_index: usize, 
    /// The index in the original sequence at which this seed begins.
    pub local_index: usize,
    /// The start positions of sub-seeds (for instance, strobemers contain multiple strobes).
    pub word_indices: Vec<usize> 
}
impl PartialEq for SeedObject {
    /// Think about this some more!
    fn eq(&self, other: &Self) -> bool {
        (self.identifier == other.identifier) &
        (self.local_index == other.local_index) &
        (self.word_indices == other.word_indices)
    }
}

impl core::hash::Hash for SeedObject {
    /// Used for HashSets, not comparisons 
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.identifier.hash(state);
    }
}

// --------------------------------------------------
// Kmer code

pub trait U2 {
    /// Returns the u2-encoding of a utf_8-encoded nucleotide.
    fn to_u2(&self) -> u8;
}

impl U2 for u8 {
    /// Adapted from Robert Hubley.
    fn to_u2(&self) -> u8 {
        match self {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0, // Default for unsupported characters
        }
    }
}

/// Allows for easy manipulation of slices of objects that represent nucleotides.
pub trait VecU2 {
    /// Returns a nucleotide-u2 encoding of the utf_8 encoded [`u8`].
    fn encode_u2(&self) -> Vec<u8>;
    /// Assumes the [`u8`] is u2-encoded. Returns the [`u8`] bit-packed as a u64.
    fn pack(&self) -> u64;
}

impl VecU2 for [u8] {
    fn encode_u2(&self) -> Vec<u8> {
        self
            .iter()
            .map(|nt| nt.to_u2())
            .collect()
    }

    fn pack(&self) -> u64 {
        self
            .into_iter()
            .fold(0, |packed, nt_u2| packed << 2 | *nt_u2 as u64)
    }
}

#[cfg(test)]
mod vec_u2_tests {
    use crate::U2;

    use super::VecU2;

    use pretty_assertions::assert_eq;

    #[test]
    fn test_to_vec_u2_demonstration() {
        let gold_vec: Vec<u8> = b"ACGTx".iter().map(|c| c.to_u2()).collect();
        assert_eq!(gold_vec, b"ACGTx".encode_u2());
        assert_eq!(Vec::<u8>::new(), b"".encode_u2());
    }

    #[test]
    fn test_to_u64_u2_demonstration() {
        let seq = b"ACGT";
        let seq_u2 = 0b00011011;

        assert_eq!(seq_u2, seq.encode_u2().pack());
    }
}

pub trait Packed {
    /// Packs a u2-encoded nucleotide onto &self in-place.
    /// 
    /// Step 1: Shift the bit-string 2 bits to the left
    /// 
    /// Step 2: Slap on the 2-bit encoded nucleotide
    /// 
    /// Step 3: Cut off the extra 2 bits
    fn pack_on(&mut self, nt: u8, mask: &u64);
}

impl Packed for u64 {
    fn pack_on(&mut self, nt: u8, mask: &u64) {
        // make room for a new nt with bit-shift-left.
        *self = self.wrapping_shl(2); 
        // turn the first 2 bits into the nt using bitwise OR
        *self |= nt as u64; 
        // erase bits over the length of k with bitwise AND
        *self &= mask; 
    }
}

pub fn mask(k: usize) -> u64 {
    match k {
        1..=31 => (1 << (k * 2)) - 1,
        32 => u64::MAX,
        _ => 0,
    }
}

#[cfg(test)]
mod mask_tests {
    use super::mask;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_mask_basic() {
        assert_eq!(0b11, mask(1));
        assert_eq!(0b1111, mask(2));
        assert_eq!(0b111111, mask(3));
    }
    
    #[test]
    fn test_mask_edge_cases() {
        assert_eq!(0, mask(0));
        assert_eq!(u64::MAX, mask(32));
        assert_eq!(0, mask(33))
    }
}

/// Turn a \[[u8]\] of u2-encoded, u8-represented nucleotides into a [Vec]\<[`SeedObject`]\>.
pub fn seq_to_kmers(sequence: &[u8], kmer_args: &KmerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut kmers = vec![];
    if (kmer_args.k == 0) | (sequence.len() < kmer_args.k) {
        return Ok(kmers);
    }

    let mask = mask(kmer_args.k.clone());
    let mut kmer_bits = 0_u64;
    let mut l = 0; // how many characters we have "loaded" (must be >= k to be a valid kmer)

    for i in 0..sequence.len() { // i is the beginning of a k-mer
        let nt = *sequence.get(i)
            .ok_or(anyhow!("Index exceeds sequence length (this should never happen)."))?;
        kmer_bits.pack_on(nt, &mask);
        l += 1;
        if l >= kmer_args.k {
            kmers.push(
                SeedObject {
                    identifier: kmer_bits.clone(),
                    ref_index: kmer_args.ref_index.clone(),
                    local_index: i + 1 - kmer_args.k,
                    word_indices: vec![],
                }
            );
        }
    }
    Ok(kmers)
}

// --------------------------------------------------
// Tests for kmer methods
#[cfg(test)]
mod kmer_tests {
    use anyhow::Result;
    use crate::{mask, seq_to_kmers, KmerSpecificArgs, Packed, SeedObject, U2, VecU2};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_nt_encoded_as_u2() {
        let seq = b"ACGTacgtNx";

        assert_eq!(0o0, seq[0].to_u2());
        assert_eq!(0o1, seq[1].to_u2());
        assert_eq!(0o2, seq[2].to_u2());
        assert_eq!(0o3, seq[3].to_u2());
        assert_eq!(0o0, seq[4].to_u2());
        assert_eq!(0o1, seq[5].to_u2());
        assert_eq!(0o2, seq[6].to_u2());
        assert_eq!(0o3, seq[7].to_u2());
        assert_eq!(0o0, seq[8].to_u2());
    }

    #[test]
    fn test_pack_on() {
        let mask = &mask(3);
        let mut kmer_bits = 0_u64;

        kmer_bits.pack_on(0b00, mask);
        assert_eq!(kmer_bits, b"A".encode_u2().pack());

        kmer_bits.pack_on(0b01, mask);
        assert_eq!(kmer_bits, b"AC".encode_u2().pack());

        kmer_bits.pack_on(0b10, mask);
        assert_eq!(kmer_bits, b"ACG".encode_u2().pack());

        kmer_bits.pack_on(0b11, mask);
        assert_eq!(kmer_bits, b"CGT".encode_u2().pack());

        kmer_bits.pack_on(0b11, mask);
        assert_eq!(kmer_bits, b"GTT".encode_u2().pack());
    }
    
    #[test]
    // This test demonstrates that my kmer implementation works.
    fn test_seq_to_kmers_demonstration() -> Result<()>{
        let kmers_returned =seq_to_kmers(
            &b"ACTG".encode_u2(), 
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

/// Hashes a bitpacked seed for use in minimization methods.
/// 
/// Not lexicographically ordered like the bitpacked integer representation.
/// 
/// mask is a series of k 1's so the proper seed is hashed (a u64 always represents
/// 32 bitpacked nucleotides, so we need to mask the u64 so hash64 knows which nucleotides
/// to pay attention to).
/// 
/// Adapted from Kristoffer Sahlin: https://github.com/ksahlin/strobemers/blob/main/strobemers_cpp/index.cpp.
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

/// Map seq's indices to the [`hash64`] hash value of the index's k-mer.
/// 
/// seq is a u2-encoded, u8-represented nucleotide. 
pub fn kmer_hash_table(seq: &[u8], k: &usize) -> HashMap<usize, u64> {
    let mut hash_table = HashMap::new();
    let mask = mask(k.clone());
    let mask = match mask {
        0 => {return hash_table;}, // I feel kinda scummy about this
        _ => mask,
    };

    let mut kmer_bits = 0_u64;
    
    let mut l = 0;
    for idx in 0..seq.len() {
        let nt = seq.get(idx)
            .ok_or(anyhow!("Index exceeds sequence length (this should never happen).")).ok().unwrap();
        kmer_bits.pack_on(*nt, &mask);
        l += 1;
        if l >= *k {
            let hash = hash64(&kmer_bits, &mask);
            hash_table.insert(idx + 1 - k, hash);
        }
    }
    
    hash_table
}


// --------------------------------------------------
// Tests for hashing methods
#[cfg(test)]
mod hashing_tests {
    use std::collections::HashMap;
    use crate::{kmer_hash_table, hash64, mask, VecU2};
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
        let seq_u2: Vec<u8> = b"ACGT".encode_u2();
        let hash_table = kmer_hash_table(&seq_u2, &2);
        let mut true_hash_table = HashMap::new();
        true_hash_table.insert(0, hash64(&0b0001, &0b1111)); // AC
        true_hash_table.insert(1, hash64(&0b0110, &0b1111)); // CG
        true_hash_table.insert(2, hash64(&0b1011, &0b1111)); // GT
        assert_eq!(hash_table, true_hash_table);
    }
    #[test]
    fn test_kmer_hash_table_seqlen33_k32(){
        let seq_u2: Vec<u8> = b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA".encode_u2();
        let table = kmer_hash_table(&seq_u2, &32);
        let mut true_table = HashMap::new();
        true_table.insert(0, hash64(
            &b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".encode_u2().pack(),
            &mask(32))
        );
        true_table.insert(1, hash64(
            &b"AAACCCCGGGGTTTTAAAACCCCGGGGTTTTA".encode_u2().pack(), 
            &mask(32))
        );
        assert_eq!(table, true_table);
    }

    #[test]
    fn test_kmer_hash_table_k33(){
        let seq_u2 = b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA".encode_u2();
        let table = kmer_hash_table(&seq_u2, &33);
        let true_table = HashMap::new();
        assert_eq!(table, true_table);
    }

    #[test]
    fn test_kmer_hash_table_k0(){
        let seq_u2 = b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA".encode_u2();
        let table = kmer_hash_table(&seq_u2, &0);
        let true_table = HashMap::new();
        assert_eq!(table, true_table);
    }
}


// --------------------------------------------------
// Strobemers code

/// This enum makes our life easier for pattern matching on strobemer protocol types.
/// 
/// There are 3 strobemer protocols: Minstrobes, Randstrobes and Hybridstrobes.
#[derive(Debug, Clone)]
pub enum Protocol {
    Rand,
    Min,
    Hybrid,
}

impl FromStr for Protocol { // for pattern matching
    type Err = anyhow::Error;

    fn from_str(input: &str) -> Result<Protocol, Self::Err> {
        match input {
            "rand" => Ok(Protocol::Rand),
            "min" => Ok(Protocol::Min),
            "hybrid" => Ok(Protocol::Hybrid),
            _ => Err(anyhow!("Invalid protocol: {}", input))
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

/// Generate a sequence's set of strobemers. 
pub fn seq_to_strobemers(seq: &[u8], args: &StrobemerSpecificArgs) -> Result<Vec<SeedObject>> {
    let mut strobemers: Vec<SeedObject> = Vec::new();

    // if seq.len() < args.w_max {
    //     return Ok(strobemers);
    // }    

    // pre-compute hashes for l-mers at every index
    let lmer_hash_at_index: HashMap<usize, u64> = kmer_hash_table(&seq, &args.strobe_length);

    let seq_len = lmer_hash_at_index.len() - 1 + args.strobe_length;
    let strobemer_len = args.strobe_length * args.order;
    if seq_len < strobemer_len {
        return Ok(strobemers)
    }
    for strobemer_start_idx in 0..=seq_len - strobemer_len { // One strobemer for each strobe start index
        let strobemer = create_strobemer(
            &args,
            &lmer_hash_at_index,
            strobemer_start_idx
        )?;
        strobemers.push(strobemer);
    }
    Ok(strobemers)
}

/// Generate a sequence's strobemer at a specified start index.
/// 
/// Instead of taking the sequence, this method takes a HashMap mapping indices to kmer hashes
/// to be minimized over.
pub fn create_strobemer(args: &StrobemerSpecificArgs, lmer_hash_at_index: &HashMap<usize, u64>, strobemer_start_idx: usize) -> Result<SeedObject> {
    // Select first strobe
    let mut strobe_indices = vec![strobemer_start_idx];
    let mut strobemer_id = lmer_hash_at_index[&strobemer_start_idx];

    // let seq_length = lmer_hash_at_index.len() - 1 + args.strobe_length;
    let w_max_alt = match args.order {
        1 => 0, // does not matter, won't be used
        _ => (lmer_hash_at_index.len() - strobemer_start_idx)/(args.order-1) // comes into play at the end
        // ^ modification from Sahlin's algorithm
    };
    let w_u = min(args.w_max, w_max_alt);
    // From Sahlin's algorithm:
    // w_u = min(args.w_max, ((lmer_hash_at_index.len() + args.strobe_length - 1)-strobemer_start_idx)/(args.order - 1));
    
    let conditional_w_min = match args.w_min < args.w_max - w_u {
        true => 0, // prevent subtraction with overflow
        false => args.w_min - (args.w_max - w_u)
    };
    // From Sahlin's program:
    // w_l = max(args.strobe_length, args.w_min - (args.w_max - w_u));
    let w_l = max(conditional_w_min, args.strobe_length);
    
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
    Ok(
        SeedObject {
            identifier: strobemer_id,
            ref_index: args.ref_index.clone(),
            local_index: strobemer_start_idx,
            word_indices: strobe_indices,
        }
    )

}

/// Returns the index corresponding to the smallest hash.
/// 
/// mask is a bitstring which deterministically disrupts the l-mer's hash value (useful for
/// yet-unimplemented randstrobes)
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
    return match min_key {
        Some(key) => Ok(key),
        None => Err(anyhow!("min_key was never assigned: likely that strobe_window has length 0"))
    }

}

// --------------------------------------------------
// Strobemer tests
// --------------------------------------------------

#[cfg(test)]
mod strobemer_tests {
    use crate::{argmin, create_strobemer, hash64, kmer_hash_table, mask, seq_to_strobemers, SeedObject, StrobemerSpecificArgs, VecU2};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_argmin_demonstration() {
        let seq_u2= b"ACGTACGT".encode_u2();
        let lmer_hash_at_index = kmer_hash_table(&seq_u2, &4);
        let min_lmer_idx = argmin(&lmer_hash_at_index, 0..5, &0).ok().unwrap();
        
        let min_hash = [b"ACGT", b"CGTA", b"GTAC", b"TACG", b"ACGT"]
            .iter()
            .map(|kmer| hash64(&kmer.encode_u2().pack(), &mask(4)))
            .min()
            .unwrap();
        assert_eq!(lmer_hash_at_index[&min_lmer_idx], min_hash)
    }

    #[test]
    fn test_argmin_basic() {
        let seq = b"ACTG".encode_u2();
        let lmer_hash_at_index = kmer_hash_table(&seq, &4);
        assert_eq!(lmer_hash_at_index[&0], hash64(&seq.pack(), &mask(4)));
    }

    #[test]
    fn test_argmin_failure_out_of_range() {
        let lmer_hash_at_index = kmer_hash_table(&b"ACGTACGT".encode_u2(), &4);
        let min_lmer_idx = argmin(&lmer_hash_at_index, 0..6, &0);
        assert!(min_lmer_idx.is_err());
    }

    #[test]
    fn test_argmin_failure_empty_window() {
        let seq_u2 = b"ACGTACGT".encode_u2();
        let lmer_hash_at_index = kmer_hash_table(&seq_u2, &4);
        let min_lmer_idx = argmin(&lmer_hash_at_index, 0..0, &0);
        assert!(min_lmer_idx.is_err());
    }

    #[test]
    fn test_kmer_hash_table_reflects_seq_len() {
        let lmer_hash_at_index = kmer_hash_table(&b"AACCTTGG".encode_u2(), &5);
        assert_eq!(lmer_hash_at_index.len(), 8-5+1);


        let lmer_hash_at_index_2 = kmer_hash_table(&b"AACCTTGGAAAA".encode_u2(), &2);
        assert_eq!(lmer_hash_at_index_2.len(), 12-2+1);
    }

    #[test]
    fn test_create_strobemer_simple() {
        // order 1, length 5 string, strobe_length = 4. 2 strobemers, each is only 1 strobe.
        let args = StrobemerSpecificArgs{
            protocol: crate::Protocol::Min,
            order: 1,
            strobe_length: 4,
            w_min: 2,
            w_max: 10,
            ref_index: 0
        };

        let seq: Vec<u8> = b"ACGTA".encode_u2();
        let lmer_hash_at_index = &kmer_hash_table(&seq, &args.strobe_length);
        let mask = &mask(4);
        println!("Index 1: {}", hash64(&0b00011011, mask));

        let strobemer_idx1_gold = SeedObject {
            identifier: hash64(&0b00011011, mask),
            ref_index: args.ref_index.clone(),
            local_index: 0,
            word_indices: vec![0],
        };
        
        let strobemer_idx2_gold = SeedObject {
            identifier: hash64(&0b01101100, mask),
            ref_index: args.ref_index.clone(),
            local_index: 1,
            word_indices: vec![1],
        };
        

        let strobemer_idx1 = create_strobemer(&args, lmer_hash_at_index, 0).ok().unwrap();
        let strobemer_idx2 = create_strobemer(&args, lmer_hash_at_index, 1).ok().unwrap();

        assert_eq!(strobemer_idx1, strobemer_idx1_gold);
        assert_eq!(strobemer_idx2, strobemer_idx2_gold);
    }

    #[test]
    fn test_create_strobemer_order2() {
        // length 12 string, strobe length 4, order 2. Second strobe window range from 6 to 8.
        // It just fits because lmer_hash_at_index[8] needs a 12-length string
        // there must be a 4-mer at index 8, which requires a 12-length string.
        let args = StrobemerSpecificArgs{
            protocol: crate::Protocol::Min,
            order: 2,
            strobe_length: 4,
            w_min: 6,
            w_max: 8,
            ref_index: 0
        };

        let seq = b"ACGTACGTACGT".encode_u2();
        let lmer_hash_at_index = &kmer_hash_table(&seq, &args.strobe_length);

        let strobe_1_hash = hash64(&0b00011011, &mask(4));

        let strobe_2_window = args.w_min..args.w_max;
        let strobe_2_index = argmin(&lmer_hash_at_index, strobe_2_window, &0).ok().unwrap();
        let strobe_2_hash = lmer_hash_at_index[&strobe_2_index];

        let strobemer_idx1_gold = SeedObject {
            identifier: {
                strobe_1_hash / 3 +
                strobe_2_hash / 2
            },
            ref_index: args.ref_index.clone(),
            local_index: 0,
            word_indices: vec![0, strobe_2_index,],
        };

        let strobemer_idx1 = create_strobemer(&args, lmer_hash_at_index, 0).ok().unwrap();
        assert_eq!(strobemer_idx1, strobemer_idx1_gold);
    }

    #[test]
    fn test_create_strobemers_order1() {
        // length 12 string, strobe length 4, order 2. Second strobe window range from 6 to 8.
        // It just fits because lmer_hash_at_index[8] needs a 12-length string
        // there must be a 4-mer at index 8, which requires a 12-length string.
        let args = StrobemerSpecificArgs{
            protocol: crate::Protocol::Min,
            order: 1,
            strobe_length: 4,
            w_min: 6,
            w_max: 8,
            ref_index: 0
        };

        let seq = b"ACGTACGT".encode_u2();
        let lmer_hash_at_index = &kmer_hash_table(&seq, &args.strobe_length);
        let strobemers_gold: Vec<SeedObject> = (0..5).into_iter()
            .map(|i|
                SeedObject {
                    identifier: lmer_hash_at_index[&i],
                    ref_index: args.ref_index.clone(),
                    local_index: i,
                    word_indices: vec![i],
                }
            )
            .collect();

        let strobemers = seq_to_strobemers(&seq, &args).ok().unwrap();
        assert_eq!(strobemers, strobemers_gold);
    }

    #[test]
    fn test_create_strobemers_order2() {
        // length 12 string, strobe length 4, order 2. Second strobe window range from 6 to 8.
        // It just fits because lmer_hash_at_index[8] needs a 12-length string
        // there must be a 4-mer at index 8, which requires a 12-length string.
        let args = StrobemerSpecificArgs{
            protocol: crate::Protocol::Min,
            order: 2,
            strobe_length: 4,
            w_min: 6,
            w_max: 8,
            ref_index: 0
        };

        let seq = b"ACGTACGTACGT".encode_u2();
        let lmer_hash_at_index = &kmer_hash_table(&seq, &args.strobe_length);
        let strobemers_gold: Vec<SeedObject>= [
            (b"ACGT", 0, 6..8),
            (b"CGTA", 1, 7..9),
            (b"GTAC", 2, 7..9),
            (b"TACG", 3, 7..9),
            (b"ACGT", 4, 8..9),
        ]
            .iter()
            .map(|(strobe1, strobe1_index, strobe2_window)|
                {
                    let strobe2_index = argmin(&lmer_hash_at_index, strobe2_window.clone(), &0).ok().unwrap();
                    SeedObject {
                        identifier: {
                            hash64(&strobe1.encode_u2().pack(), &mask(4)) / 3 +
                            lmer_hash_at_index[&strobe2_index] / 2
                        },
                        ref_index: args.ref_index.clone(),
                        local_index: strobe1_index.clone(),
                        word_indices: vec![*strobe1_index, strobe2_index,]
                    }
                }
            )
            .collect();

        let strobemers = seq_to_strobemers(&seq, &args).ok().unwrap();
        assert_eq!(strobemers, strobemers_gold);
    }
}


// --------------------------------------------------
// Similarity methods
// --------------------------------------------------

/// Jaccard similarity of multi-sets. Relies on SeedObject implementing PartialEq. 
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

/// A struct for managing matching-seed related code. Unimplemented.
pub struct Match<'a> {
    pub ref_id: &'a str,
    pub query_id: &'a str,
    pub ref_idx: usize,
    pub query_idx: usize
}

/// Find non-overlapping approximate matches given two [`Vec`]\<[`SeedObject`]>'s. Unimplemented.
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