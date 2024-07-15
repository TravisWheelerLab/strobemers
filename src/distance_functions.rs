use anyhow::Result;
use seahash::SeaHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::DistanceFunction;

// --------------------------------------------------
pub fn match_distance_function(
    base_item_bag: Vec<&[char]>,
    mod_item_bag: Vec<&[char]>,
    distance_function: DistanceFunction,
) -> Result<f64> {
    match distance_function {
        DistanceFunction::EuclideanDistance => {
            euclidean_distance(base_item_bag, mod_item_bag)
        }
        DistanceFunction::CosineSimilarity => {
            cosine_similarity(base_item_bag, mod_item_bag)
        }
        DistanceFunction::JaccardSimilarity => {
            jaccard_similarity(base_item_bag, mod_item_bag)
        }
        DistanceFunction::MinHash => {
            minhash_similarity(base_item_bag, mod_item_bag)
        }
    }
}

// --------------------------------------------------
pub fn frequency_vectors<'a>(
    base_item_bag: Vec<&'a [char]>,
    mod_item_bag: Vec<&'a [char]>,
) -> (HashMap<&'a [char], f64>, HashMap<&'a [char], f64>) {
    let mut base_item_counts: HashMap<&[char], f64> = HashMap::new();
    let mut mod_item_counts: HashMap<&[char], f64> = HashMap::new();
    for item in base_item_bag {
        base_item_counts
            .entry(item)
            .and_modify(|count| *count += 1.0)
            .or_insert(1.0);
        mod_item_counts.entry(item).or_insert(0.0);
    }
    for item in mod_item_bag {
        mod_item_counts
            .entry(item)
            .and_modify(|count| *count += 1.0)
            .or_insert(1.0);
        base_item_counts.entry(item).or_insert(0.0);
    }
    (base_item_counts, mod_item_counts)
}

// --------------------------------------------------
pub fn euclidean_distance(
    base_item_bag: Vec<&[char]>,
    mod_item_bag: Vec<&[char]>,
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) =
        frequency_vectors(base_item_bag, mod_item_bag);

    let mut distance = 0.0;
    for (item, base_count) in base_item_counts.iter() {
        let mod_count = mod_item_counts.get(item).unwrap_or(&0.0);
        distance += (base_count.max(*mod_count) - base_count.min(*mod_count))
            .powf(2.0);
    }
    Ok(distance.sqrt())
}

// --------------------------------------------------
/* This function calculates the cosine similarity between two strings' kmeer
 * vectors.  cosine similarity = (v1 * v2) / (||v1|| * ||v2||)
 */
pub fn cosine_similarity(
    base_item_bag: Vec<&[char]>,
    mod_item_bag: Vec<&[char]>,
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) =
        frequency_vectors(base_item_bag, mod_item_bag);

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

// --------------------------------------------------
pub fn jaccard_similarity(
    base_item_bag: Vec<&[char]>,
    mod_item_bag: Vec<&[char]>,
) -> Result<f64> {
    let (base_item_counts, mod_item_counts) =
        frequency_vectors(base_item_bag, mod_item_bag);

    let mut intersection = 0.0;
    let mut union = 0.0;
    for (item, base_count) in base_item_counts.iter() {
        let mod_count = mod_item_counts.get(item).unwrap();
        union += base_count + mod_count;
        intersection += base_count.min(*mod_count);
    }
    Ok(intersection / union)
}

// --------------------------------------------------
pub fn minhash_similarity(
    base_item_bag: Vec<&[char]>,
    mod_item_bag: Vec<&[char]>,
) -> Result<f64> {
    fn hash_function<T: Hash + ?Sized>(item: &T, seed: u64) -> u64 {
        let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
        item.hash(&mut hasher);
        hasher.finish()
    }
    let mut agreement_hash = 0;
    let mut total_hash = 0;
    for i in 0..64 {
        let base_signature = base_item_bag
            .iter()
            .map(|item| hash_function(item, i))
            .min();
        let mod_signature =
            mod_item_bag.iter().map(|item| hash_function(item, i)).min();
        total_hash += 1;
        if base_signature == mod_signature {
            agreement_hash += 1
        }
    }
    Ok(agreement_hash as f64 / total_hash as f64)
}
