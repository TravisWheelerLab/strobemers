extern crate alignment_free_methods;

#[cfg(test)]
mod set_similarity_tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;
    
    use alignment_free_methods::{jaccard_similarity, seq_to_kmers};

    #[test]
    fn test_jaccard_similarity() -> Result<()> {
        let result = jaccard_similarity(
            &seq_to_kmers(b"AAAA", 3, 0)?,
            &seq_to_kmers(b"AAAA", 3, 0)?
        )?;
        let expected = 1.0;
        assert_eq!(result, expected);
        Ok(())
    }
}