extern crate alignment_free_methods;

#[cfg(test)]
mod set_similarity_tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;
    
    use alignment_free_methods::utils::similarity_methods;

    #[test]
    fn test_jaccard_similarity() -> Result<()> {
        let result = similarity_methods::jaccard_similarity(
            &vec![vec!['A'], vec!['A']],
            &vec![vec!['A'], vec!['A']]
        )?;
        let expected = 1.0;
        assert_eq!(result, expected);
        Ok(())
    }
}

#[cfg(test)]
mod something {
    use pretty_assertions::assert_eq;
    use anyhow::Result;
    
    use alignment_free_methods::utils::similarity_methods;

    #[test]
    fn test_jaccard_similarity() -> Result<()> {
        let result = similarity_methods::jaccard_similarity(
            &vec![vec!['A'], vec!['A']],
            &vec![vec!['A'], vec!['A']]
        )?;
        let expected = 1.0;
        assert_eq!(result, expected);
        Ok(())
    }
}