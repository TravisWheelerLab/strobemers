extern crate alignment_free_methods;

#[cfg(test)]
mod kmer_tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;
    
    use std::collections::HashMap;

    use alignment_free_methods::utils::representation_methods;

    #[test]
    fn test_generate_kmers() -> Result<()> {
        let kmers = representation_methods::generate_kmers(
            &vec!['A', 'C', 'T', 'G'],
            2
        )?;
        let expected = vec![&['A', 'C'], &['C', 'T'], &['T', 'G']];
        assert_eq!(kmers, expected);
        Ok(())
    }

    #[test]
    fn test_generate_kmers_with_step() -> Result<()> {
        let kmers = representation_methods::generate_kmers_with_step(
            &vec!['A', 'C', 'T', 'G'],
            2,
            2
        )?;
        let expected = vec![&['A', 'C'], &['T', 'G']];
        assert_eq!(kmers, expected);
        Ok(())
    }

    #[test]
    fn test_generate_g_spaced_kmers_with_step1() -> Result<()> {
        let gapmers = representation_methods::generate_g_spaced_kmers_with_step(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T'],
            3,
            1,
            1
        )?;
        let mut expected = HashMap::new();
        expected.insert(
            vec!['1', '0', '1', '1'],
            vec![
                vec!['A', 'T', 'G'],
                vec!['C', 'G', 'A'],
                vec!['T', 'A', 'C'],
                vec!['G', 'C', 'T']
            ]
        );
        expected.insert(
            vec!['1', '1', '0', '1'],
            vec![
                vec!['A', 'C', 'G'],
                vec!['C', 'T', 'A'],
                vec!['T', 'G', 'C'],
                vec!['G', 'A', 'T']
            ]
        );
        assert_eq!(gapmers, expected);

        Ok(())
    }

    #[test]
    fn test_generate_g_spaced_kmers_with_step2() -> Result<()> {
        let gapmers = representation_methods::generate_g_spaced_kmers_with_step(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'],
            3,
            2,
            1
        )?;
        let mut expected = HashMap::new();
        expected.insert(
            vec!['1', '0', '0', '1', '1'],
            vec![
                vec!['A', 'G', 'A'],
                vec!['C', 'A', 'C'],
                vec!['T', 'C', 'T'], 
                vec!['G', 'T', 'G']
            ]
        );
        expected.insert(
            vec!['1', '0', '1', '0', '1'],
            vec![
                vec!['A', 'T', 'A'],
                vec!['C', 'G', 'C'],
                vec!['T', 'A', 'T'],
                vec!['G', 'C', 'G']
            ]
        );

        expected.insert(
            vec!['1', '1', '0', '0', '1'],
            vec![
                vec!['A', 'C', 'A'],
                vec!['C', 'T', 'C'],
                vec!['T', 'G', 'T'],
                vec!['G', 'A', 'G']
            ]
        );

        assert_eq!(expected, gapmers);
        Ok(())
    }

    #[test]
    fn test_generate_masked_seeds_1() -> Result<()> {
        let mask = vec!['1', '0', '1', '1'];
        let gapmers = representation_methods::generate_masked_seeds(
            &vec!['A', 'C', 'T', 'G', 'A'],
            &mask,
            1
        )?;
        let expected = vec![&['A', 'T', 'G'], &['C', 'G', 'A']];
        assert_eq!(gapmers, expected);
        Ok(())
    }


    #[test]
    fn test_generate_masked_seeds_2() -> Result<()> {
        let mask = vec!['1', '0', '1', '1'];
        let gapmers = representation_methods::generate_masked_seeds(
            &vec!['A', 'C', 'T', 'G', 'A', 'C'],
            &mask,
            2
        )?;
        let expected = vec![&['A', 'T', 'G'], &['T', 'A', 'C']];
        assert_eq!(gapmers, expected);
        Ok(())
    }
}

#[cfg(test)]
mod minimizer_tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;

    use alignment_free_methods::utils::representation_methods;

    fn my_test_hash_function (input: &Vec<char>, _seed: u64) -> u64 {
        match input[0..2] {
            ['A', 'C'] => 1,
            ['C', 'T'] => 0,
            ['T', 'G']=> 2,
            _ => 100
        }
    }


    #[test]
    fn test_generate_single_minimizer1() -> Result<()> {
        let minimizer = representation_methods::generate_single_minimizer(
            &vec!['A', 'C', 'T', 'G'],
            2,
            my_test_hash_function,
            69
        )?;
        let expected = &['C', 'T'];
        assert_eq!(minimizer, expected);
        Ok(())
    }

    #[test]
    fn test_generate_single_minimizer2() -> Result<()> {
        let minimizer = representation_methods::generate_single_minimizer(
            &vec!['A', 'C', 'C', 'G'],
            2,
            my_test_hash_function,
            69
        )?;
        let expected = &['A', 'C'];
        assert_eq!(minimizer, expected);
        Ok(())
    }

    #[test]
    fn test_generate_minimizers() -> Result<()> {
        let minimizer = representation_methods::generate_minimizers(
            &vec!['A', 'C', 'C', 'T'],
            2,
            3,
            1,
            my_test_hash_function,
            69
        )?;
        let expected = vec![&['A', 'C'], &['C', 'T']];
        assert_eq!(minimizer, expected);
        Ok(())
    }

}


#[cfg(test)]
mod strobemer_tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;

    use alignment_free_methods::utils::representation_methods;

    fn my_test_hash_function (input: &Vec<char>, _seed: u64) -> u64 {
        match input[0..2] {
            ['A', 'C'] => 0,
            ['C', 'T'] => 1,
            ['T', 'G'] => 2,
            ['G', 'A'] => 3,
            _ => 100
        }
    }

    #[test]
    fn test_generate_single_strobemer() -> Result<()> {
        let minimizer = representation_methods::generate_single_strobemer(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'],
            3, // order
            2, // strobe len
            0, // window gap
            3, // window len
            my_test_hash_function
        )?;
        let expected = &['A', 'C', 'T', 'G', 'C', 'T'];
        assert_eq!(minimizer, expected);
        Ok(())
    }

    #[test]
    fn test_generate_strobemers1() -> Result<()> {
        let minimizer = representation_methods::generate_strobemers(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'],
            3, // order
            2, // strobe len
            0, // window gap
            3, // window len
            1, // step
            my_test_hash_function
        )?;
        let expected = vec![&['A', 'C', 'T', 'G', 'C', 'T']];
        assert_eq!(minimizer, expected);
        Ok(())
    }

    #[test]
    fn test_generate_strobemers2() -> Result<()> {
        let minimizer = representation_methods::generate_strobemers(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'],
            2, // order
            2, // strobe len
            0, // window gap
            3, // window len
            2, // step
            my_test_hash_function
        )?;
        let expected = vec![&['A', 'C', 'T', 'G'], &['T', 'G', 'A', 'C']];
        assert_eq!(minimizer, expected);
        Ok(())
    }

}