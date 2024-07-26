GEN_COMP = cargo run --bin generate_comparison_data --

small_kmer_eucidean_distance:
	DB_SIZE = 3
	DB = tests/inputs/sequences_$(DB_SIZE).csv
	METHOD = kmer_euclidean_distance
	K = 6
	STEP = 3
	OUTPATH = tests/outputs/data_$(DB_SIZE)/$(METHOD)/$(K)/$(STEP)/data.csv
	mkdir -p $(OUTPATH)
	$(GEN_COMP) -s $(SMALL_DB) -e METHOD -o $(SMALL_OUT)KED_k6_step3.csv --k 6 --step 3

big_strobemer:
	DB_SIZE = 2001
	DB = tests/inputs/sequences_$(DB_SIZE).csv
	METHOD = strobemer_euclidean_distance
	K = 6
	STEP = 3
	OUTPATH = tests/outputs/data_$(DB_SIZE)/$(METHOD)/$(K)/$(STEP)/data.csv
	mkdir -p $(OUTPATH)
	$(GEN_COMP) -s $(SMALL_DB) -e METHOD -o $(SMALL_OUT)KED_k6_step3.csv --k 6 --step 3

test:
	cargo run --bin generate_comparison_data -- -s tests/inputs/sequences_data_2001.csv -e strobemer_euclidean_distance --strobemer_order 0 --strobe_length 0 --strobe_window_gap 0 --strobe_window_length 0 -o 