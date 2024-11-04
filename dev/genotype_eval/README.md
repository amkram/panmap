# Produce simulated data for testing genotyping

Only evaluate SNPs for now

## Evaluate genotyping from simulated mutations

The reference genome is irrelevant for simulated mutations, so we use a random SARS strain `USA/CA-CDPH-500076119/2022|OP769690.1|2022-08-27` from the tree.

Command line options (all required):
1. path to tree \<string\>
2. number of mutations \<int\>
3. coverage to evaluate \<vector<int>\>
4. number of replicates for each coverage\<int\>
5. prefix for output files \<string\>

Example command below evaluates 5 SNPs relative to the reference genome at 1x, 3x, 5x, 7x, and 9x coverage, with 10 replicates for each coverage, and outputs the results to `5_snps_eval/prefix`:

`./evaluate_simulated_mutations.sh ../panman/sars_20000.panman 5 1,3,5,7,9 10 5_snps_eval/prefix`




