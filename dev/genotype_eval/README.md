# Produce simulated data for testing genotyping

Only evaluate SNPs for now

## Evaluate genotyping from simulated mutations

The reference genome is irrelevant for simulated mutations, so we use a random SARS strain `USA/CA-CDPH-500076119/2022|OP769690.1|2022-08-27` from the tree.

Command line options (all required):
1. path to tree \<string\>
2. coverage to evaluate \<vector<int>\>
3. number of replicates for each coverage\<int\>
4. prefix for output files \<string\>


Example command below simulates 20 mutant replicates using the mutation spectrum model of the panman file and evalutes them at 1-10x coverage. The reads are simulated using the NovaSeq model. The output is saved to `snps_eval_novaseq/snps`. Read simulation model can be modified in `evaluate_simulated_mutations.sh` and be sure to also modify the reads per coverage in the script.

```
./evaluate_simulated_mutations.sh ../panman/sars_20000.panman 1,2,3,4,5,6,7,8,9,10 20 snps_eval_novaseq/snps
```



Output directory structure:

`prefix_refFasta/` contains the original reference genome sequence, `USA/CA-CDPH-500076119/2022|OP769690.1|2022-08-27`

`prefix_varFasta/` contains the reference genome sequence with the simulated mutations

`prefix_reads/` contains the simulated reads of the mutated reference genomes

`prefix_bam/` contains the BAM and SAM files of the simulated reads aligned to the mutated reference genomes

`prefix_vcfTrue/` contains the true VCF file of the simulated mutations

`prefix_bcf_vcf/` contains the VCF file of the simulated mutations produced by `bcftools mpileup`

~~`prefix_panmap_output/` contains the VCF file of the simulated mutations produced by `panmap`~~ (not used)

`prefix_roc_bcf_vcf/` contains the ROC file of the simulated mutations produced by `vcfroc`

~~`prefix_roc_panmap/` contains the ROC file of the simulated mutations produced by `vcfroc`~~ (not used)


### Plot ROC curves

To plot the ROC curves, run `plot_roc.py` with the appropriate prefix.

Example command below plots the ROC curves for the SNPs evaluated for 20 replicates at 1-10x coverage using the NovaSeq model. The scaling factor for the mutation matrix is 16.5.

```
python3 plot_roc.py "USA_CA-CDPH-500076119_2022|OP769690.1|2022-08-27" snps_eval_novaseq/snps 20 16.5
```