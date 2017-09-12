# snpAnalysis

/ABBA_BABA.pl <input.snps> <output.dstats> <num_subsamples> <num_bootstrap>
	input.snps = The input file in .snps format
		Note: for this script to work correctly, your outgroup pop. must be the LAST population in the file
	output.dstats = The name of an output file to save all D statistics
	num_subsamples = The number of replicate trees to calculate 1 D statistic
	num_bootstrap = The number of bootstrap reps to check variation in
	D stats


##### AFS
```perl
	AFS.pl <input.maf> <output_afs.bins>
	input.maf = A text file with minor allele frequencies, created by snps2maf or getDerived
	output.afs.bins = A text file with the allele counts for each freq. bin in each pop.
```
##### filterMAF
```perl
filterMAF.pl <input.maf> <output.filtered> <min_coverage> <min_alleleFreq> <npops>
	input.maf = The input file with maf values, created by snps2maf
	output.filtered = An output file with maf values for sites that pass filter threshold
	min_coverage = The minimum number of non-missing individuals per pop to keep a site
	min_alleleFreq = The minimum minor allele frequency to keep a site
	npops = The number of populations in the file
```

##### FST_by_MAF
```perl
FST_by_MAF_Uni.pl <infile.maf> <outfile.fst> <npops>
	infile.maf = The input file with minor allele frequencies; created by snps2maf
	outfile.fst = The output file with Fst calculated for every pair of populations
	npops = The number of populations in the file
	```
	
##### getDerived
```perl
/getDerived.pl <input.filtered> <output.derived> <OUTGROUP_NAME>
	input.filtered = An input text file with minor allele frequences (.maf)
	output.derived = An output file with ancestral and derived frequencies
	OUTGROUP_NAME = The name of the population to use as an outgroup
	```
/getSitelist.pl <input.snps> <output.list>
	input.snps = The input file in .snps format
	output.list = A filtered list of sites that meet dadi or IM requirements


/loci2IM.pl <input.loci> <output.IM.infile> <pops> <num.loci> <pop sizes> <mutation model> <inheritance scalar>
	input.loci = The input file created by makeLoci.pl
	output.IM.infile = The file to create that will be in IM input format
	pops = A comma delimited list of population names
	pop_sizes = A comma delimited list of the # individuals in each pop
	mutation model = See IM documentation
	inheritance scalar = See IM documentation

makeLoci.pl <sites.list> <consensus.fasta> <output.loci> <names.txt>
	sites.list = The sites file created by getSitelist.pl
	consensus.fasta = The reference fasta file
	output.loci = The name of the output file to create
	names.txt = A file with 1 row per population, and a comma
	delimited list of names on each row

makeRADmigrate.pl <input.snps> <output_migrate.in>

##### permuteFst
	```perl
permuteFst.pl <input.snps> <output.perm> <num_iterations> <coverage_cutoff> <MAF_minimum>
	input.snps = The input file in .snps format
	output.perm = The output file to save all permuted Fst values
	num.iterations = The number of permutations to perform
	coverage_cutoff = The min. number individuals to keep a permuted population
	MAF_minimum = The minimum MAF to keep a permuted population
```

shared_v_cov2.pl <input.counts> <output.shared>
	input.counts = A text file of shared snp types created by SNP_categories.pl
	output.shared = A output file with a line for every coverage value, and the proportion of shared sites at that cov

SNP_Categories_Uni.pl <input.maf> <output.cats> <output.counts>
	input.maf = The input file with minor allele frequencies
	output.cats = A text file with the total counts for each SNP category
	output.counts = The full SNP list with a coverage category

snps2dadi.pl <input.snps> <output.dadi.txt> <pop1,pop2,...,popN>
	input.snps = The input file in .snps format
	output.dadi.txt = The name of the data input file to create
	A comma delimited list of population names


snps2hmp.pl <input.snps> <output.hmp> <names.txt> <min_cov> <min_MAF> <Haploid or Diploid>
	input.snps = The input file in .snps format
	output.hmp = The name of a hapmap file (used by TASSEL) to create
	names.txt = A file with sample names.
		One line per POPULATION, with a comma delimited list of sample names in each line
	min_cov = The minimum number of non-missing individuals to keep a site
	min_MAF = The minimum minor allele frequency to keep a site
	Haploid or Diploid = Whether to use Heterozygot codes (diploid) or print 2 haploid strings for each individual


##### snps2maf
```perl
	snps2maf.pl <input.snps> <output.maf> <Pop_Names>
		input.snps = The input file in .snps format create by RAD data pipeline
		output.maf = A text file with the following columns:
			CHROM.POS = The physical location of the SNP
			P = The most frequent allele
			Q = The minor allele
			Pop1_MAF = One column per population with the minor allele frequency in that pop
			Pop1_Alleles = One column per pop. with the total (non-missing) allele count 
			Total_MAF = The overall minor allele frequency for the data set
	Pop_Names: Must be a comma delimited list of names; used in column
	headers
```

snps2nexus.pl <input.snps> <output.nex> <names.txt>
	input.snps = The input file in .snps format
	output.nex = The name of a nexus output file to create
	names.txt = A file with 1 row per population, and a comma delimited list of names on each row

snps2structure.pl <input.snps> <output.structure> <pop_list> <names.file> <min_coverage> <min_maf>
	input.snps = The input file in .snps format
	output.structure = The name of an output file in STRUCTURE input format
	pop_list = A numeric list of which populations to include (start at 0).  Each number corresponds to a row in:
	names.file = A text file with 1 ROW per POP, with comma delimited sample names in each row
	min_coverage = The minimum coverage per population to keep a site
	min_maf = The minimum minor allele freq. to keep a site

subSample_nex.pl <input.nex> <output.nex> <List of Subsample sizes>
	input.nex = The name of the full input nexus file
	output.nex = The name of the subsampled output nexus file
	List of Subsample Sizes = A comma delimited list of the # individuals to subsample from each population


./subsetSNPs.pl <input.snps> <output.snps> <pop1.ind1,pop1.ind2> <pop2.ind1,pop2.ind2>
	input.snps = The input file in .snps format
	output.snps = The name of an output file also in .snps format
	 A separate list of individuals to include in each subset population
		Population and Individual Numbers should be 0-based
		Any number of new output populations may be specified
		Individual IDs may be re-used
