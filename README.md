# snpAnalysis
## Perl scripts for the analysis of dRAD SNP data

### Example Workflow (to run the analyses described in ![this Molecular Ecology paper](http://onlinelibrary.wiley.com/doi/10.1111/mec.14116/abstract))

#### Population Differentiation:

##### Calculate Minor Allele Frequencies:

```perl
perl snps2maf.pl "All.snps" "All.maf" "Pop1,Pop2,…,PopN"
```
##### Categorize SNPs as Shared or Private Polymorphisms 
In the paper we did this before filtering, to find the right coverage threshold, and then again after filtering for final results.

`perl SNP_Categories_Uni.pl "All.maf" "All.cats" "All.counts"`

##### Calculate Proportion Shared as Function of Pop-level:

`perl shared_v_cov2.pl "input.counts" "output.txt"`

##### Filter based on coverage and minor allele frequency:

`perl filterMAF.pl "All.maf" "Filtered.maf" MinCov MinMAF NumPops`

##### Polarize SNPs (determine if they are ancestral or derived).  

`perl getDerived.pl "Filtered.maf" "Derived.maf"`

##### Calculate the Derived Frequency Spectrum 

`perl AFS.pl "Derived.maf" "AFS_bins.txt"`

##### Calculate site-by-site Fst:

`perl FST_by_MAF_Uni.pl <infile.maf> <outfile.fst> <npops>`

##### Run Fst permutation test.  
You will need to check paths to other scripts within permuteFst.pl!

`perl permuteFst.pl <input.snps> <output.perm> <num_iterations> <coverage_cutoff> <MAF_minimum>`

##### Split the SNP file into different populations or subsets.  
This script was used to get within population Fst, but could really be used to get a subset of the `.snps` file for any reason.  In the below example, the output will be a new snps file with individuals 1-7 from population 1 in group 1, and individuals 8-10 from population 1 in group 2.  The code “0.0” refers to “pop1.indiv1.”  To get individual 4 from population 3, use the code “2.3”. The first number is pop, the number after the decimal is individual.  Note that Perl number starts at 0, not 1!  After subsetting the snps file, the same scripts for calculating MAF and Fst were re-run to get the new results.

`perl subsetSNPs.pl All.snps Subset.snps 0.0,0.1,0.2,0.3,0.4,0.5,0.6  0.7,0.8,0.9`

##### Create a STRUCTURE input file

`perl snps2structure.pl <input.snps> <output.structure> <pop_list> <names.file> <min_coverage> <min_maf>`

*  pop_list = A comma delimited list of numeric codes for populations, starting at 0, so ‘0,1,2’ refers to populations 1-3 in the .snps file
*  names.file = A text file containing all individual names for every population.  Each LINE must correspond to a single population, and consist of a comma-delimited list of individual names.  Note that the naming format must meet STRUCTURE’s requirements. Below is an example for 3 populations with 4 individuals each:
```
	FO3048,FO3162,FO3171,FO3402
	FU0261,FU0262,FU0278,FU0286
	MA3607,MA3626,MA3656,MA3665
```
##### Convert snps format to hapmap.  
The hapmap format is another text format, with rows as positions and columns as individuals.  It is very similar to the .snps file, but is designed to be read by the program ![TASSEL](http://www.maizegenetics.net/#!tassel/c17q9), which has several useful features and analytical functions.  Among other things, TASSEL will output in various file formats used by other software, including Plink.

```perl
perl snps2hmp.pl All.snps Out.hmp names.file minCov minMAF Haploid/Diploid

	names.file = The same kind of names file used for Structure script
	minCov = Population level coverage required for keeping a SNP (integer)
	minMAF = minimum overall MAF for including a SNP (given as decimal)
	Ploidy = Specify Haploid or Diploid (or just H or D); this affects the 	OUTPUT, not input.  If 	you specify haploid, the output will have 2 lines per individual instead of 1, such that there are no heterozygous SNP codes
```

#### Species Tree

##### Create a nexus format file from the SNP data (2 steps).  
Step 1 gets a list of sites from the SNP file that meet the requirements of SNAPP (i.e. no missing data)  The “names.file” is the same format described above.  

Step 2 converts those sites to Nexus format:

`perl getSitelist.pl "All.snps" "Filtered.list" "names.file"`

`perl snps2nexus.pl "Filtered.list" "Filtered.nex"`

##### Sub-sample the nexus file 
(i.e. get a random set of individuals from each population to reduce computation time).  
The list of sample sizes should be a comma delimited list of integers; desired sizes should be given in the same order of the populations in the file (e.g. 3,3,3,4 means get 3 individuals each from pops 1,2, and 3, and take 4 individuals from pop 4):

`perl subsample_nex.pl "Filtered.nex" "Sample.nex" <List of Sample Sizes>`

#### ABBA-BABA test

Perform the ABBA-BABA test using the filtered list of SNPs with no missing data.  This script will do the following: 1) randomly pick 4 individuals (1 from each population/species), 2) Determine if the site is informative, and classify the SNP as ABBA or BABA, 3) count the total number of ABBA and BABA SNPs and use this to calculate the D-statistic, 4) report a D-statistic for each sub-sampling run with bootstrap values (optional) to determine significance

`perl ABBA_BABA.pl "Filtered.list" "Out.dstats" <Num_subsamples> <Num_Bootstrap>`

*  Num_subsamples = The number of sub-sampling steps to perform
*  Num_Bootstrap = The number of bootstrap reps to perform (if desire to calculate Z-score; this is not actually something we did in the  manuscript, but has been described by others)

#### Migration Analyses

##### Create the SNPs input format file specific to dadi.  
This script assumes the last population in the SNPs file can be treated as an outgroup, and will filter for SNPs fixed in this group.  “PopList” is a comma delimited list of population names as you want them to appear in the dadi file.

`perl snps2dadi.pl All.snps dadiSNPs.txt <PopList>`

Following the dadi manual instructions, make the frequency spectrum files as follows (as done in iPython).  The remaining dadi analyses are performed in the Python scripts provided separately:

```
ipython -pylab
dd = dadi.Misc.make_data_dict('dadiSNPs.txt')
fs = dadi.Spectrum.from_data_dict(dd, ['POP1','POP2','POP3'], [24,24,24])
fs.to_file("dadiSNPs.fs")
```

##### Create an IMa2 format input file:

`perl makeLoci.pl Filtered.list Ref.fasta IMa2.loci`

`perl loci2IM.pl IMa2.loci IM.infile <pops> <num.loci> <pop sizes> <seq.length> <mutation model> <inheritance scalar>`

##### Create a Migrate input file

`perl im2migrate.pl IM.infile Migrate.infile` (Locus format)  OR
`perl makeRADmigrate.pl  All.snps Migrate.infile` (“old” SNP format)

To create LAMARC input files, it is possible to use TASSEL to export into Phylip format, which is readable by LAMARC

### List of all scripts with corresponding arguments and descriptions
1.  ABBA_BABA

```
ABBA_BABA.pl <input.snps> <output.dstats> <num_subsamples> <num_bootstrap>
	input.snps = The input file in .snps format
		Note: for this script to work correctly, your outgroup pop. must be the LAST population in the file
	output.dstats = The name of an output file to save all D statistics
	num_subsamples = The number of replicate trees to calculate 1 D statistic
	num_bootstrap = The number of bootstrap reps to check variation in D stats
```
2. AFS
```perl
	AFS.pl <input.maf> <output_afs.bins>
	input.maf = A text file with minor allele frequencies, created by snps2maf or getDerived
	output.afs.bins = A text file with the allele counts for each freq. bin in each pop.
```
3. filterMAF
```perl
filterMAF.pl <input.maf> <output.filtered> <min_coverage> <min_alleleFreq> <npops>
	input.maf = The input file with maf values, created by snps2maf
	output.filtered = An output file with maf values for sites that pass filter threshold
	min_coverage = The minimum number of non-missing individuals per pop to keep a site
	min_alleleFreq = The minimum minor allele frequency to keep a site
	npops = The number of populations in the file
```
4.  Fst_by_MAF
```perl
FST_by_MAF_Uni.pl <infile.maf> <outfile.fst> <npops>
	infile.maf = The input file with minor allele frequencies; created by snps2maf
	outfile.fst = The output file with Fst calculated for every pair of populations
	npops = The number of populations in the file
```
5. getDerived
```perl
getDerived.pl <input.filtered> <output.derived> <OUTGROUP_NAME>
	input.filtered = An input text file with minor allele frequences (.maf)
	output.derived = An output file with ancestral and derived frequencies
	OUTGROUP_NAME = The name of the population to use as an outgroup
```
6.  getSitelist
```perl
getSitelist.pl <input.snps> <output.list>
	input.snps = The input file in .snps format
	output.list = A filtered list of sites that meet dadi or IM requirements
```
7.  loci2IM
```perl
loci2IM.pl <input.loci> <output.IM.infile> <pops> <num.loci> <pop sizes> <mutation model> <inheritance scalar>
	input.loci = The input file created by makeLoci.pl
	output.IM.infile = The file to create that will be in IM input format
	pops = A comma delimited list of population names
	pop_sizes = A comma delimited list of the # individuals in each pop
	mutation model = See IM documentation
	inheritance scalar = See IM documentation
```
8.  makeLoci
```perl
makeLoci.pl <sites.list> <consensus.fasta> <output.loci> <names.txt>
	sites.list = The sites file created by getSitelist.pl
	consensus.fasta = The reference fasta file
	output.loci = The name of the output file to create
	names.txt = A file with 1 row per population, and a comma
	delimited list of names on each row
```
9.  makeRADmigrate
```perl
makeRADmigrate.pl <input.snps> <output_migrate.in>
```
Note that this script was written specifically for the data set analyzed in the paper, and will NOT automatically work with any data set (with different population numbers and sample sizes).

10.  permuteFst
```perl
permuteFst.pl <input.snps> <output.perm> <num_iterations> <coverage_cutoff> <MAF_minimum>
	input.snps = The input file in .snps format
	output.perm = The output file to save all permuted Fst values
	num.iterations = The number of permutations to perform
	coverage_cutoff = The min. number individuals to keep a permuted population
	MAF_minimum = The minimum MAF to keep a permuted population
```
11.  shared_v_cov
```perl
shared_v_cov2.pl <input.counts> <output.shared>
	input.counts = A text file of shared snp types created by SNP_categories.pl
	output.shared = A output file with a line for every coverage value, and the proportion of shared sites at that cov
```
12. SNP_Categories
```perl
SNP_Categories_Uni.pl <input.maf> <output.cats> <output.counts>
	input.maf = The input file with minor allele frequencies
	output.cats = A text file with the total counts for each SNP category
	output.counts = The full SNP list with a coverage category
```
13.  snps2dadi
```perl
snps2dadi.pl <input.snps> <output.dadi.txt> <pop1,pop2,...,popN>
	input.snps = The input file in .snps format
	output.dadi.txt = The name of the data input file to create
	A comma delimited list of population names
```
14. snps2hmp
```perl
snps2hmp.pl <input.snps> <output.hmp> <names.txt> <min_cov> <min_MAF> <Haploid or Diploid>
	input.snps = The input file in .snps format
	output.hmp = The name of a hapmap file (used by TASSEL) to create
	names.txt = A file with sample names.
		One line per POPULATION, with a comma delimited list of sample names in each line
	min_cov = The minimum number of non-missing individuals to keep a site
	min_MAF = The minimum minor allele frequency to keep a site
	Haploid or Diploid = Whether to use Heterozygot codes (diploid) or print 2 haploid strings for each individual
```
15.  snps2maf
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
16. snps2nexus
```perl
snps2nexus.pl <input.snps> <output.nex> <names.txt>
	input.snps = The input file in .snps format
	output.nex = The name of a nexus output file to create
	names.txt = A file with 1 row per population, and a comma delimited list of names on each row
```
17. snps2structure
```perl
snps2structure.pl <input.snps> <output.structure> <pop_list> <names.file> <min_coverage> <min_maf>
	input.snps = The input file in .snps format
	output.structure = The name of an output file in STRUCTURE input format
	pop_list = A numeric list of which populations to include (start at 0).  Each number corresponds to a row in:
	names.file = A text file with 1 ROW per POP, with comma delimited sample names in each row
	min_coverage = The minimum coverage per population to keep a site
	min_maf = The minimum minor allele freq. to keep a site
```
18. subSample_nex
```perl
subSample_nex.pl <input.nex> <output.nex> <List of Subsample sizes>
	input.nex = The name of the full input nexus file
	output.nex = The name of the subsampled output nexus file
	List of Subsample Sizes = A comma delimited list of the # individuals to subsample from each population
```
19. subsetSNPs
```perl
subsetSNPs.pl <input.snps> <output.snps> <pop1.ind1,pop1.ind2> <pop2.ind1,pop2.ind2>
	input.snps = The input file in .snps format
	output.snps = The name of an output file also in .snps format
	 A separate list of individuals to include in each subset population
		Population and Individual Numbers should be 0-based
		Any number of new output populations may be specified
		Individual IDs may be re-used
```
