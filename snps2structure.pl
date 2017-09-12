#!/usr/bin/perl -w
#
# snps2structure.pl
#
# Get a STRUCTURE input file from the snps data file
#
# October 7, 2013
# Liz Cooper

use strict;

# Get from the command line arguments:
# The input and output file names,
# The list of populations to include, (0-based)
# A file with the individual names for each population
# The min. coverage per population, 
# and the min. minor allele frequency

my ($USAGE) = "\n$0 <input.snps> <output.structure> <pop_list> <names.file> <min_coverage> <min_maf>
\tinput.snps = The input file in .snps format
\toutput.structure = The name of an output file in STRUCTURE input format
\tpop_list = A numeric list of which populations to include (start at 0).  Each number corresponds to a row in:
\tnames.file = A text file with 1 ROW per POP, with comma delimited sample names in each row
\tmin_coverage = The minimum coverage per population to keep a site
\tmin_maf = The minimum minor allele freq. to keep a site\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $poplist, $namefile, $mincov, $minmaf) = @ARGV;
my @pops = split(/,/, $poplist);

# Create a hash keyed by individual number to save all the genotypes at every position
# Create an array to save all of the position IDs
# Create an array to save all of the alternative alleles
my %genotypes = ();
my @positions = ();
my @alleles = ();

# Save the previous tag number outside of the loop, to avoid linked sites
my $tag = -1;

# Open the input file and start processing the data
open (IN, $input) || die "\nUnable to open the file $input!\n";

FILELOOP: while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);

  # Check that the site is not on the same tag as the previous site
  my ($t_id, $p_id) = split(/\./, $info[0]);
  $t_id =~ s/^JH//g;
  if ($t_id == $tag) {
    next FILELOOP;
  } else {
    $tag = $t_id;
  }

  # Get the base strings for each population specified by the user
  # Check that each population meets the minimum coverage requirement
  my @strings = ();
  my $cov_check = 0;
  foreach my $pop (@pops) {
    my $basestring = $info[$pop+3];
    push (@strings, $basestring);
    my $tempstring = $basestring;
    $tempstring =~ s/N//g;
    if ((length $tempstring) < $mincov) {
      $cov_check++;
    }
  }
  unless ($cov_check == 0) {
    next FILELOOP;
  }

  # Check that the Minor Allele Frequency meets the minimum requirement
  my $longstring = join('', @strings);
  my $snp_count = 0;
  for (my $p = 0; $p < length $longstring; $p++) {
    if (substr($longstring, $p, 1) =~ /$info[2]/) {
      $snp_count += 2;
    } elsif (substr($longstring, $p, 1) =~ /[MRWSYK]/) {
      $snp_count += 1;
    }
  }
  my $maf = $snp_count/(length($longstring));
  if (($maf < $minmaf) || ($maf > (1-$minmaf))) {
    next FILELOOP;
  }

  # If the site meets all of the requirements, add the genotype to the individual string in the hash
  # Add the position to the positions list
  push (@positions, $info[0]);
  push (@alleles, $info[2]);
  for (my $s = 0; $s < scalar @strings; $s++) {
    my $popstring = $strings[$s];
    for (my $i = 0; $i < length $popstring; $i++) {
      my $id = $s . '.' . $i;
      if (exists $genotypes{$id}) {
	$genotypes{$id} .= substr($popstring, $i, 1);
      } else {
	$genotypes{$id} = substr($popstring, $i, 1);
      }
    }
  }
}
close(IN);

#my @names = ('MM_102,MM_103,MM_129,MM_130,MM_132,MM_133,MM_135,MM_165,MM_166,MM_167,MM_168,MM_169,MM_170,MM_171,MM_173,MM_174,MM_175,MM_176,MM_177,MM_178,MM_179,MM_181,MM_182,MM_183', 'SA_104,SA_105,SA_106,SA_107,SA_121,SA_122,SA_123,SA_124,SA_125,SA_082,SA_084,SA_085,SA_087,SA_091,SA_093,SA_094,SA_095,SA_096,SA_097,SA_Mt6,SA_Mt7','UU_108,UU_110,UU_115,UU_116,UU_146,UU_147,UU_148,UU_149,UU_150,UU_151,UU_152,UU_153,UU_154,UU_155,UU_156,UU_157,UU_158,UU_160,UU_161,UU_162,UU_163,UU_917,UU_918,UU_919','RI_199,RI_200,RI_202,RI_203,RI_204,RI_205,RI_206,RI_207,RI_211,RI_214,RI_215,RI_216,RI_218,RI_233,RI_234,RI_235,RI_236,RI_237,RI_238,RI_239,RI_240,RI_241');

# Get the names from the names file
open (NAME, $namefile) || die "\nUnable to open the file $namefile!\n";
my @names = ();
while (<NAME>) {
  chomp $_;
  push (@names, $_);
}
close(NAME);

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Get and sort all of the individual keys from the genotypes hash
my @keys = keys %genotypes;
my @sorted= sort {$a <=> $b} @keys;

# Get the individual name, and the two haplotype strings for each key
# Print the haplotypes to the output file
foreach my $key (@sorted) {
  my ($pop_code, $ind_code) = split(/\./, $key);
  my $name_list = $names[$pop_code];
  my @n_list = split(/,/, $name_list);
  my $name = $n_list[$ind_code];

  my @hap_1 = '';
  my @hap_2 = '';

  my $genstring = $genotypes{$key};

  for (my $g = 0; $g < length $genstring; $g++) {
    my $gen = substr($genstring, $g, 1);
    if ($gen =~ /N/) {
      push (@hap_1, '-9');
      push (@hap_2, '-9');
    } elsif ($gen =~ /[MRWSYK]/) {
      push (@hap_1, '0');
      push (@hap_2, '1');
    } elsif ($gen =~ /$alleles[$g]/) {
      push (@hap_1, '1');
      push (@hap_2, '1');
    } else {
      push (@hap_1, '0');
      push (@hap_2, '0');
    }
  }
  my $h1 = join("\t", @hap_1);
  my $h2 = join("\t", @hap_2);

  print OUT $name, "_", "1", "\t", $pop_code+1, "\t", $h1, "\n";
  print OUT $name, "_", "2", "\t", $pop_code+1, "\t", $h2, "\n";
}
close(OUT);
print scalar @positions, "\n";
exit;

  
