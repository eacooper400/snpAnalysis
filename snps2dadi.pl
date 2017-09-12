#!/usr/bin/perl -w
#
# snps2dadi.pl
#
# Convert SNP format to dadi input (frequency spectra) format
#
# July 24, 2015

use strict;

# Take as input the name of the SNPs file, the name of the output file, and the population names from the command line
my ($USAGE) = "\n$0 <input.snps> <output.dadi.txt> <pop1,pop2,...,popN>
\tinput.snps = The input file in .snps format
\toutput.dadi.txt = The name of the data input file to create
\tA comma delimited list of population names\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $poplist) = @ARGV;

my @pops = split(/,/, $poplist);

# Open the input and output files
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Print a header line to the output file
print OUT "InGroup\t", "OutGroup\t", "Allele1\t";
foreach my $p (@pops) {
  print OUT $p, "\t";
}
print OUT "Allele2\t";
foreach my $p (@pops) {
  print OUT $p, "\t";
}
print OUT "Gene\t", "Position", "\n";

# Now, read through the input file, 
# determine which sites are FIXED in the Russells, 
# and convert those sites into dadi frequency format, then print to the output
while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  my $outgroup = $info[(scalar @info) - 1];

  # Determine if the outgroup is fixed at this site
  my $check = 0;
  my $ref = '';
  for (my $j = 0; $j < length $outgroup; $j++) {
    my $base = substr($outgroup, $j, 1);
    if (length $ref == 0) {
      unless (($base =~ /N/) || ($base =~ /[MRWSYK]/)) {
	$ref = $base;
      }
    } else {
      unless (($base =~ /$ref/) || ($base =~ /N/)) {
	$check++;
      }
    }
  }
  if ($check == 0) {
  #if (($new_len == $start_len) || ($new_len == 0)) {
    
    # If the site is fixed in the outgroup, count the frequencies of each allele in the other populations
    # Then print everything in the correct format to the out file
    my @in_groups = @info[3..(scalar @info -2)];
    print OUT "-", $info[1], "-", "\t";
    my $outgroup2 = $info[(scalar @info) - 1];
    $outgroup2 =~ s/[NMRWSYK]//ig;
    my $outallele = substr($outgroup2, 0, 1);
    unless (length $outallele == 0) {
      print OUT "-", $outallele, "-", "\t";
    } else {
      print OUT "-", $info[1], "-", "\t";
    }
    print OUT $info[1], "\t";
    foreach my $i (@in_groups) {
      my $count = 0;
      for (my $s = 0; $s < length $i; $s++) {
	if (substr($i, $s, 1) =~ /[MRWSYK]/) {
	  $count += 1;
	} elsif (substr($i, $s, 1) =~ /$info[1]/) {
	  $count += 2;
	}
      }
      print OUT $count, "\t";
    }
    print OUT $info[2], "\t";
    foreach my $i (@in_groups) {
      my $count = 0;
      for (my $s = 0; $s < length $i; $s++) {
	if (substr($i, $s, 1) =~ /[MRWSYK]/) {
	  $count += 1;
	} elsif (substr($i, $s, 1) =~ /$info[2]/) {
	  $count += 2;
	}
      }
      print OUT $count, "\t";
    }
    my ($gene, $position) = split(/\./, $info[0]);
    print OUT $gene, "\t", $position, "\n";
  } else {
    next;
  }
}
close(IN);
close(OUT);
exit;


  
