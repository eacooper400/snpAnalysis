#!/usr/bin/perl -w
#
# snps2maf.pl
#
# Get the minor allele frequency at every snp site
#
# September 25, 2013
# Liz Cooper

use strict;

# Get the name of the input and output files from the command line
# Get a list of population names from the command line
my ($USAGE) = "\n$0 <input.snps> <output.maf> <Pop_Names>
\tinput.snps = The input file in .snps format create by RAD data pipeline
\toutput.maf = A text file with the following columns:
\t\tCHROM.POS = The physical location of the SNP
\t\tP = The most frequent allele
\t\tQ = The minor allele
\t\tPop1_MAF = One column per population with the minor allele frequency in that pop
\t\tPop1_Alleles = One column per pop. with the total (non-missing) allele count 
\t\tTotal_MAF = The overall minor allele frequency for the data set
\tPop_Names: Must be a comma delimited list of names; used in column headers\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $popnames) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my @pops = split(/,/, $popnames);

# Print a HEADER line to the output file
print OUT "CHROM.POS\tP\tQ\t";
foreach my $pop (@pops) {
  print OUT $pop, "_MAF", "\t";
}
foreach my $pop (@pops) {
  print OUT $pop, "_Alleles", "\t";
}
print OUT "Total_MAF\n";

while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  
  
  # Calculate the minor allele counts for each population
  # Count the total number of alleles in the populations
  my %mac = ();
  my %all = ();

  for (my $i = 0; $i < scalar @pops; $i++) {
    $mac{$pops[$i]} = 0;
    $all{$pops[$i]} = 0;
    my $genstring = $info[$i+3];

    for (my $p = 0; $p < length $genstring; $p++) {
      my $base = substr($genstring, $p, 1);
      if ($base =~ /$info[1]/) {
	$all{$pops[$i]} += 2;
      } elsif ($base =~ /0/) {
	$all{$pops[$i]} += 2;
      } elsif ($base =~ /$info[2]/) {
	$mac{$pops[$i]} += 2;
	$all{$pops[$i]} += 2;
      } elsif ($base =~ /1/) {
	$mac{$pops[$i]} += 2;
	$all{$pops[$i]} += 2;
      } elsif ($base =~ /[MRWSYK]/) {
	$mac{$pops[$i]} += 1;
	$all{$pops[$i]} += 2;
      }
    }
  }

  # Calculate the allele frequency for each population, and the overall allele frequency
  my %maf = ();

  for (my $i = 0; $i < scalar @pops; $i++) {
    $maf{$pops[$i]} = 0;
    
    unless ($all{$pops[$i]} == 0) {
      $maf{$pops[$i]} = $mac{$pops[$i]}/$all{$pops[$i]};
    }
  }
  
  # Check if the alternative allele is actually the minor allele
  my $total_maf = 0;
  my $num = 0;
  my $denom = 0;

  for (my $i = 0; $i < scalar @pops; $i++) {
    $num += $mac{$pops[$i]};
    $denom += $all{$pops[$i]};
  }
  unless ($denom == 0) {
    $total_maf = $num/$denom;
  }

  # Print out the correct reference and minor allele for each population
  my $line = $info[0] . "\t";
  if ($total_maf > 0.5) {
    $line .= $info[2] . "\t" . $info[1] . "\t";
  } else {
    $line .= $info[1] . "\t" . $info[2] . "\t";
  }
  foreach my $population (@pops) {
    my $maf = $maf{$population};
    if ($total_maf > 0.5) {
      if ($all{$population} > 0) {
	my $new = 1 - $maf;
	$line .= $new . "\t";
      } else {
	$line .= "0" . "\t";
      }
    } else {
      $line .= $maf . "\t";
    }
  }
  
  # Print out the total allele counts for each population
  foreach my $population (@pops) {
    $line .= $all{$population} . "\t";
  }

  if ($total_maf > 0.5) {
    print OUT $line, 1-$total_maf, "\n"; 
  } else {
    print OUT $line, $total_maf, "\n";
  } 
}
close(IN);
close(OUT);
exit;

    
