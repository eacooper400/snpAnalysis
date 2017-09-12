#!/usr/bin/perl -w
#
# ABBA_BABA.pl
#
# Perform the ABBA BABA test on SNPs in the Monarcha data set
#
# June 5, 2014
# Liz Cooper

use strict;
srand(time|$$);

# Get the names of the input and output files from the command line
# Also get the number of replicate "trees/subsamples" to do, and the number of bootstrap replicates
my ($USAGE) = "\n$0 <input.snps> <output.dstats> <num_subsamples> <num_bootstrap>
\tinput.snps = The input file in .snps format
\t\tNote: for this script to work correctly, your outgroup pop. must be the LAST population in the file
\toutput.dstats = The name of an output file to save all D statistics
\tnum_subsamples = The number of replicate trees to calculate 1 D statistic
\tnum_bootstrap = The number of bootstrap reps to check variation in D stats\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $numReps, $numBootstrap) = @ARGV;

# Save a hash of informative SNP patterns
my %patterns = (
		'1010' => 'BABA',
		'0110' => 'ABBA',
		'1020' => 'BABA',
		'1002' => 'ABBA',
		'0120' => 'ABBA',
		'0102' => 'BABA',
		'2010' => 'BABA',
		'0210' => 'ABBA',
		'2020' => 'BABA',
		'2002' => 'ABBA',
		'0220' => 'ABBA',
		'0202' => 'BABA',
		'1220' => 'ABBA',
		'1202' => 'BABA',
		'2120' => 'BABA',
		'2102' => 'ABBA',
		'2012' => 'ABBA',
		'0212' => 'BABA',
		'1212' => 'BABA',
		'2112' => 'ABBA'
		);

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Open the input file and save the first line to use in randomly selecting individuals
open (IN, $input) || die "\nUnable to open the file $input!\n";
my $firstline = <IN>;
close(IN);
my @temp = split(/\s{1,}/, $firstline);
my @pops = @temp[2..(scalar @temp -1)];
my @sizes = ();
foreach my $pop (@pops) {
  push (@sizes, (length $pop));
}

# For each replicate, select at random 1 individual from each population
# Then read through the file, find SNP sites that fit the ABBA/BABA criteria, and classify them
# Save the total count of each outside of the file loop
while ($numReps) {
  my $abba = 0;
  my $baba = 0;

  # Save an array of all informative loci
  my @loci = ();

  my @rand_IDs = ();
  foreach my $size (@sizes) {
    my $id = int(rand($size));
    push (@rand_IDs, $id);
  }

  open (IN, $input) || die "\nUnable to open the file $input!\n";
 LOOP:while (<IN>) {
    chomp $_;
    my @info = split(/\s{1,}/, $_);
    my @pstrings = @info[2..(scalar @info - 1)];

    # Check that the outgroup is fixed 
    my $outstring = $pstrings[scalar @pstrings - 1];
    if ($outstring =~ /[MRWSYK]/) {
      next LOOP;
    }
    my $ref = substr($outstring, 0, 1);
    my $teststring = $outstring;
    $teststring =~ s/$ref//g;
    unless (length $teststring == 0) {
      next LOOP;
    }
    
    my $gen_string = '';

    for (my $i = 0; $i < scalar @rand_IDs; $i++) {
      my $col = $rand_IDs[$i];
      my $snp = substr($pstrings[$i], $col, 1);
      $gen_string .= $snp;
    }
    
    # Convert the site to a string of 0 (homozygous ref), 1 (het) and 2 (homozygous alt)
    my $binstring = '';
    for (my $g = 0; $g < length $gen_string; $g++) {
      if (substr($gen_string, $g, 1) =~ /$ref/) {
	$binstring .= '0';
      } elsif (substr($gen_string, $g, 1) =~ /[MRWSYK]/) {
	$binstring .= '1';
      } else {
	$binstring .= '2';
      }
    }
    
    # Now use the hash to check if the string is informative, and whether the SNP is ABBA or BABA like
    if (exists $patterns{$binstring}) {
      my $type = $patterns{$binstring};
      push (@loci, $type);
      if ($type =~ /ABBA/) {
	$abba++;
      } elsif ($type =~ /BABA/) {
	$baba++;
      }
    }
    
  }
  close(IN);
  
  my $D = 0;
  unless (($abba + $baba) == 0) {
    $D = ($abba - $baba) / ($abba + $baba);
  }
  my $total = $abba + $baba;
  print $total, "\n";
  print OUT $D;

  # Now, do bootstrap replicates of the informative loci (with replacement)
  my $n = $numBootstrap;
  while ($n) {
    my @resampled = ();
    my $num_informative = scalar @loci;
    while ($num_informative) {
      my $randLocus = $loci[int(rand(scalar @loci))];
      push (@resampled, $randLocus);
      $num_informative--;
    }
    my $bootstrap_abba = 0;
    my $bootstrap_baba = 0;

    foreach my $locus (@resampled) {
      if ($locus =~ /ABBA/) {
	$bootstrap_abba++;
      } else {
	$bootstrap_baba++;
      }
    }
    my $bootstrap_D = 0;
    unless (($bootstrap_abba + $bootstrap_baba) == 0) {
      $bootstrap_D = ($bootstrap_abba - $bootstrap_baba) / ($bootstrap_abba + $bootstrap_baba);
    }
    print OUT ",", $bootstrap_D;
    $n--;
  }
  print OUT "\n";
    

  $numReps--;
}
close(OUT);
exit;
      
 
