#!/usr/bin/perl -w
#
# makeLoci_haps.pl
#
# Make a file with locus information (instead of SNPs) for individuals
# This should be able to be used with IMa2 or Migrate with minor modifications
#
# October 2, 2013
# Liz Cooper

use strict;
use Data::Dumper;

# Get the sites list file and the consensus fasta file names from the command line
# Also get the name of the output file and the file with the individual names
my ($USAGE) = "\n$0 <sites.list> <consensus.fasta> <output.loci> <names.txt>
\tsites.list = The sites file created by getSitelist.pl
\tconsensus.fasta = The reference fasta file
\toutput.loci = The name of the output file to create
\tnames.txt = A file with 1 row per population, and a comma delimited list of names on each row\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($sitesfile, $tagfile, $output, $namesfile) = @ARGV;

# Do an initial read through of the sites list to get all of the 
# tags and the positions associated with them into an array

my %sites = ();
my %genotypes = ();

open (SITES, $sitesfile) || die "\nUnable to open the file $sitesfile!\n";

while (<SITES>) {
  chomp $_;
  my @info = split(/\t/, $_);
  $info[0] =~ s/^JH//g;

  if (exists $sites{$info[0]}) {
    my @temp = @{$sites{$info[0]}};
    push (@temp, $info[1]);
    @{$sites{$info[0]}} = @temp;

    my @populations = @{$genotypes{$info[0]}};
    my @new_pops = ();
    for (my $i = 0; $i < scalar @populations; $i++) {
      my @individuals = split(/,/, $populations[$i]);
      my @new_inds = ();
      for (my $k = 0; $k < scalar @individuals; $k++) {
	my $new_string = $individuals[$k] . substr($info[$i+2], $k, 1);
	push (@new_inds, $new_string);
      }
      my $new_ind_string = join(",", @new_inds);
      push (@new_pops, $new_ind_string);
    }
    @{$genotypes{$info[0]}} = @new_pops;
      
  } else {
    my @temp = ();
    push (@temp, $info[1]);
    @{$sites{$info[0]}} = @temp;

    my @populations = @info[2..((scalar @info) -1)];
    my @new_populations = ();
    foreach my $p (@populations) {
      my @new_p = split('', $p);
      my $new_string = join(",", @new_p);
      push (@new_populations, $new_string);
    }
    @{$genotypes{$info[0]}} = @new_populations;
  }
}
close(SITES);

#print Dumper (\%genotypes);

# Save the names of the individuals in seperate arrays stored in a hash
my %names = ();

open (NAMES, $namesfile) || die "\nUnable to open the file $namesfile!\n";
my $pop_number = 0;

while (<NAMES>) {
  chomp $_;
  my @pop_array = split(/,/, $_);
  @{$names{$pop_number}} = @pop_array;
  $pop_number++;
}
close(NAMES);

# Read through the consensus file and save the tag sequence for every listed tag
my %tags = ();

open (FASTA, $tagfile) || die "\nUnable to open the file $tagfile!\n";
my $id = -1;
my $seq = '';

while (<FASTA>) {
  chomp $_;
  if ($_ =~ /^>/) {
    my @info = split(/\|/, $_);
    $info[0] =~ s/^>//;
    unless (length $seq == 0) {
      if (exists $sites{$id}) {
	$tags{$id} = $seq;
      }
    }
    $id = $info[0];
    $seq = '';
  } else {
    $seq = $_;
  }
}
if (exists $sites{$id}) {
  $tags{$id} = $seq;
}
close(FASTA);

# Now open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Start going through the list of positions in the sites array
my @site_keys = keys %sites;
my @sorted_keys = sort {$a <=> $b} @site_keys;
print scalar @sorted_keys, "\n";

foreach my $key (@sorted_keys) {
  
  # Print out the locus name before printing the individuals
  print OUT "Locus ", $key, "\n";

  # Get the list of positions and the list of genotypes for this site
  # Get the consensus sequence for the site
  my @snp_positions = @{$sites{$key}};
  my @ind_genotypes = @{$genotypes{$key}};
  my $consensus = $tags{$key};

  # Get a list of the reference bases at every SNP position
  my @ref_bases = ();
  foreach my $snp_pos (@snp_positions) {
    my $ref = substr($consensus, ($snp_pos-1), 1);
    push (@ref_bases, $ref);
  }

  # Now, go through every individual in each population, and print 2 DNA strings
  my @populations = keys %names;
  my @pop_sort = sort {$a <=> $b} @populations;

  foreach my $population (@pop_sort) {
    my @individuals = split(/,/, $ind_genotypes[$population]);
    my @name_array = @{$names{$population}};

    for (my $i = 0; $i < scalar @name_array; $i++) {
      my $string1 = $consensus;
      my $string2 = $consensus;
    
      for (my $j = 0; $j < scalar (@snp_positions); $j++) {
	my $new_base = substr($individuals[$i], $j, 1);
	if ($new_base =~ /[MRWSYK]/) {
	  my $base2 = get_snp_bases($new_base, $ref_bases[$j]);
	  substr($string1, ($snp_positions[$j] - 1), 1) = $ref_bases[$j];
	  substr($string2, ($snp_positions[$j] - 1), 1) = $base2;
	} else {
	  substr($string1, ($snp_positions[$j] - 1), 1) = $new_base;
	  substr($string2, ($snp_positions[$j] - 1), 1) = $new_base;
	}
      }
      print OUT $name_array[$i], "_1", " ", " ", $string1, "\n";
      print OUT $name_array[$i], "_2", " ", " ", $string2, "\n";
    }
  }
}

  
close(OUT);
exit;

##########################################################################
# Subroutine(s)
##########################################################################
sub get_snp_bases {
  my ($het, $ref) = @_;
  my $new = '';

  if ($ref =~ /A/) {
    if ($het =~ /M/) {
      $new = 'C';
    } elsif ($het =~ /R/) {
      $new = 'G';
    } elsif ($het =~ /W/) {
      $new = 'T';
    } else {
      print "ERROR: HET CODE IS: $het and REF BASE IS: $ref!\n";
    }
  } elsif ($ref =~ /C/) {
    if ($het =~ /M/) {
      $new = 'A';
    } elsif ($het =~ /S/) {
      $new = 'G';
    } elsif ($het =~ /Y/) {
      $new = 'T';
    } else {
      print "ERROR: HET CODE IS: $het and REF BASE IS: $ref!\n";
    }
  } elsif ($ref =~ /G/) {
    if ($het =~ /R/) {
      $new = 'A';
    } elsif ($het =~ /S/) {
      $new = 'C';
    } elsif ($het =~ /K/) {
      $new = 'T';
    } else {
      print "ERROR: HET CODE IS: $het and REF BASE IS: $ref!\n";
    }
  } elsif ($ref =~ /T/) {
    if ($het =~ /W/) {
      $new = 'A';
    } elsif ($het =~ /Y/) {
      $new = 'C';
    } elsif ($het =~ /K/) {
      $new = 'G';
    } else {
      print "ERROR: HET CODE IS: $het and REF BASE IS: $ref!\n";
    }
  }

  return($new);
}
  
