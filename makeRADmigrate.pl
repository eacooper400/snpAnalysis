#!/usr/bin/perl -w
#
# makeRADmigrate_in.pl
#
# Make a Migrate input file using the RAD SNP data
#
# September 27, 2013

use strict;

# Get the input SNP file and the name of the outfile 
my ($USAGE) = "$0 <input.snps> <output_migrate.in>\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Store the code names for the individuals in each population
my @Mak = ('MM_102', 'MM_103', 'MM_129', 'MM_130', 'MM_132', 'MM_133', 'MM_135', 'MM_165', 'MM_166', 'MM_167', 'MM_168', 'MM_169', 'MM_170', 'MM_171', 'MM_173', 'MM_174', 'MM_175', 'MM_176', 'MM_177', 'MM_178', 'MM_179', 'MM_181', 'MM_182', 'MM_183');
my @SA = ('SA_104','SA_105','SA_106','SA_107','SA_121','SA_122','SA_123','SA_124','SA_125','SA_082','SA_084','SA_085','SA_087','SA_091','SA_093','SA_094','SA_095','SA_096','SA_097','SA_Mt6','SA_Mt7');
my @Ugi = ('UU_108', 'UU_110','UU_115','UU_116','UU_146','UU_147','UU_148','UU_149','UU_150','UU_151','UU_152','UU_153','UU_154','UU_155','UU_156','UU_157','UU_158','UU_160','UU_161','UU_162','UU_163','UU_917','UU_918','UU_919');
my @Rus = ('RI_199','RI_200','RI_202','RI_203','RI_204','RI_205','RI_206','RI_207','RI_211','RI_214','RI_215','RI_216','RI_218','RI_233','RI_234','RI_235','RI_236','RI_237','RI_238','RI_239','RI_240','RI_241');

# Open a temporary output file for each population
open (MAK, ">makira.migrate") || die "\nUnable to open the file makira.migrate!\n";
open (SA, ">santaana.migrate") || die "\nUnable to open the file santaana.migrate!\n";
open (UGI, ">ugi.migrate") || die "\nUnable to open the file ugi.migrate!\n";
open (RUS, ">russells.migrate") || die "\nUnable to open the file russells.migrate!\n";

# Create an array and a flag to temporarily store lines from the same tag
my $tag = 1000000000000;
my @lines = ();

# Open the input file and start reading through it
open (IN, $input) || die "\nUnable to open the file $input!\n";

LOOP: while (<IN>) {
  chomp $_;
  my ($id, $p, $q, $mstring, $sstring, $ustring, $rstring) = split(/\t/, $_);

  # Skip indels
  if (($mstring =~ /(0|1)/) || ($sstring =~ /(0|1)/) || ($ustring =~ /(0|1)/) || ($rstring =~ /(0|1)/)) {
    next LOOP;
  }

  # Skip sites with overall MAF less than 10%
  my $temp_string = join('', $mstring, $sstring, $ustring, $rstring);
  $temp_string =~ s/A/AA/ig;
  $temp_string =~ s/C/CC/ig;
  $temp_string =~ s/G/GG/ig;
  $temp_string =~ s/T/TT/ig;
  
  $temp_string =~ s/M/AC/ig;
  $temp_string =~ s/R/AG/ig;
  $temp_string =~ s/W/AT/ig;
  $temp_string =~ s/S/CG/ig;
  $temp_string =~ s/Y/CT/ig;
  $temp_string =~ s/K/GT/ig;

  my $mac = 0;
  for (my $k = 0; $k < length $temp_string; $k++) {
    if (substr($temp_string, $k, 1) =~ /$q/) {
      $mac++;
    }
  }
  my $maf = $mac/(length $temp_string);
  if ($maf < 0.1) {
    next LOOP;
  }
  
  # Check that all (required) individuals are present
  my $check = 0;
  my @m =split('', $mstring);
  for (my $i = 0; $i < scalar @m; $i++) {
    if ($m[$i] =~ /N/) {
      unless (($i == 14) || ($i == 15) || ($i == 16)) {
	$check++;
      }
    }
  }
  if ($check > 0) {
    next LOOP;
  }
  my @s =split('', $sstring);
  for (my $i = 0; $i < scalar @s; $i++) {
    if ($s[$i] =~ /N/) {
      $check++;
    }
  }
  if ($check > 0) {
    next LOOP;
  }
  my @u =split('', $ustring);
  for (my $i = 0; $i < scalar @u; $i++) {
    if ($u[$i] =~ /N/) {
      unless (($i == 17) || ($i == 20) || ($i == 23)) {
	$check++;
      }
    }
  }
  if ($check > 0) {
    next LOOP;
  }
  my @r =split('', $rstring);
  for (my $i = 0; $i < scalar @r; $i++) {
    if ($r[$i] =~ /N/) {
      unless ($i == 19) {
	$check++;
      }
    }
  }
  if ($check > 0) {
    next LOOP;
  }

  # Now, check if the position is on the same tag as the previous position
  my ($t, $pos) = split(/\./, $id);
  if ($t == $tag) {
    push (@lines, $_);
  } else {
    unless ((scalar @lines) == 0) {

      # Process the information into migrate format output
      for (my $i = 0; $i < scalar @Mak; $i++) {
	if (($i == 14) || ($i == 15) || ($i == 16)) {
	  next;
	} else {
	  my $gen1 = '';
	  my $gen2 = '';
	  foreach my $line (@lines) {
	    my @info = split(/\t/, $line);
	    my $base = substr($info[3], $i, 1);
	    if ($base =~ /[MRWSYK]/) {
	      $gen1 .= $info[1];
	      $gen2 .= $info[2];
	    } else {
	      $gen1 .= $base;
	      $gen2 .= $base;
	    }
	  }
	  print MAK $Mak[$i], "_1", " ", " ", $gen1, "\n";
	  print MAK $Mak[$i], "_2", " ", " ", $gen2, "\n";
	}
      }

      for (my $i = 0; $i < scalar @SA; $i++) {
	my $gen1 = '';
	my $gen2 = '';
	foreach my $line (@lines) {
	  my @info = split(/\t/, $line);
	  my $base = substr($info[4], $i, 1);
	  if ($base =~ /[MRWSYK]/) {
	    $gen1 .= $info[1];
	    $gen2 .= $info[2];
	  } else {
	    $gen1 .= $base;
	    $gen2 .= $base;
	  }
	}
	print SA $SA[$i], "_1", " ", " ", $gen1, "\n";
	print SA $SA[$i], "_2", " ", " ", $gen2, "\n";
      }

      for (my $i = 0; $i < scalar @Ugi; $i++) {
	if (($i == 17) || ($i == 20) || ($i == 23)) {
	  next;
	} else {
	  my $gen1 = '';
	  my $gen2 = '';
	  foreach my $line (@lines) {
	    my @info = split(/\t/, $line);
	    my $base = substr($info[5], $i, 1);
	    if ($base =~ /[MRWSYK]/) {
	      $gen1 .= $info[1];
	      $gen2 .= $info[2];
	    } else {
	      $gen1 .= $base;
	      $gen2 .= $base;
	    }
	  }
	  print UGI $Ugi[$i], "_1", " ", " ", $gen1, "\n";
	  print UGI $Ugi[$i], "_2", " ", " ", $gen2, "\n";
	}
      }

      for (my $i = 0; $i < scalar @Rus; $i++) {
	if ($i == 19) {
	  next;
	} else {
	  my $gen1 = '';
	  my $gen2 = '';
	  foreach my $line (@lines) {
	    my @info = split(/\t/, $line);
	    my $base = substr($info[6], $i, 1);
	    if ($base =~ /[MRWSYK]/) {
	      $gen1 .= $info[1];
	      $gen2 .= $info[2];
	    } else {
	      $gen1 .= $base;
	      $gen2 .= $base;
	    }
	  }
	  print RUS $Rus[$i], "_1", " ", " ", $gen1, "\n";
	  print RUS $Rus[$i], "_2", " ", " ", $gen2, "\n";
	}
      }
    }
    $tag = $t;
    @lines = ();
    push (@lines, $_);
  }
}
close(IN);

# Process the last set of tags
for (my $i = 0; $i < scalar @Mak; $i++) {
  if (($i == 14) || ($i == 15) || ($i == 16)) {
    next;
  } else {
    my $gen1 = '';
    my $gen2 = '';
    foreach my $line (@lines) {
      my @info = split(/\t/, $line);
      my $base = substr($info[3], $i, 1);
      if ($base =~ /[MRWSYK]/) {
	$gen1 .= $info[1];
	$gen2 .= $info[2];
      } else {
	$gen1 .= $base;
	$gen2 .= $base;
      }
    }
    print MAK $Mak[$i], "_1", " ", " ", $gen1, "\n";
    print MAK $Mak[$i], "_2", " ", " ", $gen2, "\n";
  }
}

for (my $i = 0; $i < scalar @SA; $i++) {
  my $gen1 = '';
  my $gen2 = '';
  foreach my $line (@lines) {
    my @info = split(/\t/, $line);
    my $base = substr($info[4], $i, 1);
    if ($base =~ /[MRWSYK]/) {
      $gen1 .= $info[1];
      $gen2 .= $info[2];
    } else {
      $gen1 .= $base;
      $gen2 .= $base;
    }
  }
  print SA $SA[$i], "_1", " ", " ", $gen1, "\n";
  print SA $SA[$i], "_2", " ", " ", $gen2, "\n";
}

for (my $i = 0; $i < scalar @Ugi; $i++) {
  if (($i == 17) || ($i == 20) || ($i == 23)) {
    next;
  } else {
    my $gen1 = '';
    my $gen2 = '';
    foreach my $line (@lines) {
      my @info = split(/\t/, $line);
      my $base = substr($info[5], $i, 1);
      if ($base =~ /[MRWSYK]/) {
	$gen1 .= $info[1];
	$gen2 .= $info[2];
      } else {
	$gen1 .= $base;
	$gen2 .= $base;
      }
    }
    print UGI $Ugi[$i], "_1", " ", " ", $gen1, "\n";
    print UGI $Ugi[$i], "_2", " ", " ", $gen2, "\n";
  }
}

for (my $i = 0; $i < scalar @Rus; $i++) {
  if ($i == 19) {
    next;
  } else {
    my $gen1 = '';
    my $gen2 = '';
    foreach my $line (@lines) {
      my @info = split(/\t/, $line);
      my $base = substr($info[6], $i, 1);
      if ($base =~ /[MRWSYK]/) {
	$gen1 .= $info[1];
	$gen2 .= $info[2];
      } else {
	$gen1 .= $base;
	$gen2 .= $base;
      }
    }
    print RUS $Rus[$i], "_1", " ", " ", $gen1, "\n";
    print RUS $Rus[$i], "_2", " ", " ", $gen2, "\n";
  }
}

close(MAK);
close(SA);
close(UGI);
close(RUS);

# Need to count the number of loci and the length of each locus
my $first = $Mak[0];
$first .= '_1';
my $num_loci = 0;
my @lengths = ();

open (MAK, 'makira.migrate') || die "\nUnable to open the file makira.migrate!\n";
while (<MAK>) {
  chomp $_;
  my ($name, $locus) = split(/\s{1,}/, $_);
  if ($name =~ /^$first/) {
    $num_loci++;
    my $len = length $locus;
    push (@lengths, $len);
  } else {
    next;
  }
}
close(MAK);

# Open the output file to print the full data set
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

print OUT "N", "\t", "4\t", $num_loci, "\n";
print OUT "@lengths\n";

print OUT "42\t", "Mak\n";
open (MAK, 'makira.migrate') || die "\nUnable to open the file makira.migrate!\n";
while (<MAK>) {
  print OUT $_;
}
close(MAK);

print OUT "42\t", "SA\n";
open (SA, 'santaana.migrate') || die "\nUnable to open the file santaana.migrate!\n";
while (<SA>) {
  print OUT $_;
}
close(SA);

print OUT "42\t", "Ugi\n";
open (UGI, 'ugi.migrate') || die "\nUnable to open the file ugi.migrate!\n";
while (<UGI>) {
  print OUT $_;
}
close(UGI);

print OUT "42\t", "Rus\n";
open (RUS, 'russells.migrate') || die "\nUnable to open the file russells.migrate!\n";
while (<RUS>) {
  print OUT $_;
}
close(RUS);
close(OUT);

exit;
