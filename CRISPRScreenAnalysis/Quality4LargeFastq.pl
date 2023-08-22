#!/usr/bin/perl -w
use warnings;
use strict;
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Linux: This script is called in Quality4LargeFastq.R
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
#####################################################################################################
#                                            SETTINGS

# Barcodes (all of the same length)
my @barcode = ("CGTGAT", "ACATCG", "GCCTAA", "TGGTCA", "CACTGT", "ATTGGC", 
		"GATCTG", "TCAAGT", "CTGATC", "AAGCTA", "GTAGCC", "TACAAG");
		
# Number of allowed mismaches in barcodes (Use the same number as you used in FastqToCountTable.pl)
my $BCmm = 1;
#####################################################################################################
# Make hash with barcode numbers
my %BCnum;
foreach my $i (1 .. scalar @barcode) {
  $BCnum{$barcode[$i-1]} = $i;
}
# Make hash with barcode pieces 
my %BCpieces;
my $BCsize = length $barcode[0];
my $index = 1;
foreach my $bc (@barcode) {
  my $bcn = $index++;
  foreach my $nt (1 .. $BCsize) {
    $BCpieces{"$bcn"."."."$nt"} = substr($bc,$nt-1,1);   
  }  
}

while (my $filename = shift) {
  # Make multidimensional array with read sequences
  my @readTable = ();
  my @allSeqs = ();
  my @countTable = ();
  # Read Fastq file and write reads in the right barcode file
  open FASTQ, "cat $filename |" or die "Could not open sequence file $_";
  while (<FASTQ>) { 
    my $seq = <FASTQ>; 
    chomp $seq;
    <FASTQ>;
    <FASTQ>;
    
    # Identify barcode
    my $BCnumber = 0;
    my $BCreads = substr($seq,0,$BCsize);
    if (grep(/^$BCreads$/, @barcode)) {
      # When barcode has no mismatches
      $BCnumber = $BCnum{$BCreads};
    } else {if ($BCmm > 0) {       
                # Cut read BC in pieces
                my %BCread = ();
                foreach my $i (1 .. $BCsize) {
                  $BCread{"s".$i} = substr($seq,$i-1,1);
                }
                # Calculate similarities with barcode hash
                foreach my $ibc (1 .. scalar @barcode) {
                  my $m1 = 0;
                  foreach my $i (1 .. $BCsize) {
                    if (($BCread{"s".$i}) eq ($BCpieces{"$ibc"."."."$i"})) {$m1+=1};
                  }
                  if ($m1 >= ($BCsize-$BCmm)) {$BCnumber = $ibc};
                  last if ($m1 >= ($BCsize-$BCmm));
                }
            }
      }
      # Make multidimensional array with read sequences     
      push @{$readTable[$BCnumber]},substr($seq,$BCsize,(length $seq)-($BCsize));
      push @allSeqs,substr($seq,$BCsize,(length $seq)-($BCsize));
  }  
  close FASTQ; 

  # Take unique reads
  my @uniqueReads = "Read_without_barcode";
  my %uniqueSeqs = ();
  foreach my $read (@allSeqs) {
    ++$uniqueSeqs{$read};
  }
  $uniqueSeqs{$_}++ foreach @allSeqs;
  foreach my $keys (sort keys %uniqueSeqs) {
    push(@uniqueReads,$keys);
  }
  push @{$countTable[0]},@uniqueReads;
  
  # Produce count tables per barcode
  foreach my $i (0 .. scalar @barcode) {
    my @reads = @{$readTable[$i]};
    # Calculate counts per unique read
    my %ReadCounts = ();
      foreach my $read (@reads) {
         ++$ReadCounts{$read};
      }         
    # Add to count table
    push @{$countTable[$i+1]}, @ReadCounts{@uniqueReads};
  }
  
  # Name barcoded samples
  $countTable[1][0] = "Unidentified";
  foreach my $i (1 .. scalar @barcode) {
    $countTable[$i+1][0] = $barcode[$i-1];
  }
  # Replace undefined elements for 0 counts
  my $numberguides = @uniqueReads;
  foreach my $i (0 .. scalar @barcode) {
    for (my $a=0; $a<$numberguides; $a++) {
      if (!defined($countTable[$i+1][$a])) {
          $countTable[$i+1][$a] = 0;
      }
    }
  }
  # Transpose table
  my @transposed = ();
  for my $row (@countTable) {
    for my $column (0 .. $#{$row}) {
      push(@{$transposed[$column]}, $row->[$column]);
    }
  }
  # Write count table file
  open RT, ">ReadTable\_".$filename.".csv" or die "Could not open read table file $_";
  for my $new_row (@transposed) {
    for my $new_col (@{$new_row}) {
        print RT $new_col, ",";
    }
    print RT "\n";
  }
  close RT;
}
