#!/usr/bin/perl -w
use warnings;
use strict;
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Nextera PCR primers can be found in 'PCR design Nextera-NovaSeq.xlsx'
# Linux: run in command line: ./FastqToCountTable.pl LibraryX.csv *.fastq.gz
# Windows: Install Strawberry Perl, Unpack fastq.gz file with 7zip, replace 'zcat' (line 59) for 'type',
# and run in command prompt: perl FastqToCountTable.pl LibraryX.csv data1.fastq data2.fastq
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
#####################################################################################################
#                                            SETTINGS

# Barcodes (all the same lentgh) and allowed barcode mismaches and max indels allowed in upstream region (primers)
# If a barcode is not present in your FASTQ file it could give an error
my @barcode = ("CGTGAT", "ACATCG", "GCCTAA", "TGGTCA", "CACTGT", "ATTGGC", 
                "GATCTG", "TCAAGT", "CTGATC", "AAGCTA", "GTAGCC", "TACAAG");
my ($BCmm, $indel) = (1, 3);

# Location CRISPR seq-1, length CRISPR, and upstream nucleotides (lentiCRISPRv2/LentiGuide: 41, 20, "CACCG")
my ($loc, $CRISPRsize, $upseq)  = (41, 20, "CACCG");
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

# Load library and make count table precursor
my $library = shift;
my (@libraryGuides, @libraryGenes, @librarySeqs);
open (my $lib, '<', $library) or die "Could not open library file $_";
while (my $line = <$lib>) {
  chomp $line;
  my @fields = split ",", $line;
  push @libraryGuides, $fields[0];
  push @libraryGenes, $fields[1];
  push @librarySeqs, $fields[2];
}
close $lib;
# Array of arrays (multidimensional array)
my @countTablePre = ([@libraryGuides], [@libraryGenes], [@librarySeqs]);

while (my $filename = shift) {
  # Make multidimensional array with read sequences 
  my @readTable = ();
  
  # Read Fastq file and write reads in the right barcode file
  open FASTQ, "zcat $filename |" or die "Could not open sequence file $_";
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
    
    # Identify CRISPR sequence
    if ($BCnumber > 0) {
      my $guide = substr($seq,$loc,$CRISPRsize);
      # When there is an indel in constant region, change CRISPR seq location 
      my $upseqlength = length $upseq;
      my $seqcheck = substr($seq,$loc-$upseqlength,$upseqlength);
      if ($upseq ne $seqcheck && $indel > 0) {
        foreach my $i (-$indel .. $indel) {
           my $seqcheck = substr($seq,$loc-$upseqlength+$i,$upseqlength);
           if ($upseq eq $seqcheck) {$guide = substr($seq,$loc+$i,$CRISPRsize)};
           last if ($upseq eq $seqcheck);
         }
       }
      # Make multidimensional array with read sequences     
      push @{$readTable[$BCnumber-1]},$guide;
    }
    
  }
  close FASTQ;
    
  # Start new count table
  my @countTable = @countTablePre;
  
  # Produce count tables per barcode
  foreach my $i (1 .. scalar @barcode) {
    my @reads = @{$readTable[$i-1]};
    # Calculate counts per unique read
    my %readCounts = ();
      foreach my $read (@reads) {
         ++$readCounts{$read};
      }            
    # Add to count table
    push @{$countTable[$i+2]},@readCounts{@librarySeqs};
  }
  
  # Name barcoded samples
  foreach my $i (1 .. scalar @barcode) {
    $countTable[$i+2][0] = $barcode[$i-1];
  }
  # Replace undefined elements for 0 counts
  my $numberguides = @librarySeqs;
  foreach my $i (1 .. scalar @barcode) {
    for (my $a=0; $a<$numberguides; $a++) {
      if (!defined($countTable[$i+2][$a])) {
          $countTable[$i+2][$a] = 0;
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
  open CT, ">CountTable\_".$filename.".csv" or die "Could not open count table file $_";
  for my $new_row (@transposed) {
    for my $new_col (@{$new_row}) {
        print CT $new_col, ",";
    }
    print CT "\n";
  }
  close CT;
}
