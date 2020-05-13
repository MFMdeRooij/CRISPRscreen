#!/usr/bin/perl -w
use warnings;
use strict;
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Linux: run in command line: ./Quality.pl data1.fastq.gz data2.fastq.gz
# Windows: Install Strawberry Perl, Unpack fastq.gz file with 7zip, replace 'zcat' (line 45) for 'type',
# and run in command prompt: perl Quality.pl data1.fastq data2.fastq
# After running this code (and FastqToCountTable.pl), go further with the R script (Quality.R)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
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
my $index = 1;
foreach my $i (@barcode) {
  my $bcn = $index++;
  $BCnum{$i} = $bcn;
}
# Make hash with barcode pieces 
my %BCpieces;
my $BCsize = length $barcode[0];
$index = 1;
foreach my $i (@barcode) {
  my $bcn = $index++;
  foreach my $n (1 .. $BCsize) {
    $BCpieces{"$bcn"."."."$n"} = substr($i,$n-1,1);   
  }  
}

while (my $filename = shift) {
  # Make multidimensional array with read sequences
  my @readTable = ();
  my @allSeqs = ();
  my @countTable = ();
  # Read Fastq file and write reads in the right barcode file
  open FASTQ, "zcat $filename |" or die "Could not open sequence file $_";
  while (my $name = <FASTQ>) { 
    my $seq = <IN>; 
    chomp $seq;
    my $plusline = <IN>;
    my $qual = <IN>;
    
    # Determine barcode
    my $BCnumber = 0;
    my $BCreads = substr($seq,0,$BCsize);
    if (grep(/^$BCreads$/, @barcode)) {
      # When barcode has no mismatches
      $BCnumber = $BCnum{$BCreads};
    } else {  if ($BCmm > 0) {       
                # Cut read BC in pieces
                my %BCread = ();
                foreach my $i (1 .. $BCsize) {
                  $BCread{"s".$i} = substr($seq,$i-1,1);
                }
                # Calculate similarities with barcode hash
                $index = 1;
                foreach my $ibc (@barcode) {
                  my $n = $index++;     
                  my $m1 = 0;
                  foreach my $i (1 .. $BCsize) {
                    if (($BCread{"s".$i}) eq ($BCpieces{"$n"."."."$i"})) {$m1 = $m1+1};
                  }
                  if ($m1 >= ($BCsize-$BCmm)) {$BCnumber = $n};
                  last if ($m1 >= ($BCsize-$BCmm));
                }
              }
      }
      # Make multidimensional array with read sequences     
      push @{$readTable[$BCnumber]},substr($seq,$BCsize,(length $seq)-$BCsize);
      push @allSeqs,substr($seq,$BCsize,(length $seq)-$BCsize);
  }  
  close IN; 

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
  foreach my $i (0 .. scalar @barcode) {
    $countTable[$i+1][0] = "BC".$i;
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
  open FILE, ">ReadTable\_".$filename.".csv" or die "Could not open read table file $_";
  for my $new_row (@transposed) {
    for my $new_col (@{$new_row}) {
        print FILE $new_col, ",";
    }
    print FILE "\n";
  }
  close FILE;
}
