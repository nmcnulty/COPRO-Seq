#!/usr/bin/perl

# Usage: perl get_full_length_NM_seqs.pl [NM FASTA] [original SCARF] [barcode length]

# Take as input the "no match" FASTA file produced by the COPRO-Seq pipeline, as well as the
#	SCARF file from which it was derived, and return the full-length sequences for use in a 
#	BLAST search
# This script is necessary because the trimmed Illumina reads are usually too short to give
#	significant hits against large databases like nt

use strict;
use warnings;

# Get list of relevant headers that can be passed to grep via the -f option
open(INPUT, $ARGV[0]) || die "Can't open input FASTA file $ARGV[0]!\n";
open(PATTERNS, ">patterns.txt") || die "Can't open output file patterns.txt!\n";
while(<INPUT>) {
	chomp;
	my $a = $_;
	$a =~ s/>//;
	if(/^>/) {
		print PATTERNS $a."\n";
	}
}
close INPUT;
close PATTERNS;

# Grab all relevant lines from original SCARF file and pass to temporary text file
system("grep -F -f patterns.txt $ARGV[1] > tmp.txt");

# Cycle through each line of temporary text file, trimming barcode from sequence
open(TMP, "tmp.txt") || die "Can't open temporary file tmp.txt!\n";
open(OUT, ">out.txt") || die "Can't open output file out.txt!\n";
while(<TMP>) {
	my @line = split(/\t/, $_);
	my $seq = substr($line[2], $ARGV[2]);
	print OUT "$line[0]\n$seq\n";
}
close TMP;
close OUT;

exit;