#!/usr/bin/perl

# Inputs: Sequence file (SCARF, FASTA or FASTQ format), a text file of sequence headers to include/exclude in output, output filename, barcode length, mode, sequence file format

# This script looks for a list of sequences as defined by a text file of headers and parses an original sequencing results file,
#	either reporting all full-length sequences matching the headers provided ('include' mode) or all full-length sequences
#	not matching the headers provided ('exclude' mode).

# It is often used in conjunction with the COPRO-Seq pipeline to get a FASTA file of full-length sequences for reads
#	that failed to match any of the reference genomes in an analysis (i.e., those assigned the 'NM' code)

# Usage: perl get_full_length_NM_seqs.pl [text file of headers] [original SCARF/FASTQ sequencing results] [output file] [barcode length] [mode ('i' for include or 'e' for 'exclude') [sequence file format ('fastq' or 'scarf')]

use strict;
use warnings;
use Getopt::Long;

my ($headers_file, $seq_file, $output_file, $bc_length, $mode, $seq_format);

GetOptions (	'b|bclength=i'			=> \$bc_length,
				'f|format=s'			=> \$seq_format,	# 'fastq', 'fasta' or 'scarf'
				'h|headers=s'			=> \$headers_file,
				'm|mode=s'				=> \$mode,			# 'i' (include mode) or 'e' (exclude mode)
				'o|output=s'			=> \$output_file,
				's|seq_results=s'		=> \$seq_file
) or die usage();

# Check a few parameters (eventually consolidate these lines into a subroutine)
if (defined ($seq_format)) {
	unless (($seq_format eq 'scarf') || ($seq_format eq 'fastq') || ($seq_format eq 'fasta')) { die "Sequencing file format must be specified as 'scarf', 'fasta' or 'fastq'.\n\n"; }	
}
else { die "You must define the sequence file format using -f ('scarf', 'fasta' or 'fastq').\n\n"; }

unless (($mode eq 'i') || ($mode eq 'e')) {	die "Mode must be specified as either 'e' (exclude entries from results) or 'i' (include entries in results).\n\n"; }
unless (-e $seq_file) { die "Can't find sequence file $seq_file.\n\n"; }
unless (-e $headers_file) { die "Can't find text file of sequences to search for named $headers_file.\n\n"; }

my %unique_headers;	# Will store a unique list of all sequence headers from the headers text file provided
# Cycle through headers
open(HEADERS, $headers_file) || die "Can't open headers text file $headers_file!\n";
print "Opening headers file...\n";
while(<HEADERS>) {
	chomp;
	$unique_headers{$_} = 1;
}
print "Done grabbing headers.\n";
close HEADERS;

open(SEQ, $seq_file) || die "Can't open input (original sequences) file $seq_file!\n";
open(OUTPUT, ">$output_file") || die "Can't open output file $output_file!\n";
#	If you find a header that was in your headers file (using the hash as a lookup):
#		If you're in 'include' mode: report the header (add ">") and 2nd line (sequence), after trimming barcode (if length defined), to your output file, then ignore 3rd and 4th lines
#		If you're in 'exclude' mode: ignore this line and the next three lines without doing anything
#	If you find a header that wasn't in your headers file:
#		If you're in 'include' mode: ignore this line and the next three lines without doing anything
#		If you're in 'exclude' mode: report the header (add ">") and 2nd line (sequence), after trimming barcode (if length defined), to your output file, then ignore 3rd and 4th lines
my $trackprogress = 0;
while (<SEQ>) {
	$trackprogress++;
	if ($trackprogress % 10000 == 0) { print "Done searching $trackprogress entries in your original sequencing results...\n"; }
	if ($seq_format eq 'fasta') {
		if ($_ =~ /^>(.+)/) {
			my $h = $1;
			$h =~ s/\s/_/;	# .NM FASTAs reported by COPRO-Seq pipeline have had all whitespace in original results converted to underscores
			if (defined $unique_headers{$h}) {
				if ($mode eq 'i') {
					my $untrimmedseq = <SEQ>;
					print OUTPUT ">$h\n";
					if (defined $bc_length) {
						print OUTPUT substr($untrimmedseq, $bc_length);
					}
					else {
						print OUTPUT $untrimmedseq;
					}
				}
				elsif ($mode eq 'e') {
					# Ignore next line (sequence)
					<SEQ>;
				}
				else {
					die "Your mode (-m) appears to be set incorrectly.\n";
				}
			}
			else {
				if ($mode eq 'i') {
					# Ignore next line (sequence)
					<SEQ>;
				}
				elsif ($mode eq 'e') {
					my $untrimmedseq = <SEQ>;
					print OUTPUT ">$h\n";
					if (defined $bc_length) {
						print OUTPUT substr($untrimmedseq, $bc_length);
					}
					else {
						print OUTPUT $untrimmedseq;
					}
				}
				else {
					die "Your mode (-m) appears to be set incorrectly.\n";
				}
			}
		}
		else { die	"There seems to be a problem with the format of your FASTA sequence file or with this script's efforts to parse it.\n"; }
	}
	if ($seq_format eq 'fastq') {
		if ($_ =~ /^@(.+)/) {
			my $h = $1;
			$h =~ s/\s/_/;	# .NM FASTAs reported by COPRO-Seq pipeline have had all whitespace in original results converted to underscores
			if (defined $unique_headers{$h}) {
				if ($mode eq 'i') {
					my $untrimmedseq = <SEQ>;
					print OUTPUT ">$h\n";
					if (defined $bc_length) {
						print OUTPUT substr($untrimmedseq, $bc_length);
					}
					else {
						print OUTPUT $untrimmedseq;
					}
					<SEQ>;
					<SEQ>;
				}
				elsif ($mode eq 'e') {
					# Ignore next three lines
					<SEQ>;
					<SEQ>;
					<SEQ>;
				}
				else { die "Your mode (-m) appears to be set incorrectly.\n"; }
			}
			else {
				if ($mode eq 'i') {
					# Ignore next three lines 
					<SEQ>;
					<SEQ>;
					<SEQ>;
				}
				elsif ($mode eq 'e') {
					my $untrimmedseq = <SEQ>;
					print OUTPUT ">$h\n";
					if (defined $bc_length) {
						print OUTPUT substr($untrimmedseq, $bc_length);
					}
					else {
						print OUTPUT $untrimmedseq;
					}
					<SEQ>;
					<SEQ>; 					
				}
				else { die "Your mode (-m) appears to be set incorrectly.\n"; }
			}
		}
		else { die	"There seems to be a problem with the format of your FASTA sequence file or with this script's efforts to parse it.\n"; }
	}
	if ($seq_format eq 'scarf') {
		# Don't yet have any validation of the sequence file format for SCARF -- should add at some point
		my @line = split(/:/, $_);
		my $h = "$line[0]:$line[1]:$line[2]:$line[3]:$line[4]";
		$h =~ s/\s/_/;	# .NM FASTAs reported by COPRO-Seq pipeline have had all whitespace in original results converted to underscores
		if (defined $unique_headers{$h}) {
			if ($mode eq 'i') {
				my $untrimmedseq = $line[5];
				print OUTPUT ">$h\n";
				if (defined $bc_length) {
					print OUTPUT substr($untrimmedseq, $bc_length);
				}
				else {
					print OUTPUT $untrimmedseq;
				}
			}
			elsif ($mode eq 'e') {
				# Do nothing (ignore line)
			}
			else { die "Your mode (-m) appears to be set incorrectly.\n"; }
		}
		else {
			if ($mode eq 'i') {
				# Do nothing (ignore line)
			}
			elsif ($mode eq 'e') {
				my $untrimmedseq = $line[5];
				print OUTPUT ">$h\n";
				if (defined $bc_length) {
					print OUTPUT substr($untrimmedseq, $bc_length);
				}
				else {
					print OUTPUT $untrimmedseq;
				}
			}
			else { die "Your mode (-m) appears to be set incorrectly.\n"; }
		}			
	}
}

close OUTPUT;
close SEQ;

exit;