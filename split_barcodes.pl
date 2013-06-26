#!/usr/bin/perl -w

# Note that this script is pretty slow -- just to read such a large file takes awhile.
# It would be hard to get much faster in perl, as I've already used a lot of efficiency tricks to
# make it faster (e.g. using a hash of file handles)
#
# Consider rewriting in C sometime if this is likely to be used often.....

use strict;
use Getopt::Long qw(:config no_ignore_case);

my $SCARF = 1;
my $SCARF_ASCII = 2;

my ($barcode_file, $sequence_file, $barcode_stats_flag, $trim_flag, $compress_output);

GetOptions(
	"b|barcodes=s"				=> \$barcode_file,
	"c|compress"				=> \$compress_output,
	"s|seqs=s"					=> \$sequence_file,
	"r|report"					=> \$barcode_stats_flag,
	"t|trim"					=> \$trim_flag
) or die "Usage: split_barcodes.pl [-b barcode file] [-s sequence file] [-c (compress outputs using gzip)] [-r (report barcode stats)] [-t (trim barcodes off after demultiplexing)]\n";

my $barcode_out = $sequence_file . ".barcode_stats" if $barcode_stats_flag;

die "Usage: split_barcodes.pl [-b barcode file] [-s sequence file] [-c (compress outputs using gzip)] [-r (report barcode stats)] [-t (trim barcodes off after demultiplexing)]\n" unless $barcode_file && $sequence_file;

my ($barcodes, $barcode_to_name, $files_created) = read_barcodes($barcode_file);
split_sequences($barcodes, $sequence_file, $barcode_out, $barcode_to_name, $files_created);

exit;

sub split_sequences {
	my ($barcodes, $in, $barcode_out, $barcode_to_name, $paths_by_bc)=@_;
	
	my ($num_read, $num_missed) = 0;
	my %counts;

	my $barcode_len = get_barcode_len($barcodes);

	# Make best guess about what type of sequence file is being split
	my $seq_format = infer_seq_format($sequence_file);
	my %files_created;	# Used to keep track of all files created so they can be compressed later
	
	if ($seq_format eq 'scarf') {
		my $scarf_type = 0;
		open(IN, $in) || die "Can't open input file $in!\n";
		while (<IN>) {
			chomp;
			$scarf_type = infer_scarf_type($_) unless $scarf_type;	

			my $missed_barcode = 1;
			my @pieces = split /:/, $_;
			my $bc = substr($pieces[5], 0, $barcode_len);
			if (my $fh = $barcodes->{$bc}) { # match barcode
				$files_created{$paths_by_bc->{$bc}} = 1;
				$missed_barcode = 0;
				$counts{$bc}++;
				# Do trimming if -t flag was used; otherwise, do not modify sequence or quality scores
				if ($trim_flag) {
					$pieces[5] = substr($pieces[5], $barcode_len);
					if ($scarf_type == $SCARF_ASCII) {
						$pieces[6]=~ s/^\s*\D{$barcode_len}//;
					}
					elsif ($scarf_type == $SCARF) {
						$pieces[6] =~ s/^\s*([-\d]+\s){$barcode_len}//
					}
				}
				my $seq = join ":", @pieces;
				print $fh "$seq\n";
			}

			$num_missed++ if $missed_barcode;
			$num_read++;
		}
		close IN;
	}
	# Note: FASTQ splitting is under development and not yet rigorously tested
	elsif ($seq_format eq 'fastq') {
		open(IN, $in) || die "Can't open input file $in!\n";
		while(<IN>) {
			my $line1 = $_;		# Header line
			my $line2 = <IN>;	# Sequence line
			my $line3 = <IN>;	# "+" line
			my $line4 = <IN>;	# Quality scores line
			my $missed_barcode = 1;
			my $bc = substr($line2, 0, $barcode_len);
			if (my $fh = $barcodes->{$bc}) {
				$files_created{$paths_by_bc->{$bc}} = 1;
				$missed_barcode = 0;
				$counts{$bc}++;
				if ($trim_flag) {
					$line2 = substr($line2, $barcode_len);
					$line4 = substr($line4, $barcode_len);
				}
				print $fh $line1.$line2.$line3.$line4;
			}
			$num_missed++ if $missed_barcode;
			$num_read++;
		}
		close IN;
	}
	else {
		die "Could not split file (neither FASTQ nor SCARF format detected for your input file).\n";
	}
	my $num_matched = $num_read-$num_missed;
	printf STDERR "\n\t\tSUMMARY STATS\n";
	printf STDERR "%d of %d reads (%.1f%%) had a barcode match\n", $num_matched, $num_read, ($num_matched/$num_read)*100;
	printf STDERR "#%d\t%d\t%f\n", $num_matched, $num_read, ($num_matched/$num_read);
	if ($barcode_out) {
		open (OUT, ">$barcode_out") or die "can't open $barcode_out: $!\n";
	}	
	printf OUT "%d\t%d\t%f\n", $num_matched, $num_read, ($num_matched/$num_read) if $barcode_out;
	
	my @s = sort {$counts{$b} <=> $counts{$a}} keys %counts;
	for my $s (@s) {
		printf STDERR "$s\t%.1f%% ($counts{$s} of $num_matched)\n", ($counts{$s}/$num_matched) * 100;
		if ($barcode_to_name && $barcode_to_name->{$s}) {
			printf OUT "$barcode_to_name->{$s}\t$s\t$counts{$s}\t$num_matched\n", ($counts{$s}/$num_matched) * 100 if $barcode_out;
		}
		else {
			printf OUT "$s\t$counts{$s}\t$num_matched\n", ($counts{$s}/$num_matched) * 100 if $barcode_out;
		}
	}
	close OUT if $barcode_out;
	
	# Compress outputs if -c flag turned on
	if ($compress_output) {
		foreach my $k (keys %files_created) {
			print "Compressing file $k...\n";
			`gzip $k`;
			print "Done.\n";
		}
	}
}

# This is a clunky, ugly piece of code, but it should work as a hack for now
# To be detected as FASTQ formatted, a file must:
#	1) Have a first line that starts with "@" (as the header lines for FASTQ sequences do)
#	2) Have a second line with some sequence of one or more A/T/C/G/N bases
#	3) Have a third line that consists only of the character '+'
# To be detected as SCARF formatted, a file must:
#	1) Have a first line that consists of 7 terms of any character type, separated by colons (i.e., 6 colons total)
sub infer_seq_format {
	my $filename = shift;
	open(TEST, $filename) || die "Couldn't open input file $filename when attempting to detect file type!\n";
	my $firstline = <TEST>;
	if ($firstline =~ m/^@/) {
		my $secondline = <TEST>;
		if ($secondline =~ /[ATCGN]+/) {
			my $thirdline = <TEST>;
			if ($thirdline =~ m/^\+$/) {
				return 'fastq';		
			}
		}
	}
	elsif ($firstline =~ /(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)/) {
		return 'scarf';
	}
	else {
		die "Your input file $filename does not appear to conform to either the SCARF or FASTQ formats. This script only supports these two file types.\n";
	}
	close TEST;
}

sub get_barcode_len {
	my $barcodes = shift;

	my @keys = keys %$barcodes;
	my $first = shift @keys;

	my $len = length($first);

	for my $k (@keys) {
		if ($len != length($k)) {
			die "Error: barcodes must all be the same length (found $len and " . length($k)  . " bp\n";
		}
	}

	return $len;
}

sub infer_scarf_type {
	shift;

	if(/^\S+:\d+:\d+:\d+:\d+.?.?.?.?:\D+$/) {
		return $SCARF_ASCII;
	}
        elsif(/^\S+:\d+:\d+:\d+:\d+/) {
		return $SCARF;
	}

}

sub read_barcodes {
	my $in=shift;	# $in = barcode mapping file path
	open (IN, $in) or die "can't open $in: $!\n";

	my %barcodes;
	my %seq_to_barcode_name;
	my %file_paths_by_barcode;
	while (<IN>) {
		chomp;
		my (@stuff) = split /\t/;

		my ($barcode_name, $barcode, $filename);
		if (scalar @stuff == 3) {
			($barcode_name, $barcode, $filename) = @stuff;
		}
		else {
			($barcode, $filename) = @stuff;
		}
	 	open my $out, ">$filename" or die "can't open $filename: $!\n";
		$file_paths_by_barcode{$barcode} = $filename;
		$barcodes{$barcode} = $out; # barcode points to the file handle
		if ($barcode_name) {
			$seq_to_barcode_name{$barcode} = $barcode_name;
		}
	}

	close IN;
	
	return \%barcodes, \%seq_to_barcode_name, \%file_paths_by_barcode;
}
