#!/usr/bin/perl -w

# note this is pretty slow; just to read such a large file in perl takes a while
# would be hard to gain much more performance in perl, as I've already used a lot of 
# efficiency tricks to make it faster (e.g. use a hash of file handles)
#
# might be worth rewriting in C sometime.....

use strict;
use Getopt::Long qw(:config no_ignore_case);

my $SCARF = 1;
my $SCARF_ASCII = 2;

my ($barcode_file, $solexa_file, $barcode_stats_flag, $trim_flag);

GetOptions(
	"b|barcodefile=s"			=> \$barcode_file,
	"s|solexafile=s"			=> \$solexa_file,
	"r|report"					=> \$barcode_stats_flag,
	"t|trim"					=> \$trim_flag
) or die "Usage: split_barcodes.pl [-b barcode file] [-s sequence file] [-r (report barcode stats)] [-t (trim barcodes off after demultiplexing)]\n";

my $barcode_out = $solexa_file . ".barcode_stats" if $barcode_stats_flag;

die "Usage: split_barcodes.pl [-b barcode file] [-s sequence file] [-r (report barcode stats)] [-t (trim barcodes off after demultiplexing)]\n" unless $barcode_file && $solexa_file;

my ($barcodes, $barcode_to_name) = read_barcodes($barcode_file);
split_solexa($barcodes, $solexa_file, $barcode_out, $barcode_to_name);

sub split_solexa {
	my ($barcodes, $in, $barcode_out, $barcode_to_name)=@_;
	open (IN, $in) or die "can't open $in: $!\n";

	my $num_read = 0;
	my $num_missed = 0;
	my %counts;
	my $format_type = 0;

	my $barcode_len = get_barcode_len($barcodes);

	if ($barcode_out) {
		open (OUT, ">$barcode_out") or die "can't open $barcode_out: $!\n";
	}

	while (<IN>) {
		chomp;
		$format_type = infer_format_type($_) unless $format_type;	

		my $missed_barcode = 1;
		my @pieces = split /:/, $_;
		my $bc = substr($pieces[5], 0, $barcode_len);

		if (my $fh = $barcodes->{$bc}) { # match barcode
			$missed_barcode = 0;
			$counts{$bc}++;
			# Do trimming if -t flag was used; otherwise, do not modify sequence or quality scores
			if ($trim_flag) {
				$pieces[5] = substr($pieces[5], $barcode_len);
				if ($format_type == $SCARF_ASCII) {
					$pieces[6]=~ s/^\s*\D{$barcode_len}//;
				}
				elsif ($format_type == $SCARF) {
					$pieces[6] =~ s/^\s*([-\d]+\s){$barcode_len}//
				}
			}
			my $seq = join ":", @pieces;
			print $fh "$seq\n";
		}

		$num_missed++ if $missed_barcode;
		$num_read++;
	}
	my $num_matched = $num_read-$num_missed;
	printf STDERR "\n\t\tSUMMARY STATS\n";
	printf STDERR "%d of %d (%.1f%%) of reads had a barcode match\n", $num_matched, $num_read, ($num_matched/$num_read)*100;
	printf STDERR "#%d\t%d\t%f\n", $num_matched, $num_read, ($num_matched/$num_read);
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

sub infer_format_type {
	shift;

	if(/^\S+:\d+:\d+:\d+:\d+.?.?.?.?:\D+$/) {
		return $SCARF_ASCII;
	}
        elsif(/^\S+:\d+:\d+:\d+:\d+/) {
		return $SCARF;
	}

}

sub read_barcodes {
	my $in=shift;
	open (IN, $in) or die "can't open $in: $!\n";

	my @fileHandles;
	my %barcodes;
	my %seq_to_barcode_name;
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
		$barcodes{$barcode} = $out; # barcode points to the file handle
		if ($barcode_name) {
			$seq_to_barcode_name{$barcode} = $barcode_name;
		}
	}

	close IN;
	
#	return \@fileHandles, \@barcodes;
	return \%barcodes, \%seq_to_barcode_name;
}
