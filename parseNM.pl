#! /usr/bin/perl

use strict;
use warnings;

# $ARGV[0] = .elandout input file
# $ARGV[1] = desired filename for NM sequences
# $ARGV[2] = desired filename for non-NM sequences (optional)

my $prefix;
if ($ARGV[0] =~ /(.+)\.elandout$/) {
	$prefix = $1;
}
else {
	die "Your input file must be of type .elandout\n";
}

if ($ARGV[1]) {
	system("grep '^.*[[:space:]].*[[:space:]]NM[[:space:]]' $ARGV[0] > $prefix\_NM.temp");
	make_fasta_from_grep("$prefix\_NM.temp", $ARGV[1]);
	system("rm $prefix\_NM.temp");
}
else {
	die "You must provide an output filename for NM sequences when callign parseNM.pl\n";
}

if ($ARGV[2]) {
	system("grep -v '^.*[[:space:]].*[[:space:]]NM[[:space:]]' $ARGV[0] > $prefix\_non-NM.temp");
	make_fasta_from_grep("$prefix\_non-NM.temp", $ARGV[2]);
	system("rm $prefix\_non-NM.temp");
}

# Inputs: .temp output file created by grep, output filename
sub make_fasta_from_grep {
	my ($tempfile, $output) = @_;
	open(TEMP, $tempfile) || die "Can't open $tempfile\n";
	open(OUTPUT, ">$output") || die "Can't open $output\n"; 	
	while(<TEMP>) {
		my @line = split(/\t/, $_);
		print OUTPUT "$line[0]\n$line[1]\n";
	}
	close TEMP;
	close OUTPUT;
}