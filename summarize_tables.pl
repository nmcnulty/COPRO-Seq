#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $input;				# May be a comma-separated list of files or a directory path
my $col_for_row_names;
my $col_for_values;

GetOptions 	(	
	'i|input=s'			=> \$input,
	'r|row=i'			=> \$col_for_row_names,
	'v|values=i'		=> \$col_for_values
) or die "Failed to initiate script.  Please check your parameters.";

my @filelist;

if (-d $input) {		# If a directory path was passed...
	opendir (DIR, $input) || die "Error: Could not open directory $input.\n";
	my @templist = readdir(DIR);
	foreach my $f (@templist) {
		unless (($f eq ".") || ($f eq "..")) {
			push(@filelist, $input . "\/" . $f);
		}
	}
	closedir DIR;
}

else {
	@filelist = split(/,/, $input);
	foreach(@filelist) {	s/\s//;	}	# Eliminate any whitespace between elements of comma-separated list
}

summarize_tables(\@filelist, $col_for_row_names, $col_for_values);

exit;

sub summarize_tables {
	my ($files_ref, $col_for_rows, $col_for_data, $headers_ref) = @_;		# Use of optional headers is not yet implemented...
	$col_for_rows--;
	$col_for_data--;
	my @sortedfilenames = sort {$a cmp $b} @$files_ref;
	my %HoH;			# Structure: filename -> table row name -> table value
	my %all_rownames_seen;	# Keep track of every row name (e.g. species name, barcode, uniqueness code, etc.) across all files
	my %col_headers;
	my $tag_present;
	foreach my $f (@sortedfilenames) {			# Files will be opened in sorted order
		($tag_present, $col_headers{$f}) = get_source_tag($f);
		#print $col_headers{$f} . "\n";
		open (TEMP, $f) || die "Error: Can't open $_.\n";		# Find file in file system
		if($tag_present == 1) { <TEMP>; }						# If there's a tag, need to ignore it for reading in info
		while (my $templine = <TEMP>) {
			$templine =~ s/[\n\r]//;	# Eliminate newline character at end of line
			my @tempval = split (/\t/, $templine);
			# Check to see if current genome is included in the running tally of all genomes encountered
			if (!defined $all_rownames_seen{$tempval[$col_for_rows]}) { $all_rownames_seen{$tempval[$col_for_rows]} = 1; }
			$HoH{$col_headers{$f}}{$tempval[$col_for_rows]}=$tempval[$col_for_data];
		}
		close TEMP;
	}

	my @colheaders;
	foreach (keys %col_headers) {	push(@colheaders, $col_headers{$_});	}
	my @sortedcolheaders = sort {lc($a) cmp lc($b)} (@colheaders);

	my $num_of_columns = scalar @sortedcolheaders;
	my @sortedrownames = sort {lc($a) cmp lc($b)} (keys %all_rownames_seen);
	
	# Output HoH contents to new summary table of all results with appropriate column and row labels 
	foreach(@sortedcolheaders) {	print "\t$_";	}	# Print column headers
	print "\n";
	foreach (@sortedrownames) {
		print "$_\t";									# Print row label
		# Print all values separated by tabs, except the one from the last column (which needs a newline after it) 	
		for (my $i = 0; $i < ($num_of_columns - 1); $i++) {
			if (defined $HoH{$sortedcolheaders[$i]}{$_}) { print "$HoH{$sortedcolheaders[$i]}{$_}\t"; }
			else { print "0\t"; }
		}
		# Print last value for the row
		if (defined $HoH{$sortedcolheaders[$num_of_columns-1]}{$_}) { print "$HoH{$sortedcolheaders[$num_of_columns-1]}{$_}\n"; }
		else { print "0\n"; }
	}
}

# Input: a table file path
# Output: if a source tag is detected, return source name (e.g. for <source = 1-1>, would return string "1-1"); otherwise, return file name
sub get_source_tag {
	my $path = shift;
	open(INPUT, "$path") || die "Error: Can't open file $path.\n";
	my $firstline = <INPUT>;
	close INPUT;
	if ($firstline =~ /<(.+?)>/) {
		return (1, $1);
	}
	else {
		return (0, $path);
	}
}