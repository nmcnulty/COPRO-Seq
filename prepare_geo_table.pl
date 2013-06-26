#! /usr/bin/perl
# prepare_geo_table.pl
# Author: Nate McNulty
# This script is intended to be run after the COPRO-Seq pipeline (invoked using the -G switch) 
#	has completed its operations (using the -G switch)
#	Important: Cannot run this after cleanup.sh has been run (mention this in manual)

# Improvements needed: currently, this script is not equipped to consider genomes specified in
#	the 'external genome paths' field of project.info when specifying 'organism' in the GEO table
#	file created; this may take a fair bit of effort to implement, and for now, I would recommend
#	just using workarounds (i.e., manually adjusting the table at the end)

# Usage example:
#	perl prepare_geo_table.pl -g A -m ./inputs/NM1200.mapping -s ./inputs/project.info

# Inputs:
#	Location of project.info file (spreadsheet) which has been previously downloaded
#	Group designation (should be same as group used in COPRO-Seq analysis)
#	Location of .mapping file (same as used in previous COPRO-Seq analysis)
#		Copy .mapping file to folder "/GEO/inputs" when running batch_coproseq.pl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use FindBin;
use lib "$FindBin::Bin/lib";	# Tell perl where to find barcodes.pm, genomecodes.pm (looks in folder where prepare_geo_table.pl is stored)
use barcodes;
use genomecodes;
use Digest::MD5 qw(md5_hex);

my ($group, $mapping_file_path, $output_file_path, $project_file_path);

GetOptions(
	"g|group=s"					=> \$group,
	"m|mapping=s"				=> \$mapping_file_path,
	"o|output=s"				=> \$output_file_path,
	"s|spreadsheet=s"			=> \$project_file_path
) or die usage();

my $groups = get_groups($group);
my $project_data_hash = table_to_hash($project_file_path);
my $filtered_data_hash = filter_by_group($project_data_hash, $groups);
my $mapping_hash = load_mapping_file($mapping_file_path);

# Overall process:
# =================
# Open output table file
# Print the header line

# Cycle through each entry in the filtered_data_hash and:
#	i) Look for sample-specific hits file
#	ii) Calculate md5 checksum for hits file
#	iii) Look for sample-specific SCARF file
#	iv) Calculate md5 checksum for SCARF file

#	v) Print "run"[run#]_lane[lane#]_[barcode sequence]_[sample name]
#	vi) Print [sample name]
#	vii) Print blank column
#	viii) Print genome names
#	ix) Print 3 blank columns
#	x) Print "genomic DNA"
#	xi) Print contents of "GEO_description" field in Google doc
#	x) Print name of .hitratios file
#	xi) Print blank column
#	xii) Print "txt"
#	xiii) Print MD5 checksum for .hitratios file
#	xiv) Print name of .scarf file
#	xv) Print "Illumna_native (SCARF)
#	xvi) Print the md5 checksum for .scarf file

# Close output table file

open(OUTPUT, ">$output_file_path") || die "Can't open output file $output_file_path!\n";
print OUTPUT "Sample name\ttitle\tsource name\torganism\tcharacteristics: tag\tcharacteristics: ",
	"tag\tcharacteristics: tag\tmolecule\tdescription\tprocessed data file name\t",
	"processed data file type\tprocessed data file MD5 checksum\traw file name\traw file type\traw file",
	" MD5 checksum\tinstrument model\tread length\n";
for my $p (@$filtered_data_hash) {		# For each row of the spreadsheet (i.e. each run/lane combo)...
	my $filteredsamples = get_filtered_sample_list($p, $mapping_hash);
	for my $s (@$filteredsamples) {		# For each sample in this row's run/lane combo...
		my $prefix = "machine".$p->{machine}."_run".$p->{run}."_lane".$p->{lane}."_".$$mapping_hash{$p->{pool}}{$s}."_".$s;
		unless (-f "$prefix\.hitratios\.gz") {	die "Cannot find one of your compressed .hitratios files: $prefix\.hitratios\.gz\n"; }
		my $hitsmd5 = md5_hex("$prefix\.hitratios\.gz");
		# Eventually, change this to look for a .scarf.gz and .fastq.gz file (after specifying that split_scarfs.sh pass everything through gzip first
		my ($seqmd5, $seq_format);
		if (-f "$prefix\.fastq\.gz") {
			$seqmd5 = md5_hex("$prefix\.fastq\.gz");
			$seq_format = 'fastq';
		}
		elsif (-f "$prefix\.scarf\.gz") {
			$seqmd5 = md5_hex("$prefix\.scarf\.gz");
			$seq_format = 'scarf';
		}
		else {
			die "Cannot find a compressed, demultiplexed sequence file for prefix $prefix!\n";
		}
 		
 		print OUTPUT	"machine".$p->{machine}."_run".$p->{run}."_lane".$p->{lane}."_".$$mapping_hash{$p->{pool}}{$s}."_".$s."\t",	# "Sample name"
 						$s."\t",																									# "title"
 						'Defined bacterial assemblage from the feces of a gnotobiotic mouse'."\t",									# "source"
 						process_genome_field($p->{internal_genome_accessions})."\t",												# "organism"
 						"\t",																										# "characteristics: tag"
 						"\t",																										# "characteristics: tag"
 						"\t",																										# "characteristics: tag"
 						'genomic DNA'."\t",																							# "molecule"
 						"\t",																										# "description (complete within Excel later)"
 						"$prefix\.hitratios\.gz"."\t",																				# (processed data file) "file name"
 						'txt'."\t",																									# (processed data file) "file type"
 						$hitsmd5."\t",																								# (processed data file) "file checksum"
 						$prefix."\.".$seq_format."\.gz\t",																				# (raw data file) "file name"
 						$seq_format."\t",																							# (raw data file) "file type"
 						$seqmd5."\t",																								# (raw data file) "file checksum"
						$p->{platform}."\t",																						# "instrument model"
						"\t\n";																										# "read length"
	}
}

close OUTPUT;
exit;

sub process_genome_field {
	my $genomes_rough = shift;
	my $genomes_reformatted;
	my %abbrev_to_inhouseacc = genomecodes::declare_inhouseacc_for_genomeabbrev(); # Key = abbreviation (e.g. BACCAC), value = in-house accession #
	my %inhouseacc_to_name = genomecodes::declare_speciesname_for_inhouseacc();
	my @tmp = split(/,/, $genomes_rough);
	for my $k (@tmp) {
		my $name;
		$k =~ s/^\s//;	# Eliminate leading whitespace if included after commas e.g. BACAC, BACTHE, ...
		if (defined $abbrev_to_inhouseacc{$k}) {	# i.e. If an abbreviation was provided in the spreadsheet
			$name = $inhouseacc_to_name{$abbrev_to_inhouseacc{$k}};
		}
		elsif (defined $inhouseacc_to_name{$k}) {
			$name = $inhouseacc_to_name{$k};
		}
		else {	die "One of the genomes in your spreadsheet, $k, does not appear to be an accepted",
					" abbreviation or microbialomics accession number!\n"; 
		} 
		$name =~ s/_/ /g;
		if (defined $genomes_reformatted) {		# If this is not the first entry
			$genomes_reformatted .= ", $name";
		}
		else {
			$genomes_reformatted = $name;
		}
	}
	return $genomes_reformatted;
}

sub load_mapping_file {
 	my $file_path = shift;
 	my %barcode_abbrev_to_seq = barcodes::declare_barcodes();
 	my %mapping_table;		# Will actually be a HoH with structure $HoH{pool name}{sample name}=literal barcode sequence
 	my %sample_names_seen;	# To help verify there is no sample name duplication within a mapping file (which could wreak havoc on later steps of the pipeline)
 	my $bclength;
 	open (MAPPING, $file_path) || die "Can't open mapping file at $mapping_file_path!\n";
 	while (<MAPPING>) {
 		chomp;
 		my @line = split(/\t/,$_);
# 		# Make everything in barcode designation uppercase (e.g. 4bc4 = 4BC4, atag = ATAG) to ensure matching is okay later
 		$line[2] =~ tr/[a-z]/[A-Z]/;
 		if ($barcode_abbrev_to_seq{$line[2]}) {				# i.e. if a barcode abbreviation was used and can be found in barcodes.pm
 			$line[2] = $barcode_abbrev_to_seq{$line[2]};	# Swap the barcode abbreviation for the appropriate barcode sequence
 		}
 		if (! exists $sample_names_seen{$line[1]}) {					# Ensures no duplication of sample names (which guarantees a particular pool/sample combo has not been defined already) 
 			$mapping_table{$line[0]}{$line[1]} = $line[2];
 			$sample_names_seen{$line[1]} = 1;
 		}
 		else {
 			die "Error: There is more than one instance of sample name $line[1] in your mapping file. ",
 				"All sample names within your mapping file must be unique to avoid problems in downstream data summarization.\n";
 		}
 	}
 	close MAPPING;
 	return \%mapping_table;
}


# Convert table in project.info to a reference to an array of references to hashes
#	Each array reference corresponds to a row in the spreadsheet
#	Each hash key = a column header
#	Each hash value = spreadsheet value
sub table_to_hash {
	my $in = shift;
	open (IN, $in) or die "Can't open project information file $in: $!\n";
	my $line = <IN>;	# First line is tab-separated headers
	chomp $line;
	my @headers = split /\t/, $line;
	for (@headers) {	
		$_ = lc($_);
		$_ =~ s/\s/_/g;	# Convert any whitespace in headers to underscores (allows spreadsheet to be a little prettier)
	}
	# Headers will be "group", "pool", "run", "lane", etc.
	my @table;		# Array of references to hashes containing sample info
	while (<IN>) {
		chomp;
		my @cols = split /\t/;
		my %hash;
		for my $i (0 .. $#cols) {				# $# is the subscript of the last element of the array
			$hash{$headers[$i]} = $cols[$i];	# Key: Google spreadsheet column header; Value: spreadsheet value
												# e.g. $hash{"lane"} = 8;
		}
		push @table, \%hash;	# Add reference to new hash (sample) to array of references, @table
	}
	close IN;
	return \@table;
}

sub get_filtered_sample_list {
	
	my ($datahash, $mapping_hash) = @_;
	
	my (%allsamples, %samples_to_remove);
	my @filteredsamples;
	my @excludedsamples = split(/,/, $datahash->{samples_to_exclude});
	
	foreach(keys %{$$mapping_hash{$datahash->{pool}}}) { $allsamples{$_} = 1; }
	
	foreach(@excludedsamples) {	
		s/\s//; 	# Eliminate any spaces that were between elements in comma-separated list
		$samples_to_remove{$_} = 1;
	}
	
	foreach(keys %allsamples) {
		if (! exists $samples_to_remove{$_}) {
			push(@filteredsamples, $_);
		}
	}
	
	return \@filteredsamples;
}

# Take $project_data_hash and return an array of references that only covers rows from the spreadsheet
#	with the proper group code
sub filter_by_group {
 	my ($ref_to_spreadsheet_hash, $ref_to_group_hash) = @_;
 	my @filteredtable;
 	foreach my $p (@$ref_to_spreadsheet_hash) {	# i.e. foreach sample in the spreadsheet...
 		# If "group" is blank for a given sample on the spreadhseet or the value in the "group" cell
 		# is not one of the acceptable group codes, skip it
 		if (!$p->{group} || !$$ref_to_group_hash{$p->{group}}) {
			next;
		}
 		push @filteredtable, $p;
 	}
 	return \@filteredtable;
}

sub get_groups {
	my %groups;
	my @a = split(/,/, $_[0]);
	foreach (@a) {
		$groups{$_} = 1;
	}
	return \%groups;
}

sub usage {
	print "Usage: perl prepare_geo_table.pl -g [group] -m [.mapping file] -o [output file] -s [project.info file]\n";
}
