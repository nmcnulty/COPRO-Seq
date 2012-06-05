#!/usr/bin/perl
# Written by: Nate McNulty
# Last updated: Oct-2010
# First attempt at a script to go from eland output to a set of tracks file (one per species)
#	that can be visualized easily in R or LWGV to identify any problems in the data.  
# Examples of such problems might include:
#	Genome 'jackpots' (e.g. high-copy plasmid sequences which are highly over-represented in the 
#		sequencing results)
#	Segments of the genome that do not get any reads and therefore should perhaps not be included 
#		in the draft genome

# Future improvements:
# Add a check when people invoke -z to be sure all genomes that need to be zeroed out (i.e. that 
#	are in the eland output) are in the genomes folder
# Create a hash with species => # describing scale needed for LWGV when -l is invoked
#	Need to keep track as HoA is being populated for each species how big the biggest # of hits is
# Change the way -z option is executed; instead of requiring a folder of genomes, it would be
#	much more straightforward to require a table of genome lengths (one file is much easier to 
#	keep track of and keep updated than a folder of FASTA sequences)

# IMPORTANT: Looks like species names that contain hyphens may not work (due to LWGV formatting 
#	requirements)!

use strict;
use warnings;
use Getopt::Long;

my ($genomesdir, $samplename, $lwgv, $suffix, $zeros);
GetOptions (
	'genomesdir=s'	=> \$genomesdir,
	'lwgv'			=> \$lwgv,
	'name=s'		=> \$samplename,
	'suffix=s'		=> \$suffix,
	'table=s'		=> \$table_of_genome_lengths,
	'zeros'			=> \$zeros,
	'modzeros'		=> \$modzeros
);
my $check = options_ok();
unless ($check eq 'ok') { usage($check); }

####################################################################################################
my %tracks; # HoA; Key = species name, Array indices represent positions in genome, value = # times position hit during alignment
if ($zeros) { zero_out_HoA(\%tracks, get_fasta_files()); }
elsif ($modzeros) { zero_out_HoA_mod(\%tracks, get_genome_lengths_from_table_file()); }
populate_hash_from_eland_output(\%tracks);
print_HoA_to_file(\%tracks, $lwgv);
exit;
####################################################################################################

# Every line in output file is: position<tab>hit count
sub print_HoA_to_file
{
	my ($hashtoprint_ref, $lwgv_format) = @_;
	my %HoA = %$hashtoprint_ref;
	for my $a (sort keys %HoA) { # For each genome...
		my $highest = 0;	# Use to keep track of highest value of any coordinate
		# Create a separate file for each species
		if ($suffix) { open(TABLE, ">$a\.$suffix") || die "Can't open output file $a\.$suffix!\n"; }
		else { open(TABLE, ">$a\.tracks") || die "Can't open output file $a\.tracks!\n"; }
 		if ($lwgv) {
			my $total = 0;
			# First cycle through all values to tally total hits (for _e# value in beginning of file)
			for $b (0 .. $#{$HoA{$a}}) {		# For each position in the genome 0->n...
				if (defined ${HoA{$a}}[$b]) { 
					$total += $HoA{$a}[$b]; 
					if (${HoA{$a}}[$b] > $highest) {
						$highest = ${HoA{$a}}[$b];
					}
				}
			}
			if ($total > 0) {
				$a =~ s/-/_/;	# Substitute - with _ to keep the LWGV software happy
				if (defined $samplename) {
					print TABLE "graph $a\_$samplename\_e".round(log($highest))." addPoints(";
				}
				else {
					print TABLE "graph $a"."_e".round(log($highest))." addPoints(";
				}
				my $coords = "";
				for $b (0 .. $#{$HoA{$a}}) {		# For each position in the genome 0->n...				
					if (defined ${HoA{$a}}[$b]) {
						my $value = sprintf '%.2f', log($HoA{$a}[$b]);	#I think J's scripts may require numbers with two decimal places 
						$coords .= "$b\:$value,";
					}
				}
				$coords =~ s/,$//;	# Remove trailing "," from last element
				print TABLE $coords;
				print TABLE ")";
			}
 		}
 		else {
			for $b (0 .. $#{$HoA{$a}}) {		# For each position in the genome 0->n...
				if (defined ${HoA{$a}}[$b]) {
					print TABLE "$b\t$HoA{$a}[$b]\n";
				}
			}
		}
 		close TABLE;
 		print "Finished writing results for species $a\n";
 	}
}

# Defines (with zeros) every position in every genome in the genomes directory
# This is done to ensure every base ends up in the output, regardless of whether it is ever hit
#	in the Eland alignments; this allows distribution plots generated with output files to include
#	0 to give a sense for how often a base in a genome was never hit
sub zero_out_HoA {
	my ($HoA_ref, $genome_files_ref) = @_;
	foreach my $a (@$genome_files_ref) {
		open (INPUT, "$genomesdir\/$a") || die "Can't open genome file $a!\n";
		my $strippedfilename = $a;
		$strippedfilename =~ s/\.\w+//;	# Remove filename extension (e.g. .fna/.fas/.fa)
		my $sequence = '';
		while(<INPUT>) {
			unless (/^>/) {			# Ignore headers - we only care about sequence data
				$sequence .= $_;
			}
		}
		close INPUT;
		print "The length of the genome stored in $a is ".length($sequence)."bp\n";
		# Cycle through HoA key for species and set values for all positions (array elements) to zero
		for (my $i=0; $i < length($sequence); $i++) {
			$HoA_ref->{$strippedfilename}[$i] = 0;	# Not sure if this will work
		}
	}
}

# Usage: zero_out_HoA_mod(<hash containing tracks results>, <reference to hash containing genome lengths>)
sub zero_out_HoA_mod {
	my ($HoA_ref, $lengths_table_ref) = @_;
	foreach my $a (keys %{$lengths_table_ref}) {		# For each genome in the table of genome lengths
		for (my $i=0; $i < $lengths_table_ref->{$a}; $i++) {
			$HoA_ref->{$a}[$i] = 0;
		}
	}
}

sub populate_hash_from_eland_output {
	my $HoA = shift;
	open(ELANDOUT, $ARGV[0]) || die "Can't open file $ARGV[0]!";
	my $readcount = 0;
	while (<ELANDOUT>) {
		$readcount++;
		my @line = split("\t", $_);
		if (defined $line[6]) {		# Workaround for now, given that hyphens in hash key names seem to be causing problems
			$line[6] =~ s/-/_/;
		}
		if ($line[2] eq 'U0') {	# Must be a unique match for now; ultimately use whatever U# specified in the coproseq pipeline (easier if this code is integrated into the rest of the COPRO-Seq code)
			if ($line[8] eq 'F') {
				for (my $i=0; $i < length($line[1]); $i++) {
					if (defined $HoA->{$line[6]}[$line[7]+$i]) {
						$HoA->{$line[6]}[$line[7]+$i]++;	# Add one to that position to count the current read's 'hit'
					}
					else {  # If first instance at this position, value will be undefined, and so can't just use ++ shortcut 
						$HoA->{$line[6]}[$line[7]+$i] = 1;
					}
				}
			}
			elsif ($line[8] eq 'R') {
				for (my $j=0; $j < length($line[1]); $j++) {
					if (defined $HoA->{$line[6]}[$line[7]-$j]) {
						$HoA->{$line[6]}[$line[7]-$j]++;
					}
					else {	# If first instance at this position, value will be undefined, and so can't just use ++ shortcut 
						$HoA->{$line[6]}[$line[7]-$j] = 1;
					}
				}
			}
			else {
				die "There is something wrong with the formatting of your eland results.  Cannot identify",
					" the strand for the read at line $readcount of your file $ARGV[0].\n";		
			}
		}
		if ($readcount % 100000 == 0) {
			print "$readcount sequences done.\n";
		}
	}
	print "$readcount sequences done.\n";
}

# Important that genome descriptors in table file match those in ELAND output EXACTLY
sub get_genome_lengths_from_table_file {
	open(TABLE, $table_of_genome_lengths) || die "Can't open genome lengths table file $table_of_genome_lengths!\n";
	my %lengths;
	while(<TABLE>) {
		chomp;
		my @line = split(/\t/, $_);
		$lengths{$line[0]} = $line[1];
	}
	close TABLE;
	return \%lengths;
}

sub get_fasta_files {
	opendir (GENOMES, $genomesdir) || die "Can't open directory $genomesdir!\n";
	my @all = readdir(GENOMES);
	my @keep;
	foreach (@all) {
		if ((/.fna/)|(/.fas/)|(/.fa/)) {	# If a (presumably genome) fasta file; also avoids problems with "." and ".."
			push(@keep, $_);
		}
	}
	closedir GENOMES;
	return \@keep;
}

sub options_ok {
	if ($zeros) {
		unless (defined $genomesdir)	{ return 'Invoking -z requires the invokation of the -g option as well'; } 
	}
	if ($modzeros) {
		unless (defined $table_of_genome_lengths)	{ return 'Invoking -m requires specifying a table of genome lengths using -t as well'; }
	}
	return 'ok';
}

sub usage {
	my $error = shift;
	print	"\n$error\n";
	print	"\nUsage: perl eland2tracks.pl <eland output> (options)\n\n",
			"Options:\n",
			"\t-n : Name to differentiate sample from other samples within microbialomics\n",
			"\t-l : Report results in LightWeight Genome Viewer (LWGV)-compatible format\n",
			"\t-g : Path of reference genomes\n",
			"\t-z : Report positions even if they have zero hits\n\n"; 
	exit;
}

sub round {
    my($number) = shift;
    return int($number + .5);
}