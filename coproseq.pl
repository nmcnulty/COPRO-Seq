#!/usr/bin/perl
# Author: Nathan McNulty
# Last Updated: 06/29/2012

# THIS SOFTWARE IS CURRENTLY INCOMPLETE AND NOT FIT FOR DISTRIBUTION OUTSIDE OF THE LAB.

# IMPROVEMENTS FOR THE FUTURE
#	Streamline trimming of initial SCARF and hash loading (i.e. don't create intermediate trimmed file - just load straight to memory to cut down
#		on file server pain)
#	Do some validation of options passed to the script to make sure they look good (e.g. finding SquashedGenomes directory)
#	Create a subroutine for collecting stats for each SCARF input (i.e. quality scores, base composition stats,
#		how many reads, etc.); what should be returned from this subroutine?
#	Work on memory leaks; Perl appears to be using far more memory than should be needed (likely some memory is not being released back to memory pool)

use strict;
use warnings;
# FindBin variable '$FindBin::Bin' stores the path to coproseq.pl
use FindBin;
# Provide location of barcodes.pm
use lib "$FindBin::Bin/lib";
use barcodes;
use Getopt::Long;
use Cwd;
use IO::File;

my $allfiles = '';							# Default: all intermediate files not related to a used barcode are deleted/not created
my $bcseqs;									# Default: barcode sequence file name is undefined
my $input;									# Default: location of input file is undefined
my $genomesdir;								# Default: location of squashed genomes directory is undefined
my $outputdir;								# Default: output will be written to directory named prefix in folder from which script was called (declared later)
my $prefix = (time)."-".int(rand(1000));	# Default: the prefix preceding all created files will be (system time)-(random # from 1-1000)
my $readsize = 29;							# Default: program will use the first 29bp of each Illumina read unless told to use more/less
my $mismatches_allowed = 0;					# Default: program will not use any reads that eland requires 1 or 2 mismatches to map in creating summaries
my $IGS_table;								
my $single_cpu;								# Default: assume that alignment processes will be distributed across multiple CPUs in a Sun Grid Engine environment
my $tags;									# Default: don't include a source tag in the table output files (.mappingstats, .hitratios, etc.) produced

GetOptions (	'a|allfiles'		=> \$allfiles,
				'b|bcfile=s'		=> \$bcseqs,				# Must always be defined at startup (if processing a barcoded run)
				'g|genomesdir=s'	=> \$genomesdir,			# Must always be defined at startup
				'i|input=s'			=> \$input,					# Must always be defined at startup
				'l|length=i'		=> \$readsize,
				'm|mismatches=i'	=> \$mismatches_allowed,
				'n|normalization=s'	=> \$IGS_table,
				'o|outputdir=s'		=> \$outputdir,
				'p|prefix=s'		=> \$prefix,
				's|single_cpu'		=> \$single_cpu,
				't|tags'			=> \$tags
) or die usage();

if (! defined $outputdir)	{	$outputdir = $prefix;	}
#&check_options(\$bcseqs, \$genomesdir, \$input, \$readsize, \$mismatches_allowed, \$outputdir, \$prefix);

# DECLARATION OF CONSTANTS
my $WAIT_TIME = 10;				# Units = seconds; equal to number of seconds between each attempt to check if all eland alignments are finished
my $ALIGNMENT_MEMORY = '10G';	# Amount of RAM to allocate to Eland when doing alignments

# Locations of needed programs/scripts
my $parse_script_path =
	"$FindBin::Bin/parseNM.pl";
my $eland_executables_folder = 
	'/srv/cgs/local/gapipeline/latest/bin/';

# Generate directories that will be needed (if they haven't been created already; -d should possibly be -e below???)
if (! -d "$outputdir")					{	mkdir "$outputdir" || die "$!\n";					}
if (! -d "$outputdir/seqcountsbybc")	{	mkdir "$outputdir\/seqcountsbybc" || die "$!\n";	}
if (! -d "$outputdir/bcsortedseqs")		{	mkdir "$outputdir\/bcsortedseqs" || die "$!\n";		}
if (! -d "$outputdir/elandresults")		{	mkdir "$outputdir\/elandresults" || die "$!\n";		}
if (! -d "$outputdir/hitratios")		{	mkdir "$outputdir\/hitratios" || die "$!\n";		}
if (! -d "$outputdir/mappingstats")		{	mkdir "$outputdir\/mappingstats" || die "$!\n";		}
if (! -d "$outputdir/NM" )				{	mkdir "$outputdir\/NM" || die "$!\n";				}
if (! -d "$outputdir/filteredseqs")		{	mkdir "$outputdir\/filteredseqs" || die "$!\n";		}
if (! -d "$outputdir/monitoring")		{	mkdir "$outputdir\/monitoring" || die "$!\n";		}
if (! -d "$outputdir/filteredseqs/adapter")	{	mkdir "$outputdir\/filteredseqs/adapter" || die "$!\n";	}
if (! -d "$outputdir/filteredseqs/mouse")	{	mkdir "$outputdir\/filteredseqs/mouse" || die "$!\n";	}
if (! -d "$outputdir/filteredseqs/unknown")	{	mkdir "$outputdir\/filteredseqs/unknown" || die "$!\n";	}

# Open BC file, check formatting, and store file info in a hash
my $ref_to_samples_by_bc_hash = parse_bc_file($bcseqs);
my %samplesbybarcode = %$ref_to_samples_by_bc_hash;
# For %samplesbybarcode, key = barcode sequence; value = sample name
my @barcodes = keys %samplesbybarcode;

# Confirm barcodes are of consistent size and report their length
my $bcsize = check_bc_size(\@barcodes);
summarize_bcs(\@barcodes);

# Detect sequence file format, parse sequence file, and report bc tally stats
my $seq_file_type = detect_seq_file_type($input);
my ($ref_to_bc_tallies_hash, $readnumber) = parse_seq_file_by_bc($seq_file_type, \@barcodes, $bcsize);
report_bc_tallies($ref_to_bc_tallies_hash, $readnumber, \@barcodes);

# Create file of alignment jobs and either submit it to Sun Grid Engine or run it locally
my $ref_to_array_of_jobs_files = make_alignment_jobs_file(\@barcodes, $readsize, $bcsize);
if ($single_cpu) { 
	foreach (@$ref_to_array_of_jobs_files) {
		system("sh $_");
	}
}
else { 
	foreach (@$ref_to_array_of_jobs_files) {
		system("sbatch --mem=$ALIGNMENT_MEMORY $_");
	}
}

########################################################################################################################################################
# PURPOSE:	Create a jobs file for each Eland submission and pass it to nq
# NOTES:	I would have preferred to name the jobs file "$prefix_elandjobs", but in cases where $prefix begins with a number, qsub will return an error
#			Could get around this by just using some particular letter at the beginning, like "X_$prefix_elandjobs"; this would allow me to selectively
#				remove the appropriate elandjobs files at the end of the script without interrupting other submissions from other input files
#			If for some reason you suspect Eland will take longer than 1 hour to run, you must add in the "-P long" option to ensure jobs are passed
#				to the long queue
sub make_alignment_jobs_file {
	my ($barcodes_array_ref, $readsize, $bcsize) = @_;
	my $elandtype = "eland_".($readsize - $bcsize);
	my @array_of_jobs_files;
	foreach (@$barcodes_array_ref) {
		my $jobs_file_path = "$outputdir\/elandjobs\_$prefix\_$_\.job";
		push(@array_of_jobs_files, $jobs_file_path);
		open (TASKSFORBC, ">$jobs_file_path") || die "Error: Can't create $jobs_file_path\n\n";
		if (-s "$outputdir\/bcsortedseqs\/$prefix\_$_\.fas") {	# If sequence file is not empty (-s returns true for file sizes > 0)
			print TASKSFORBC
			# Align to microbial references
				"$eland_executables_folder\/$elandtype $outputdir\/bcsortedseqs\/$prefix\_$_\.fas $genomesdir $outputdir\/elandresults\/$prefix\_$_\.elandout\n",
			# Pass sequences that did not align to microbial references to a new file (use as input for alignment to adapters)
				"$parse_script_path $outputdir\/elandresults\/$prefix\_$_\.elandout NM\/$prefix\_$_\_norefs.NM\n",
			# Align remaining sequences to adapter reference, allowing 2 mismatches
			# Squashed adapter genome should be distributed with scripts
				"$eland_executables_folder\/$elandtype NM\/$prefix\_$_\_norefs.NM $FindBin::Bin/filteringrefs/adapter NM\/$prefix\_$_\_adapter.elandout\n", 
			# Pass sequences that did not align to adapter reference to a new file (use as input for alignment to mouse)
			# Also pass sequences that DID align to adapter reference to a new file
				"$parse_script_path NM\/$prefix\_$_\_adapter.elandout NM\/$prefix\_$_\_noadapters.NM filteredseqs/adapter\/$prefix\_$_\_adapter.fna\n",
			# Align remaining sequences to mouse reference, allowing 2 mismatches
			# Squashed mouse genome should be distributed with scripts
				"$eland_executables_folder\/$elandtype NM\/$prefix\_$_\_noadapters.NM $FindBin::Bin/filteringrefs/mouse NM\/$prefix\_$_\_mouse.elandout\n",
			# Pass sequences that did not align to mouse reference to a new file (these will be the final, non-matching sequences of unknown origin)
			# Also pass sequences that DID align to mouse reference to a new file
				"$parse_script_path NM\/$prefix\_$_\_mouse.elandout filteredseqs/unknown\/$prefix\_$_\_unknown.fna filteredseqs/mouse\/$prefix\_$_\_mouse.fna\n",
			# Create dummy .done file to record that all tasks are completed
				"echo COMPLETE > NM\/$prefix\_$_\.done\n";
		}
		else {
			print TASKSFORBC "echo Skipping alignments for barcode $_ because there were no sequences in your input file matching this barcode.\n";
			# Create empty .elandout file so that progress monitoring doesn't break
			print TASKSFORBC "touch $outputdir\/elandresults\/$prefix\_$_\.elandout\n";
			print TASKSFORBC "echo COMPLETE > NM\/$prefix\_$_\.done\n";
		}
		close TASKSFORBC;
		# Append new jobs here with semicolons so that multiple tasks get passed to each node but are completed as a group
		# Create .done file when all tasks are finished; change code below to monitor for presence of all .done files instead of looking for .  files
	}
	return \@array_of_jobs_files;
}

unless ($single_cpu) {
	# MONITORING BELOW WILL ULTIMATELY HAVE TO BE CHANGED TO LOOK FOR .DONE FILES CREATED IN LINES ABOVE, RATHER THAN LOOKING FOR ELANDOUT FILE FOR REF ALIGNMENTS
	
	########################################################################################################################################################
	
	########################################################################################################################################################
	# PURPOSE:	For each barcode in @barcodes, check that an eland results file was generated and has a size of > 0
	# NOTES:	Eland creates an output file once it begins (with size 0) but writes all results at once (rather than writing them incrementally 
	#				like BLAST for several queries).  So, the size of any existing results file should be either 0 or the final size when it is checked
	#			It may take awhile for a file to be created (if the wait-time on the queue is long), but it shouldn't take too long for a file's size to become
	#				greater than 0 (because Eland usually only takes 3-5m to run); it might make sense to time-out in the latter case so that the script doesn't
	#				hang if something goes wrong with one of the jobs (in this case you could just re-submit the problematic jobs from the jobs file)
	#			Use global constant WAIT_TIME to dictate how long to wait between checking for files
	
	# Sleep for 5s after all files are created with size > 0 to be sure all are done being written to (probably not necessary, but just to be safe)
	
	# When all files return Message #3 (i.e. all Eland jobs are done), proceed with generating summaries
	
	# CONSIDER RE-WRITING FILE CHECKING SUB-ROUTINE SO THAT IT RETURNS FILE SIZE INSTEAD OF AN ERROR CODE THAT HAS TO BE DECODED
	my $starttime = scalar localtime();		# Get start time before beginning so user can monitor how long they've been waiting
	my $num_running = 0;
	my $num_finished = 0;
	my $num_notstarted = 0;
	print "\nStarted checking for results of alignment to reference genomes on $starttime\n\n";		# Add in line of code to figure out how long you've been waiting
	while ($num_finished < scalar(@barcodes)) {	# While there are still Eland output files that aren't made/done yet [for alignments to reference genomes]...
		$num_running = 0;
		$num_finished = 0;
		$num_notstarted = 0;
		foreach (@barcodes) {
			my $filestring = "$prefix"."_".$_.".elandout";
			my $statuscode = check_for_file("$outputdir\/elandresults\/$filestring");
			if ($statuscode == 0) { $num_notstarted++; }
			elsif ($statuscode == 1) { $num_running++; }
			elsif ($statuscode == 2) { $num_finished++; }
		}
		print "Waiting to start: $num_notstarted\tRunning: $num_running\tFinished: $num_finished\n";		
		sleep $WAIT_TIME;
	}
	print "\nEland has finished all alignments to the reference genomes specified.\n\n";
	
	# Recycle variables above
	$num_finished = 0;
	print "Checking status of no-match sequence filtering...\n";
	while ($num_finished < scalar(@barcodes)) {	# While there are still no-match records to filter against mouse/adapter genomes...
		$num_running = 0;
		$num_finished = 0;
		$num_notstarted = 0;
		foreach (@barcodes) {
			my $filestring = "NM\/$prefix\_$_\.done";
			my $statuscode = check_for_file($filestring);
			if ($statuscode == 0) { $num_notstarted++; }
			elsif ($statuscode == 1) { $num_running++; }
			elsif ($statuscode == 2) { $num_finished++; }
		}
		print "Waiting/Running: ".($num_notstarted + $num_running)."\tFinished: $num_finished\n";		
		sleep $WAIT_TIME;
	}
	print "\nAll filtering of no-match sequences against adapter and mouse genomes complete.\n";
}

# PARSE EACH ELAND OUTPUT INTO A SUMMARY FILE
# Use squashed genomes directory to create summary file fields (i.e. for each Eland results file, create a list of possible results where each possible result 
#	is a squashed genome (this is better than only listing the results found in each Eland results file because it will create a consistent output across all 
#	checked files, which should make importing into Excel much easier if some species do not show up in the alignments
# Future improvement: rewrite code such that a single matrix output is generated for all barcodes from the input
#	Where top row is sample name for BC1, sample name for BC2, etc.
#	And where first column from top to bottom is genome 1, genome 2, etc.

# Example Eland output line:
# >HWI-EAS158:1:139:302:287#0/1	AAGCAAAAACACCATT	U0	1	5	112	B_thetaiotaomicron.fas	3168430	F	DD
# Fields separated by tabs [clusterlocation, sequence, uniqueness code, # exact matches, # 1MM, #2MM, genome matched, genome position, genome strand, N code
# First 6 fields will always be present in current format

opendir SQGEN, $genomesdir || print "\n$!\n\n";
my @genome_files = readdir SQGEN;	# Each species included in the squashed genomes directory should have a file like this: XXXXX.2bpb
my @genome_names;
print "Compiling summary and mapping files for each of your " . (scalar @barcodes) . " barcoded samples...\n";
foreach (@genome_files) {
	if ($_ =~ /\.2bpb$/) {
		my $trimmed_name;
		($trimmed_name = $_) =~ s/\.2bpb$//;
		print "Found genome file for species with name $trimmed_name\n";
		#my @file_name_parts = split (/\./, $_);
		#push(@genome_names, $file_name_parts[0]);	# Grab the genome name component of the filename (e.g. "B_caccae" from "B_caccae.fas.2bpb")
		push(@genome_names, $trimmed_name);
	}
}
foreach (@barcodes) {
	my $q = $_;
	my $filestring1 = "$prefix"."_"."$samplesbybarcode{$_}"."_specieshits_".($readsize - $bcsize)."bp_".$mismatches_allowed."MM.hitratios"; # (e.g. "99_1_AATT_specieshits_16bp_0MM.txt")
	my $filestring2 = "$prefix"."_"."$samplesbybarcode{$_}".".mappingstats";
	open (SPECIES_HITS_STATS, ">$outputdir\/hitratios\/$filestring1") || die "Can't open $filestring1\n";
	open (MAPPING_STATS, ">$outputdir\/mappingstats\/$filestring2") || die "Can't open $filestring2\n";
	if ($tags) {
		print SPECIES_HITS_STATS "<$samplesbybarcode{$_}>\n";	# This source tag will allow summarize_tables.pl to properly label master summary columns with sample names
		print MAPPING_STATS "<$samplesbybarcode{$_}>\n";			# This source tag will allow summarize_tables.pl to properly label master summary columns with sample names
	}	
	open (ELANDRESULTS, "$outputdir\/elandresults\/$prefix"."_".$_.".elandout") || die "Can't open $prefix"."_".$_.".elandout";

	my %hit_counts;	# Key = genome name; Value = raw hits to genome;
	# Initialize %hit_counts keys and set all values to zero (have to do this because some keys may never appear in Eland results file)
	foreach (@genome_names) {
		$hit_counts{$_} = 0;
	}

	my %uniqueness_counts = (	
								"NM",			0,
								"NM-adapter",	0,
								"NM-mouse",		0,
								"NM-unknown", 	0,
								"QC",			0,
								"R0",			0,
								"R1",			0,
								"R2",			0,	
								"RM",			0,
								"U0",			0,
								"U1",			0,
								"U2",			0);
	my $wc_output_adapter	= `wc -l filteredseqs\/adapter\/$prefix\_$_\_adapter.fna`;
	my @wc_array1 = split(/\s/, $wc_output_adapter);
	$uniqueness_counts{'NM-adapter'}=($wc_array1[0] / 2);
	my $wc_output_mouse 	= `wc -l filteredseqs\/mouse\/$prefix\_$_\_mouse.fna`;
	my @wc_array2 = split(/\s/, $wc_output_mouse);
	$uniqueness_counts{'NM-mouse'}=($wc_array2[0] / 2);
	my $wc_output_unknown	= `wc -l filteredseqs\/unknown\/$prefix\_$_\_unknown.fna`;
	my @wc_array3 = split(/\s/, $wc_output_unknown);
	$uniqueness_counts{'NM-unknown'}=($wc_array3[0] / 2);

	foreach (<ELANDRESULTS>) {
		my @elandfields = split(/\t/, $_);				# Store each tab-delimited Eland output field
		chomp $elandfields[2];							# In some cases, a \n will be included and mess up hash assignments later
		$uniqueness_counts{$elandfields[2]}++;
		if ($elandfields[2] =~ /^U/) {				# If Eland-reported match code suggests a unique hit was found (either "U0", "U1", or "U2")
			# Must only count genome hits if the number of mismatches required to make a hit was less than or equal to the user-defined threshold
			if (($mismatches_allowed == 0) && ($elandfields[2] eq "U0")) {
				#my @genomefield = split(/\./, $elandfields[6]);		# Grab the genome Eland hit uniquely (which may have an extraneous extension, like ".fas")
				#$hit_counts{$genomefield[0]}++;
				$hit_counts{$elandfields[6]}++;						
			}
			elsif (($mismatches_allowed == 1) && ($elandfields[2] eq "U1" || $elandfields[2] eq "U0")) {
				#my @genomefield = split(/\./, $elandfields[6]);		# Grab the genome Eland hit uniquely (which may have an extraneous extension, like ".fas")
				#$hit_counts{$genomefield[0]}++;
				$hit_counts{$elandfields[6]}++;
			}
			elsif (($mismatches_allowed == 2) && ($elandfields[2] eq "U2" || $elandfields[2] eq "U1" || $elandfields[2] eq "U0")) {
				#my @genomefield = split(/\./, $elandfields[6]);		# Grab the genome Eland hit uniquely (which may have an extraneous extension, like ".fas")
				#$hit_counts{$genomefield[0]}++;
				$hit_counts{$elandfields[6]}++;						
			}
		}
	}
	close ELANDRESULTS;

	# Get totals for each hash for reporting percentages
	my $hitstotal = 0;
	foreach my $key (keys %hit_counts) {
		$hitstotal = $hitstotal + $hit_counts{$key};
	}
	# Write results to output files, sorting each hash by key to ensure a consistent output structure for every file
	if($IGS_table) {		# If an IGS table is available for going beyond raw hit calculations... (will this ever not be the case?)
		my $IGS_HoH_ref = table_file_to_HoH($IGS_table);
		my $total_norm_hits = 0;
		foreach my $g (sort keys %hit_counts) {
			if (defined $$IGS_HoH_ref{$readsize-$bcsize}{$g}) {
				if ($$IGS_HoH_ref{$readsize-$bcsize}{$g} > 0) {
					$total_norm_hits += $hit_counts{$g} / $$IGS_HoH_ref{$readsize-$bcsize}{$g};
				}
				else {
					die "\nThe informative genome size (IGS) for genome $g is not a positive number.  Please check your IGS table.\n";
				}
			}
			else {
				die "\nYour IGS table does not include any informative genome size for genome $g.  Please calculate a new IGS table that includes this genome and re-run your analysis.\n";
			}
		}
		if ($total_norm_hits > 0) {
			foreach my $genomekey (sort keys %hit_counts) { 
				my $norm_hits;
				$norm_hits = $hit_counts{$genomekey} / $$IGS_HoH_ref{$readsize-$bcsize}{$genomekey};
				# In cases where there are very few or no reads that map uniquely to any of the references, the line below will throw a division by zero error,
				#	so a check to make sure $total_norm_hits is not still 0 at this point is required
				if ($total_norm_hits > 0) {
					printf SPECIES_HITS_STATS "%s\t%d\t%.1f\t%.8f\n", $genomekey, $hit_counts{$genomekey}, $norm_hits, $norm_hits / $total_norm_hits;
				}
			}
		}
		else {
			foreach my $genomekey (sort keys %hit_counts) {
					printf SPECIES_HITS_STATS "%s\t%d\t0\t0\n", $genomekey, $hit_counts{$genomekey};		
					#printf SPECIES_HITS_STATS "%s\t%d\t%.1f\t%s\n", $genomekey, $hit_counts{$genomekey}, $norm_hits, "Div by 0";
			}
			print "\nWarning: There are zero results in the alignment output for barcode $q that meet your mapping criteria.  This sample generated no informative data.\n";
		}
	}
	else {
		foreach my $genomekey (sort keys %hit_counts) { 
			printf SPECIES_HITS_STATS "%s\t%d\n", $genomekey, $hit_counts{$genomekey};
		}
	}
	close SPECIES_HITS_STATS;

	my $matchcodestotal = 0;
	# Which codes does the following line exclude?  How is this diff. from the commented loop that follows?
	$matchcodestotal = $uniqueness_counts{NM}+$uniqueness_counts{QC}+$uniqueness_counts{R0}+$uniqueness_counts{R1}+$uniqueness_counts{R2}+$uniqueness_counts{RM}+$uniqueness_counts{U0}+$uniqueness_counts{U1}+$uniqueness_counts{U2};
# 	foreach my $key (keys %uniqueness_counts) {
# 		$matchcodestotal = $matchcodestotal + $uniqueness_counts{$key};
# 	}
	foreach my $matchkey (sort keys %uniqueness_counts) { printf MAPPING_STATS "%s\t%d\t%7.4f%%\n", $matchkey, $uniqueness_counts{$matchkey}, ($uniqueness_counts{$matchkey}/$matchcodestotal)*100; }
	close MAPPING_STATS;
}

exit;

sub check_for_file {
	my $fileName = shift @_;	# Grab filename from parameters

	if (-e "$fileName") {		# If system can find a file with the name you're looking for...
		my $size = -s $fileName;
		if ($size > 0) {
			return 2;	# 2 = Found your file and size is > 0
		}
		else {
			return 1;	# 1 = Found your file but size is still 0
		}
	}
	else {
		return 0;		# 0 = Failed to find your file
	}
}

sub check_bc_size {
	my $bc_array_ref = shift;
	my @bcs = @{$bc_array_ref};
	my $bcsize;
	foreach (@bcs) {
		if (!defined $bcsize) {
			# Start with first BC sequence
			$bcsize = length $_;
		}	
		else {
			# Check to make sure none of the sequence lengths are different from the first
			if ($bcsize != length $_) {
				die "Your barcode sequences do not appear to be a consistent length.\n\n";
			}
		}
	}
	return $bcsize;
}

# Input: file path
# Output: Returns 1 (OK) if no instances of the carriage return character are found (the presence of these characters would imply a file is using either PC or Mac-formatted end-of-line characters (\r\n, \r, respectively)
sub check_file_for_unix_friendliness {
	my $file_path = shift;
	my $original_terminator = $/;
	my $file_info;
	undef $/;
	open(INPUT, $file_path)
		or die "Can't open file for Unix-friendly check at $file_path!\n";
	$file_info = <INPUT>;
	close INPUT;
	if ($file_info =~ /\r/) {
		if ($file_info =~ /\r\n/) {
			die "Your file, $file_path, appears to be utilizing PC-formatted end-of-line characters.  Please reformat this file to use Unix end-of-line characters.\n";
		}
		elsif ($file_info =~ /\r/) {
			die "Your file, $file_path, appears to be utilizing Mac-formatted end-of-line characters.  Please reformat this file to use Unix end-of-line characters.\n";
		}
		else {
			die "Your file, $file_path, appears to be utilizing unrecognized end-of-line characters.  Please reformat this file to use Unix end-of-line characters.\n";
		}
	}
	$/ = $original_terminator;
	return 1;
}

sub usage {
	"Usage with required parameters: perl coproseq.pl -b <file> -g <directory> -i <file>\n",
			"\t-b : Barcode (BC) file.  Plain text file of all BC sequences in sample pool (1 per line).  All BCs must be the same length.\n",
			"\t-g : Directory containing all genomes possibly represented in sample in 'squashed' format.\n",
			"\t-i : Raw input file in SCARF format generated during primary data analysis.\n\n",
			"Optional parameters:\n",
			"\t-a : Keeps 'a'll intermediate files created during analysis in place, rather than deleting them once processing is complete.  Off by default.\n",
			"\t-l : Specifies 'l'ength of each read to use.  29 by default.  Remember to allow for barcode bases when choosing a value.\n",
			"\t-m : Specifies 'm'ismatch tolerance, from 0-2, describing how many mismatches Eland will allow when doing alignments.\n",
			"\t-o : Specifies location of 'o'utput directory where data will be stored.  Set to your prefix value (within current working directory) by default.\n",
			"\t-p : Adds the designated 'p'refix to each file created to make sorting/searching easier. Set to (system time)-(1-1000, randomly) by default.\n",
			"\t-s : Designates data as being 's'ingleplex, or non-barcoded.  Off by default.\n";
}

sub table_file_to_HoH {
	my $table_file_path = shift;
	my %HoH;
	open (INPUT, $table_file_path) || die "Error Can't open table file $table_file_path.\n";
	my $firstline = <INPUT>;		# First line should be tab-separated headers with empty first cell
	chomp $firstline;
	my @headers = split (/\t/, $firstline);
	shift @headers;	# Get rid of first element of headers array representing top left cell of table
	while (<INPUT>) {
		chomp;
		my @line = split(/\t/,$_);
		my $rowname = shift @line;
		for (my $i = 0; $i < scalar @line; $i++) {
			$HoH{$headers[$i]}{$rowname} = $line[$i];
		}
	}
	close INPUT;
	return \%HoH;
}

# Open barcode sequence/sample name file and return a hash describing its
# contents (key = barcode sequence, value = sample name)
sub parse_bc_file {
	my $bc_file_path = shift;
	my %bc_abbrev = barcodes::declare_barcodes();
	my %bc_occurrence;
	my %samples_by_bc;	# Key = barcode sequence, value = sample name

	# Confirm all EOL characters in .bc file are OK for Unix
	check_file_for_unix_friendliness($bc_file_path);

	open(BC_INPUT, $bc_file_path);
	foreach(<BC_INPUT>) {
		chomp;
		my @line = split(/\t/, $_);
		# If barcode is specified in .bc file as an abbreviation (e.g., BC1)
		if (defined	$bc_abbrev{$line[0]}) {
			$bc_occurrence{$bc_abbrev{$line[0]}}++;
			$samples_by_bc{$bc_abbrev{$line[0]}} = $line[1];
		}
		else {
			$bc_occurrence{$line[0]}++;
			$samples_by_bc{$line[0]} = $line[1];
		}
		if (@line != 2) {
			die "\nThe following line in file '$bc_file_path' has the wrong".
			" number of columns (i.e., ".scalar(@line).", rather than 2):\n\n".
			$_."\n\nNote that non-obvious characters such as tabs and ".
			"whitespace may be the problem if extra columns are not evident".
			"\n\n";
		}
	}
	close BC_INPUT;
	return \%samples_by_bc;
}

sub summarize_bcs {
	my $bc_array_ref = shift;
	print "\nNumber of barcodes: ".scalar(@$bc_array_ref)."\n";
	print "Barcode size: ".check_bc_size($bc_array_ref)."\n";
	print "Barcode sequences: ";
	for (my $k = 0; $k < (scalar(@$bc_array_ref) - 1); $k++) {
		print "$$bc_array_ref[$k], ";
	}
	print $barcodes[-1];
}

# Input: path to a file/symbolic link of sequences in either SCARF or FASTQ format (only two formats supported right now)
# Output: a string describing the file format detected (i.e., 'SCARF', 'FASTQ', or 'unknown')
# Detection is based on the first 4 lines of the sequence file
sub detect_seq_file_type {
	my $seq_file_path = shift;
	my $file_type;
	my $first_four_lines;
	open(INPUT, $seq_file_path) || die "Error: Can't open sequence input file $seq_file_path!\n\n";
	for (my $k = 0; $k < 4; $k++) {
		my $line = <INPUT>;
		$first_four_lines .= $line;
	}
	close INPUT;
	if ($first_four_lines =~ /@.+\n\w+\n\+\n.+\n/) {
		$file_type = 'FASTQ';
	}
	elsif ($first_four_lines =~ /(.+\:.+\:.+\:.+\:.+\:.+\:.+\n){4}/) {
		$file_type = 'SCARF';
	}
	else {
		$file_type = 'unknown';
	}
	return $file_type;
}

# I've updated the subroutine below to create files that will not end in an empty line, based on
# my belief that empty lines at the end of the outputs was causing Eland to throw scary error messages,
# but this output needs to be double-checked still to confirm that I haven't introduced any problems
sub parse_seq_file_by_bc {
	my ($seq_file_type, $ref_to_bc_array, $bcsize) = @_;
	unless (($seq_file_type eq 'FASTQ') || ($seq_file_type eq 'SCARF')) {
		die "Error: Your sequence file appears to be in a format that is ".
			"neither FASTQ nor SCARF.  These are the only two sequence file ".
			"formats currently supported by the COPRO-Seq workflow.\n\n";
	}

	print "\n\nParsing sequence file $input (type: ".$seq_file_type.") and ".
		"trimming sequence info by barcode...\n";
	my %barcodes_seen;
	my %filehandlehash;
	foreach my $a (@$ref_to_bc_array) {
		my $fh = IO::File->new(">$outputdir\/bcsortedseqs\/$prefix\_$a\.fas");
		$filehandlehash{$a} = $fh;
		$barcodes_seen{$a} = 1;
	}
	# SCARF example: HWI-EAS158:4:1:120:636#0/1:GAAATCCGTGTGGACAGCCGTCATCTGTTCCCGTCC:aabbbaaaaaa_a^aaa
	# Note quality scores are now stored as ASCII characters which must be decoded to get numeric scores
	open (SEQ_INPUT, $input) || die "\nCan't open $input\n\n";
	my %bc_tallies;	# key = barcode, value = # of occurrences
	my $readnumber=0;	# Will keep track of number of sequences in file
	if ($seq_file_type eq 'SCARF') {
		while (<SEQ_INPUT>) {
			# If not a blank line
			if (/.+/) {
				$readnumber++;
				my @SCARFline = split(/:/, $_); # Should generate an array with 7 elements (sequence is the 6th, [5])
				my $bc = substr($SCARFline[5], 0, $bcsize);
				$bc_tallies{$bc}++;	# Tally a hit for each n-mer barcode seen in SCARF file (even those not expected)
				if ($barcodes_seen{$bc}) {
					my $fh = $filehandlehash{$bc};
					# If statements below ensure that file ends without a blank line at the end (i.e., go to next line only when a new sequence to add to the output is found)
					if ($bc_tallies{$bc} == 1) {	# i.e., this is the first sequence for this barcode
						print $fh ">$SCARFline[1]\:$SCARFline[2]\:$SCARFline[3]\:$SCARFline[4]\n".substr($SCARFline[5],$bcsize,$readsize-$bcsize);
					}
					else {
						print $fh "\n>$SCARFline[1]\:$SCARFline[2]\:$SCARFline[3]\:$SCARFline[4]\n".substr($SCARFline[5],$bcsize,$readsize-$bcsize);				
					}
				}
			}
			# (Else: line is an empty line and should be ignored)
		}
	}
	elsif ($seq_file_type eq 'FASTQ') {
		# Parse lines in groups of 4
		my $line_count = 0;
		my $latest_entry;
		while (<SEQ_INPUT>) {
			# If not a blank line
			if (/.+/) {
				$line_count++;
				$latest_entry .= $_;
				# If we've just grabbed the last/4th line of a set...
				if ($line_count % 4 == 0) {
					$readnumber++;
					my @FASTQentry = split(/\n/, $latest_entry);
					my $header = $FASTQentry[0];
					$header =~ s/^@//;
					# Convert whitespace to underscores (whitespace can wreak
					# havoc on the ability of other programs, e.g. Eland, to 
					# parse FASTA files correctly)
					$header =~ s/\s/_/;
					my $bc = substr($FASTQentry[1], 0, $bcsize);
					$bc_tallies{$bc}++;
					if ($barcodes_seen{$bc}) {
						my $fh = $filehandlehash{$bc};
						# If statements below ensure that file ends without a blank line at the end (i.e., go to next line only when a new sequence to add to the output is found)
						if ($bc_tallies{$bc} == 1) {	# i.e., this is the first sequence for this barcode
							print $fh ">$header\n".substr($FASTQentry[1],$bcsize,$readsize-$bcsize);	
						}
						else {
							print $fh "\n>$header\n".substr($FASTQentry[1],$bcsize,$readsize-$bcsize);	
						}
					}
					$latest_entry = "";
				}
			}		
		}
		# LAST ENTRY RECORDED?
	}
	close SEQ_INPUT;
	foreach(@barcodes) {
		close($filehandlehash{$_});
	}
	print "\nTotal number of reads sorted: $readnumber\n\n";
	return (\%bc_tallies, $readnumber);
}

sub report_bc_tallies {
	my ($ref_to_tallies_hash, $readnumber, $ref_to_bc_array) = @_;
	# Make a barcodes lookup hash
	my %bc_lookup_hash;
	foreach(@$ref_to_bc_array) {
		$bc_lookup_hash{$_} = 1;
	}
	my %bc_tallies = %{$ref_to_tallies_hash};	 # Dereference
	
	# Report statistics for each barcode-specific set of sequences to an output file (include barcode sequence, # occurrences, and % of all reads represented by that barcode)
	open (BCSTATS, ">$outputdir\/seqcountsbybc\/$prefix"."\.countsbybc") || die "Can't open $outputdir\/seqcountsbybc\/$prefix"."\.countsbybc";

	# Write results (sorted by most -> least # of hits per barcode) to output file
	if($tags) {	print BCSTATS "<$prefix>\n"; }
	foreach my $key (sort{$bc_tallies{$b} <=> $bc_tallies{$a}} keys %bc_tallies) {
		# If current barcode is one of the expected barcodes, print an asterisk in the fourth column to make spotting the expected barcodes easier in the results
		if (defined $bc_lookup_hash{$key}) {
			# Syntax on this printf statement should be double-checked
			printf BCSTATS "%s\t%d\t%7.4f%%%s\n", $key, $bc_tallies{$key}, ($bc_tallies{$key}/$readnumber)*100, " *";
		}
		# If not an expected barcode, just print its sequence, # of times it appeared, and its percentage of total reads
		else {
			printf BCSTATS "%s\t%d\t%7.4f%%\n", $key, $bc_tallies{$key}, ($bc_tallies{$key}/$readnumber)*100;
		}
	}
	close BCSTATS;
}

#	&check_options(\$bcseqs, \$genomesdir, \$input, \$readsize, \$mismatches_allowed, \$outputdir, \$prefix);
# sub check_options {
# 	my $bcseqs_ref = shift;
# 	my $genomesdir_ref = shift;
# 	my $input_ref = shift;
# 	my $readsize_ref = shift;
# 	my $mismatches_allowed_ref = shift;
# 	my $outputdir_ref = shift;
# 	my $prefix_ref = shift;
# 
# 	# Make sure input file exists, can be opened, and is valid
# 	# Make sure bcseqs file is valid
# 	# Make sure genomesdir exists and is valid
# 	# Make sure $readsize is a valid base # range
# 	# Make sure $mismatches_allowed is 0, 1 or 2
# 	# Make sure $outputdir does not already exist
# }
