#!/usr/bin/perl

# Last Modified 04/28/2010
# Author: Nate McNulty
# Purpose: To calculate, given two or more genomes and one or more k-mer length(s), the informative genome size (IGS) of each
#	genome
# Inputs:
#	Directory containing all genomes in FASTA nucleotide format - .fa/.fas/.fna - for use in making k-mers
#	Directory containing squashed genomes for use in Eland alignment of k-mers to all references in the analysis
#	Option to squash FASTA genomes that do not already have .2bpb and .vld files in the squashed genomes directory	
#	k-mer size(s) to consider (i.e. the read length(s) one intends to use in downstream applications such as CoPro-seq)
#		If no k-mer size(s) is/are given, the software will calculate the IGS for a large range of k-mer sizes (10-32bp)
# Outputs:
#	Squashed versions of each genome in the squashed genome directory (if -p option invoked)
#	k-mer sequence files for each genome (these can ultimately be deleted)
#	An Eland output file for each k-mer sequence file considered
#	Tab-delimited table (row: genome name, column: k-mer size) showing IGS values
#	Tab-delimited table providing genome uniqueness %age for a given k-mer size (i.e. % of all k-mers from that 
#		genome that show up only once among all genomes considered)

# NOTE: For now, the software is best for finished (closed) genomes.  As implemented, draft genomes must first be collapsed
#	to single contigs by concatenation, which will generate a very small number of false reads from the potential read space
#	of all k-mers that a genome can produce (this will occur at junctions between contigs that are concatenated).  It is a very minor
#	problem that should eventually be fixed by passing the makekmers function a hash of contigs that can be checked one key/value pair
#	at a time while building one large list of all k-mers possible

# Improvements:
#	Check that reference genomes meet all of the requirements for the underlying software:
#		References should be single-entry nucleotide FASTA files
#		References must have ONLY A,C,G,T, or N bases, or a segmentation fault will occur during mapping (there is not yet support for IUPAC ambiguity codes)
#			If incorrect bases are found, what's the best course of action?  Replace them with N's?  Kill program and report error?

use warnings;
use strict;
use Getopt::Long;
use Cwd;
use Math::Round;
use Bio::Tools::GuessSeqFormat;

my $squash_program_location = 'squashGenome';

my $cleanup;					# Default: do not clean up afterwards (i.e. leave intermediate, batch and squashed files in place)
my $fastadir = getcwd();		# Default: genomes are in directory from which script was called
my $outputdir = getcwd();		# Default: output will be written to directory from which script was called
my $squashdir = getcwd();		# Default: will look for all squashed genomes in directory from which script was called
my $kmersizesinput;				# Default: k-mer sizes are undefined
my $squash;						# Default: no genomes in $fastadir will be squashed (all are assumed to be squashed already and present in $squashdir)
GetOptions (
				'c|cleanup'				=> \$cleanup,
				'f|fastadir=s'			=> \$fastadir,
				'k|kmersizesinput=s'	=> \$kmersizesinput,		# Multiple values are split by "," below
				'o|output=s'			=> \$outputdir,
				'p|prepare'				=> \$squash,
				's|squashdir=s'			=> \$squashdir
) or die usage();

# Global constants for easy manipulation
my $MEMREQ = "8G";			# Memory that will be requested for eland alignment jobs
my $JOB_CHECK_PAUSE = 10;	# Number of seconds to wait between checking status of all Eland jobs; Err on the side of keeping this a little long (the method
							#	used for checking job status is kind of clunky and will hammer the file server if you try to do it too often for too many files)

my @kmersizes;
if ($kmersizesinput) {
	@kmersizes = split(/,/, $kmersizesinput);
}
else {	
	# If no k-mer sizes specified, do calculations for all k-mers from 15-32bp (these are sizes for which eland executables are available)
	@kmersizes = (15..32);	
}	

#Create an array of strings representing the paths of only the genome sequence files in the directory given

my @allfilenames = <$fastadir/*>;
my $fasta_files_ref = get_fasta_filenames(\@allfilenames);
my @genomefilepaths = @$fasta_files_ref;

print "\nFASTA genome files that will be considered in IGS calculations:\n";
foreach (@genomefilepaths) {
	print "$_\n";
}
print "\n";

if ($squash) {	
	system("$squash_program_location $squashdir @genomefilepaths");	
	print "\n";
}

# At this point all squashed genomes should be available regardless of options chosen at beginning of script

# For each genome:
#	For each k-mer size:
# 		Generate k-mer FASTA file needed for Eland (call makekmers sub)
#		Add that k-mers file to a batch request for eland alignment
# 		Run Eland on batch file with proper system requirements against all squashed genomes in directory
#		Parse each Eland results file to determine % of all k-mers that hit uniquely to one spot with no mismatches (U0)
#		Multiply that % by the genome's size (to give you IGS)
#		Store results in a hash of hashes
#		Print the sorted contents of the hash of hashes

# Create all k-mers files that will be needed (note, this could take up quite a bit of disk space)
my @kmerfilepaths;

foreach my $g (@genomefilepaths) {

	# Get genome sequence
	my $seq;
	open (GENOMESEQ, "<$g") || die "Can't open $g\n";
	<GENOMESEQ>;					# Skip header sequence
	while (<GENOMESEQ>) {
		chomp;
		$_ =~ s/\s//;
		$seq .= $_;	# Next line is full-length genome
	}
	close GENOMESEQ;

	# Grab genome name from file path...first get everything after the last "/", then remove extension
	my $genomename = "";
	if ($g =~ m/.*\/(.*)$/) { # Grab filename from path
		$genomename = $1;
	}
	$genomename =~ s/(\..*)$//;	# Remove extension (everything after the last ".")

	foreach my $k (@kmersizes) {
		my $kmerfilepath = $outputdir."\/".$genomename."."."$k"."bp.kmers";
		push (@kmerfilepaths, $kmerfilepath);		# Keep track of all k-mer filenames for later
		print "Creating $kmerfilepath"."...\n";
		my $fasta_output_ref = makekmers($seq, $k);
		open (KMEROUTPUT, ">$kmerfilepath") || die "Can't open $kmerfilepath\n\n";
		print KMEROUTPUT $$fasta_output_ref;
		close KMEROUTPUT;
	}
}
print "\n";

# Create batch file of all Eland jobs using k-mer filenames determined above
my @elandoutpaths;	# To keep track of all eland output filenames for later
open (ELANDBATCH, ">$outputdir\/elandIGSsubmissions.jobs");	# Modify to include unique timestamp
foreach my $n (@kmerfilepaths) {
 	my $ksize;
 	if ($n =~ m/.*\.([0-9]+)bp.kmers$/) {
 		$ksize = $1;
 	}
 	my $outfilestring = $n . ".elandout";
	push (@elandoutpaths, $outfilestring);
 	print ELANDBATCH "/srv/cgs/local/gapipeline/latest/bin/eland_$ksize $n $squashdir $outfilestring\n";
}
close ELANDBATCH;

# Submit batch file to nq
# Eventually refine this by checking for needed memory requirements by looking at total size of all squashed files
# Need some kind of status checking loop to make sure all alignments are finished - recycle coproseq code?
system("nq $outputdir\/elandIGSsubmissions.jobs | qsub -l h_vmem=$MEMREQ -l arch=lx26-amd64");

# MONITOR ELAND ALIGNMENT JOBS UNTIL ALL ARE FINISHED
# Major problem with this approach is that jobs will appear to be running indefinitely if the output file is created before the program dies
#	Could instead rely on parsing qstat output for the user, continually rechecking until there are no jobs running AND all of the files created have contents
#	To do this, would need the original jobID# (how would I know which # this is if there are multiple jobs running?), username, and ???
# Eventually collapse this into a sub that is passed references to @elandoutpaths and the # elements in @kmerfilepaths
my $notstarted = 0;
my $running = 0;
my $completed = 0;
my $starttime = scalar localtime();
print "\n";
while ($completed < scalar(@kmerfilepaths)) {	# While there are still Eland output files that aren't made/done yet...
	$notstarted = $running = $completed = 0;	# Reset variables so that all files can be checked again
	foreach (@elandoutpaths) {
		my $statuscode = &check_for_file($_);
		if ($statuscode == 0) 		{		$notstarted++;	}
		elsif ($statuscode == 1) 	{		$running++; 	}
		elsif ($statuscode == 2)	{		$completed++;	}
	}
	print	"Status of jobs as of " . scalar localtime() . " (started on $starttime):\n",
			"Total:\t" . scalar(@kmerfilepaths) . "\t",
			"Complete:\t$completed\t",
			"Running:\t$running\t",
			"Not yet started:\t$notstarted\n";
	print "(Started checking for Eland results on ". scalar localtime() . ")\n\n";		# Change code to figure out how long you've been waiting instead
	sleep $JOB_CHECK_PAUSE;
}
print "All alignments complete.\n\n";

my %allIGS= ();				# HoH; Key1 = genome, Key2 = k-mer size, value = IGS
my %allpercentages = ();	
# Load up hash of hashes with results using genome name as key1, k-mer size as key2, and IGS in Mb as value
# Consider altering the block below to be a sub instead where references to the two hashes are passed (hashes will be modified within the sub)
foreach (@elandoutpaths) {
	# Get info that will become keys for HoH
	my ($organism, $kmersize);
	# Grab everything up to "bp.kmers.elandout" (this will exist for every Eland output file created)
	if ($_ =~ m/^.*\/(.*\.[0-9]+)bp\.kmers\.elandout$/) {
		my @info = split(/\./,$1);
		$organism = $info[0];
		$kmersize = $info[1];
	}

	print "Calculating IGS for $organism at ${kmersize}bp...\n";

	# Parse every Eland results file and calculate IGS
	my ($all, $hit) = &parse_eland_by_code("U0", $_);	# Returns total lines checked, lines matching correct code for .elandout file
	my $full_length = ($all + ($kmersize - 1)) / 1000000;	# In Mb
	my $percentage = $hit / $all;
	my $IGS = $full_length * $percentage;
	$allIGS{$organism}{$kmersize} = $IGS;
	$allpercentages{$organism}{$kmersize} = $percentage;
}

open PERCENTTABLE, ">$outputdir/genome_uniqueness.table" or die "Error: Can't open genome_uniqueness.table output file!\n\n";
print_sorted_HoH_table(\%allpercentages, *PERCENTTABLE);
close PERCENTTABLE;

open IGSTABLE, ">$outputdir/IGS.table" or die "Error: Can't open IGS.table output file!\n\n";
print_sorted_HoH_table(\%allIGS, *IGSTABLE);  
close IGSTABLE;

if ($cleanup) {	cleanup(); }		# Remove intermediate files specified by user, and elandjobs error/output files (once all Eland jobs are done)

exit;

# print_sorted_HoH_table (\%, <FH>); Input: a reference to a HoH and, optionally, a file handle (will print to STDOUT if no file handle is specified)
sub print_sorted_HoH_table {
	my %HoH=();
	my $out;
	if (@_ == 2) {			# i.e. Passed a hashref and a file handle
		%HoH = %{(shift)};	# Dereference to make working with HoH easier
		$out = shift;
	}
	elsif (@_ == 1) {
		%HoH = %{(shift)};	# Dereference to make working with HoH easier
		$out = *STDOUT;
	}
	
	# Print column headers line (using a randomly selected key1 to get k-mer headers?):	"Genome\tSize#1\tSize#2\tSize#3...\n";
	my @sorted_organisms = sort {lc($a) cmp lc($b)} keys %HoH;
	my $example = shift(@sorted_organisms);
	unshift(@sorted_organisms,$example);	# Return the first random element so I can use @sorted_organisms again later
	my @sorted_sizes = sort {$a <=> $b} keys %{$HoH{$example}};
	print $out "Organism\t";
 	my $numsizes = scalar @sorted_sizes;
 	for (my $k = 0; $k < ($numsizes-1); $k++) { print $out "$sorted_sizes[$k]\t"; }
 	print $out "$sorted_sizes[($numsizes-1)]\n";

	# Print contents of HoH, sorting keys along the way
	foreach my $species (@sorted_organisms) {
		print $out "$species";
		foreach my $size (@sorted_sizes) { 
			print $out "\t".$HoH{$species}{$size}; 
#			print $out "\t".nearest(.001, $HoH{$species}{$size}); 	# Rounding eventually causes problems 
#																	(i.e. IGS-adjusted percentages only sum properly out to 2nd digit, e.g. 0.998)
		}
		print $out "\n";
	}
}

# sub makekmers($genome,$kmersize), where $genome is a complete genome sequence, 5'->3'
sub makekmers {		
	my ($genome,$grabsize) = @_;
	my $kmerfasta = "";
	my $genomelength = length($genome);
	# Line below determines the last position along genome length for which you can still grab a k-mer of the desired size
	my $stopposition = $genomelength - $grabsize;
	
	# Grab each possible k-mer from the genome string using a 1bp sliding window and write the results to a giant string that will be returned
	my $kmerseq;
	for (my $j= 0; $j <= $stopposition; $j++){		# j = Position of sliding window
		$kmerseq = substr ($genome, $j, $grabsize);		# Grab a sequence of length defined by user starting at the position of the sliding window
		# Chose to use one big string below in hopes of speeding up code - could consider a different data structure if this doesn't work well
		$kmerfasta .= ">Bases_".($j+1)."_to_".($j + $grabsize)."\n$kmerseq\n";
	}
	return \$kmerfasta;
}

sub usage {
	print "Could not initiate script.  There appears to be something wrong with your parameters.\n";
	print "Usage: perl calcIGS.pl\n",
			"\t-f : Input directory containing all nucleotide FASTA-formatted genomes you wish to include in the analysis\n",
			"\t-c : Clean up intermediate files used in IGS calculations to reduce clutter\n",
			"\t-k : k-mer sizes to include in analysis; May be one or several comma-separated values (e.g. '-k 20' or '-k 10,12,14,16,18')\n",
			"\t-o : Path of output directory where intermediate files and results tables should be written\n",
			"\t-p : Prepare (squash) genomes in FASTA genome directory (passed with -f) and store results in squashed genome directory\n",
			"\t\t(passed with -s).  If not turned on, the squashed genome files (.2bpb and .vld) for each species must be included in the\n",
			"\t\tsquashed genome directory before the script is called.\n";
			"\t-s : Squashed genome directory where pre-squashed files reside or where newly-squashed files should be deposited (when -p \n",
			"\t\toption is invoked)\n";
}

sub checkparameters {
# Rules:
#	Genomes directory must be specified and a directory that can be opened
#	Genomes directory must have at least one (two?) .fa/.fas/.fna files in it
#	If specific k-mer sizes are specified:
#		They must all be integers 15<=x<=32
#		Redundancies should be filtered out
#	If squash option is not selected, there must be at least one squashed genome present in the directory
#		Start with list of all files in directory ending in .2bpb.  Until the root of one of those files is found to also have a .vld file, 
#			continue looking.  Escape loop when a matching file is found and set error code to OK; otherwise, error code should stay at -1
#	Option to remove squashed files after program runs can only be turned on if option to squash genomes is turned on
}

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

sub cleanup {
	# Remove k-mer sequence files if specified by user
	system("rm $outputdir/*.kmers");
	# Remove elandjobs files and batch jobs file (no matter what)
	system("rm elandIGSsubmissions*");
	# Remove .elandout files?
	system("rm $outputdir/*.elandout");
}

# parse_eland_by_code ($code, $filepath)
sub parse_eland_by_code {
	my $code = shift;			# Should be a text string code (e.g. U0, U1, R0, NM, etc.)
	my $filepath = shift;
	my ($all,$hits) = (0);
	open (ELANDFILE, "$filepath") || die "Can't open $filepath during parsing of Eland results.\n\n";
	while (<ELANDFILE>) {
		$all++;
		my @line = split(/\t/,$_);
		if ($line[2] eq $code) { $hits++ };
	}
	close ELANDFILE;
	return ($all,$hits);
}

sub get_fasta_filenames {
	my $allfiles_ref = shift;
	my @fastafiles;
	foreach(@$allfiles_ref) {
		if (-f $_) {	# If this list member is a file (and not a directory, e.g.)
			# Line below guesses the format of a flat file, given a filename
			unless (($_ =~ /\.2bpb$/) || ($_ =~ /\.vld$/)) {						# Ignore any previously squashed genome outputs...
				my $guesser = new Bio::Tools::GuessSeqFormat( -file => $_ );
				if ($guesser->guess) {	# If 
					if ($guesser->guess eq 'fasta') {	# If the format appears to be fasta...
						push(@fastafiles, $_);
					}
					else {
						print "Ignoring file $_ in the FASTA genomes directory.  It does not appear to a FASTA-formatted file, and may be of type '". $guesser->guess . "'.\n";
					}
				}
				else {
					print "Ignoring file $_ in the FASTA genomes directory.  It does not appear to be a FASTA-formatted file and is of unknown type.\n";
				}
			}
		}
	}
	return \@fastafiles;
}
