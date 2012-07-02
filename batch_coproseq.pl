#!/usr/bin/perl
# Author: Nathan McNulty
# Last updated: June 2012
# Acknowledgements: J. Faith, for the re-use of some of his RNA-Seq pipeline
#  code and for his clever use of Google spreadsheets as a front-end for the
#  workflow

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use DBI;
use Bio::DB::GenBank;
# FindBin variable '$FindBin::Bin' stores the path to batch_coproseq.pl
use FindBin;	
# Provide location of barcodes.pm and genomecodes.pm
use lib "$FindBin::Bin/lib";
use barcodes;
use genomecodes;
use Cwd;

# Variables for interacting with the microbialomics server (-ncbi not invoked)
my $db_name = "microbialomics_npm_mw";
my $host = "hamlet";
my $pass = "reader";
my $user = "reader";
my $dbh = DBI->connect("DBI:mysql:$db_name:$host",$user,$pass)
	or die "can't open database: $DBI::errstr\n";

# Locations of needed programs/scripts
my $squash_exec_path =
	"/srv/cgs/local/gapipeline/GAPipeline-1.5.0/bin/squashGenome";
my $summarize_tables_script_path =
	"$FindBin::Bin/summarize_tables.pl";
my $IGS_calc_script_path =
	"$FindBin::Bin/IGScalc.pl";
my $coproseq_script_path =
	"$FindBin::Bin/coproseq.pl";

# Declarations for folder/file names where info will be stored
my $basedir = cwd;
my $genomesdir = "genomes";
my $summariesdir = "summaries";
my $GEOdir = "GEO";
my $squashedgenomesdir = "genomes/squashed";
my $project_data_path = "project.info";
my $getdata_file_path = "getdata.sh";
my $align_file_path = "align.sh";
my $cleanup_file_path = "cleanup.sh";
my $IGS_calc_file_path = "calcIGS.sh";
my $summarize_file_path = "summarize.sh";

# Options declarations
my $google_key = "0AhsSO_Vep9tqdFJFTjh1S0UyZFRmeFBHamdKU3I5RHc";
my ($allfiles, $GEO, $group, $ncbi, $mapping_file_path, $IGS_table_file, $single_cpu);
my $mismatches_allowed = 0;
my $readsize = 25;	# Length AFTER trimming off barcode

GetOptions(
	"a|allfiles"				=> \$allfiles,
	"e|errors=i"				=> \$mismatches_allowed,
	"g|group=s"					=> \$group,
	"G|GEO"						=> \$GEO,
	"i|igs=s"					=> \$IGS_table_file,
	"k|key=s"					=> \$google_key,
	'l|length=i'				=> \$readsize,
	"m|mapping=s"				=> \$mapping_file_path,
	# ncbi option below is NOT working for draft genomes
	"n|ncbi"					=> \$ncbi,
	"o|output=s"				=> \$basedir,
	"s|single_cpu"				=> \$single_cpu
) or die usage();

# Confirm all necessary options were specified correctly
check_options();

# Confirm that mapping file is formatted with the correct EOL characters
check_file_for_unix_friendliness($mapping_file_path);

# Process value passed to group (-g) option
my $groups_hash_ref = parse_groups($group);

#===============================================================================
# Google spreadsheet processing
#===============================================================================

# Copy Google spreadsheet to local text file
#download_table("https://docs.google.com/spreadsheet/ccc?key=$google_key&output=txt&gid=0", $project_data_path);
download_table("https://docs.google.com/spreadsheet/ccc?key=$google_key&output=txt", $project_data_path);
#print "Downloaded spreadsheet from location:\nhttps://docs.google.com/spreadsheet/ccc?key=$google_key&output=txt&gid=0\n";

# Load contents of 'project.info' into hierarchy of variables
my $project_data_hash = table_to_hash($project_data_path);
# Filter on group(s) specified at startup using the -g option
my $filtered_data_hash = filter_by_group($project_data_hash, $groups_hash_ref);
# See 'Appendix A' below for overview of $filtered_data_hash organization
check_for_absent_groups($filtered_data_hash, $groups_hash_ref);

#===============================================================================
# Genome download and squashing (reference database creation for aligner)
#===============================================================================
my $species_list_ref = get_species_list($filtered_data_hash);
my @species_names = prepare_references($species_list_ref);
squash_genomes($squashedgenomesdir);
print "\n";

# Prepare 'getdata.sh'
make_get_data_file($filtered_data_hash, $getdata_file_path);

# Load up mapping hash with contents of mapping file
# $mapping_hash = reference to HoH{pool}{sample}=barcode
my $mapping_hash = load_mapping_file($mapping_file_path);

# Prepare '.bc' files
my %bclength_by_pool;	# Use to keep track of barcode length in each pool
# For each lane that's part of the analysis group...
for my $p (@$filtered_data_hash) {
	my $filteredsamples = get_filtered_sample_list($p, $mapping_hash);
	$bclength_by_pool{$p->{pool}} = make_bc_file($mapping_hash, $p->{pool}, 
		$filteredsamples, $p->{machine}.'_'.$p->{run}.'_'.$p->{lane}."\.bc");
}

# Prepare '.bc_mod' files and 'split_scarfs.sh' if -G switch was turned on
# Makes preparation of barcode-split SCARF files easier for the user
# Note: split_scarfs.sh must be run AFTER getdata.sh
if ($GEO) {	
	# Declare variables for storing running tally of commands that will be
	# the contents of GEO-related .sh files
	my ($split_sh_commands, $move_sh_commands, $make_table_sh_commands);
	my $a;
	for my $p (@$filtered_data_hash) {
		my $filteredsamples = get_filtered_sample_list($p, $mapping_hash);
		my $new_bc_file = $p->{run}."_".$p->{lane}."\.bc_mod";
		my $scarf_to_split = "../".$p->{run}."_".$p->{lane}.".scarf";
		make_modified_bc_file($p->{run}, $p->{lane}, 			
			\%{$$mapping_hash{$p->{pool}}}, $filteredsamples,
			"GEO\/$new_bc_file");
		$split_sh_commands .= 
			"perl ../../split_barcodes.pl $new_bc_file $scarf_to_split\n";
		foreach(@$filteredsamples) {
			$move_sh_commands .= 
				"cp ../hitratios/".$p->{machine}.'_'.$p->{run}.'_'.$p->{lane}.'_'.$_.'_' .
				"specieshits_*bp_*MM.hitratios run".$p->{run} .
				'_lane'.$p->{lane}."_".$$mapping_hash{$p->{pool}}{$_} .
				'_'.$_.".hitratios\n";
		}
	}
	write_GEO_sh(">./GEO/split_scarfs.sh", $split_sh_commands);
	write_GEO_sh(">./GEO/move_hits.sh", $move_sh_commands);
	$make_table_sh_commands = "perl ../../prepare_geo_table.pl -g $group " .
		"-m ../$mapping_file_path -s ../$project_data_path -o GEO_table.txt\n";
	write_GEO_sh(">./GEO/prepare_GEO_table.sh", $make_table_sh_commands);
	write_GEO_readme(">./GEO/README.txt");
}

# Prepare 'align.sh'
make_align_file($filtered_data_hash, \%bclength_by_pool, $align_file_path);

# Prepare/validate 'IGS.table'
$IGS_table_file = make_igs_table_file(\@species_names);

# Prepare 'summarize.sh'
make_summarize_file($filtered_data_hash, $summarize_file_path);

# Prepare 'cleanup.sh'
make_cleanup_file($filtered_data_hash, $cleanup_file_path);

exit;



sub make_igs_table_file {
	my $names_ref = shift;
	my $IGS_table_path;
	# If user specified a pre-made IGS table...
	if ($IGS_table_file) {
		# If it fails to validate...
		if (! validate_IGS_table($IGS_table_file, $readsize, $names_ref)) {
			# Calculate a new IGS table
			make_IGS_calc_file($IGS_calc_file_path);
			# Set location of IGS table to path of newly created file
			$IGS_table_path = "$FindBin::Bin/genomes/IGS/IGS.table";
		}
	}
	# If no IGS table was specified...
	else {
		make_IGS_calc_file($IGS_calc_file_path);
		$IGS_table_path = "$FindBin::Bin/genomes/IGS/IGS.table";
	}
	return $IGS_table_path;
}

# Inputs: file path for IGS table, list of species that should be present, k-mer size that should be present
sub validate_IGS_table {
	my ($table_file_path, $kmerlength, $species_list_ref) = @_;
	my $HoH_ref = table_file_to_HoH($table_file_path);
	my %HoH = %$HoH_ref;	# Dereference to make things easier
	print "\nValidating IGS table...\n";
	# Check for correct read size
	if (! defined $HoH{$readsize}) {
		print 	"Warning: Could not detect the read size passed at startup in your IGS table.  A ",
				"new IGS table will have to be calculated by executing calcIGS.sh.\n";
		return 0;
	}
	print "Read size specified at startup detected in IGS table.\n";
	# Check that all species needed are present
	foreach(@$species_list_ref) {
		if(! defined $HoH{$readsize}{$_}) {
			print	"The IGS value for at least one genome ($_) is missing from your IGS table.  A new ",
					"IGS table will have to be calculated by executing calcIGS.sh.\n";
			return 0;
		}
	}
	print "All species specified at startup detected in IGS table.\n";
	print "User-supplied IGS table validated.  No new IGS table is required.\n\n";
	return 1;
}


# Input: file path for input table (table must have header rows for each column with first column being row names)
#	e.g.      |  A  |  A  |  A  |  A
#		-------------------------------
#		   R  |  #  |  #  |  #  |  #
#		-------------------------------
#		   R  |  #  |  #  |  #  |  #
#		-------------------------------
#		   R  |  #  |  #  |  #  |  #  
# Output: Reference to HoH with structure $HoH{column header}{row name} = table value
sub table_file_to_HoH {
	my $table_file_path = shift;
	my %HoH;
	open (INPUT, $table_file_path) || die "ERROR: Can't open table file $table_file_path.\n";
	my $firstline = <INPUT>;		# First line should be tab-separated headers with empty first cell
	chomp $firstline;
	my @headers = split (/\t/, $firstline);
	shift @headers;	# Get rid of first (empty) element of headers array representing top left cell of table
	while (<INPUT>) {
		my @line = split(/\t/,$_);
		my $rowname = shift @line;
		for (my $i = 0; $i < scalar @line; $i++) {
			$HoH{$headers[$i]}{$rowname} = $line[$i];
		}
	}
	close INPUT;
	return \%HoH;
}

# Input: reference to HoH
# Output: printed HoH table where headers and rownames are sorted
sub print_HoH_from_table_file {
	my $HoH_ref = shift;
	my %HoH = %$HoH_ref;		# Dereference to make a little easier to work with
	my @sortedheaders = sort {lc($a) cmp lc($b)} (keys %HoH);
	my @sortedrownames = sort {lc($a) cmp lc ($b)} (keys %{$HoH{$sortedheaders[0]}});
	# Choice of $sortedheaders[0] above is arbitrary; all headers should have the same associated keys (row names)
	foreach (@sortedheaders) {	print "\t$_";	}
	print "\n";
	my $num_of_columns = scalar @sortedheaders;
	foreach my $rowname (@sortedrownames) {
		print "$rowname\t";
		for (my $i=0; $i < ($num_of_columns - 1); $i++) {				# For all but the last element of the row...
			print "$HoH{$sortedheaders[$i]}{$rowname}\t";
		}
		print "$HoH{$sortedheaders[$num_of_columns - 1]}{$rowname}\n";	# Print last column with \n rather than \t

	}
}

sub make_summarize_file {
	my ($spreadsheet_hash, $filepath) = @_;
	open (SUMMARIZE, ">$filepath") || die "ERROR: Can't open $filepath!\n";
	print SUMMARIZE "mkdir $summariesdir\n";
	
	# Summarize all .countsbybc files; use spreadsheet hash to make list of all .countsbybc files
	print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 2 -i seqcountsbybc > $summariesdir/counts.bcdist\n";
	print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 3 -i seqcountsbybc > $summariesdir/percent.bcdist\n";

	# Summarize all .mappingstats files
	print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 2 -i mappingstats > $summariesdir/counts.mappingstats\n";
	print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 3 -i mappingstats > $summariesdir/percent.mappingstats\n";

	# Summarize all .hitratios files
	print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 2 -i hitratios > $summariesdir/raw_counts.profile\n";
	if($IGS_table_file) {
		print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 3 -i hitratios > $summariesdir/norm_counts.profile\n";
		print SUMMARIZE "perl $summarize_tables_script_path -r 1 -v 4 -i hitratios > $summariesdir/norm_percent.profile\n";
	}
	close SUMMARIZE;
}

sub make_IGS_calc_file {
	my $filepath = shift;
	open (IGS, ">$filepath") || die "ERROR: Can't open $filepath!\n";
	print IGS "mkdir genomes/IGS\n";
	print IGS "perl $IGS_calc_script_path -f genomes -s genomes/squashed -o genomes/IGS -k $readsize -c\n";
	close IGS;
}

sub make_cleanup_file {
	my ($spreadsheet_hash, $filepath) = @_;
	open (CLEANUP, ">$filepath") || die "ERROR: Can't open $filepath!\n";
	print CLEANUP "rm $project_data_path\n";
	print CLEANUP "rm $getdata_file_path\n";
	print CLEANUP "rm $align_file_path\n";
	print CLEANUP "rm $cleanup_file_path\n";
	print CLEANUP "if [ -f $align_file_path\* ];\nthen\nrm $align_file_path*\nfi\n";
	print CLEANUP "if [ -f $IGS_calc_file_path ];\nthen\nrm $IGS_calc_file_path\nfi\n";
 	print CLEANUP "rm $summarize_file_path\n";
#	print CLEANUP "rm -rf genomes\n";			# Main this folder in case user has manually included some of their own genomes
	for my $p (@$spreadsheet_hash) {
		my $prefix = $p->{machine}.'_'.$p->{run}.'_'.$p->{lane};
		print CLEANUP "rm $prefix.bc\n";
		print CLEANUP "rm $prefix.seq\n";
		print CLEANUP "rm -rf $prefix\n";
	}
#	print CLEANUP "if [ -f elandjobs* ];\nthen\nrm elandjobs*\nfi\n";	# Doesn't work because -f returns multiple filenames...need to rework with ls
	print CLEANUP "rm -rf bcsortedseqs\n";
	print CLEANUP "rm -rf elandresults\n";
#	print CLEANUP "rm -rf hitratios\n";			# Maintain this folder so that processed data files inside can be deposited with GEO if necessary
	print CLEANUP "rm -rf seqcountsbybc\n";
	print CLEANUP "rm -rf mappingstats\n";
	print CLEANUP "rm -rf NM\n";
	close CLEANUP;
	return 1;
}

sub make_align_file {
# perl ~/scripts/coproseq/coproseq.pl -i 207_5_SCARF -b 207_5_bcs -g ~/SquashedGenomes/NM601 -p 207_5 -m 0 -l 29
	my ($spreadsheet_hash, $bclengths_hash_ref, $filepath) = @_;
	my $cpu;
	# Set default name/location for IGS table file if it has not been specified
	# (IGS.table will always have same name/location when calcIGS.sh is run)
	unless (defined $IGS_table_file) {
		$IGS_table_file = "$basedir/genomes/IGS/IGS.table";
	}
	if ($single_cpu) {
		$cpu = ' -s';
	}
	else {
		# Leave no flag when building alignment job task
		$cpu = '';
	}
	open (ALIGNFILE, ">$filepath") || die "ERROR: Can't open $filepath!\n";
	for my $p (@$spreadsheet_hash) {
		my $prefix = $p->{machine}.'_'.$p->{run}.'_'.$p->{lane};
		print ALIGNFILE "perl $coproseq_script_path -i $prefix.seq -b $prefix.bc -g $basedir/genomes/squashed" .
						" -p $prefix -m $mismatches_allowed -o $basedir -n $IGS_table_file -t -l " . ($readsize+$$bclengths_hash_ref{$p->{pool}}) . 
						$cpu."\n";
	}
	close ALIGNFILE;
	return 1;
}

# Returns the length of all barcodes from a list, or zero if they are of different lengths
# Dependent on barcodes.pm
sub get_bc_length {
	my @barcodes = @_;
	my $bclength;
	foreach(@barcodes) {
		if (!$bclength) {
			$bclength = length($_); 
		}
		elsif ($bclength != length($_)) {
			return 0;
		}
	}
	return $bclength;
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

# Input: mapping file path
# Output: reference to a hash of hashes containing all poolnames->samples->barcodes (barcodes being literal barcode sequences)
sub load_mapping_file {
 	my $file_path = shift;
 	my %barcode_abbrev_to_seq = barcodes::declare_barcodes();
 	my %mapping_table;		# Will actually be a HoH with structure $HoH{pool name}{sample name}=literal barcode sequence
 	my %sample_names_seen;	# To help verify there is no sample name duplication within a mapping file (which could wreak havoc on later steps of the pipeline)
 	my $bclength;
 	open (MAPPING, $file_path) || die "Can't open mapping file at $mapping_file_path!\n";
 	while (<MAPPING>) {
		chomp;	# Consider replacing this with sequential checks for every EOL format (e.g., check for \r\n, then \r, then \n)
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
 			die "ERROR: There is more than one instance of sample name $line[1] in your mapping file. ",
 				"All sample names within your mapping file must be unique to avoid problems in downstream data summarization.\n";
 		}
 	}
 	close MAPPING;
 	return \%mapping_table;
}

# Pass this sub a ref to sample list, the ref to the mapping hash, and the .bc file path to write
sub make_bc_file {
	my ($mapping_hash, $pool, $sample_list, $bc_file_path) = @_;
	my @barcodes;
	my $bclength;
	foreach(@$sample_list) {
		push(@barcodes, $$mapping_hash{$pool}{$_});
	}
	test_bc_list_compatibility(\@barcodes);
	$bclength = get_bc_length(@barcodes);
	open(BCFILE, ">$bc_file_path") || die "ERROR: Can't open .bc file $bc_file_path!\n";
	foreach(@$sample_list) {
		print BCFILE "$$mapping_hash{$pool}{$_}\t$_\n";
	}
	close BCFILE;
	return $bclength;
}

# Pass it a reference to a list of barcode sequences
# Kill the program if there are any problems and report what the problem is
sub test_bc_list_compatibility {
	my $list_ref = shift;
	my @barcodes = @$list_ref;
	my %barcodes_seen;
	my $bclength;
	foreach(@barcodes) {
		# Check if this looks like an appropriate barcode (i.e. is only A,C,T,G)
		if (/[^ATCGN]/) {
			die "\nERROR: One of your barcodes ($_) has invalid characters.\n";
		}
		# Make sure this barcode is the same length as previous barcodes
		if (!$bclength) {
			$bclength = length($_); 
		}
		elsif ($bclength != length($_)) {
			die "ERROR: Your barcodes are not of uniform length.  The COPRO-Seq workflow does not",
				" currently support analysis of pool subsets containing barcodes of mixed length.\n";
		}
		if ($barcodes_seen{$_}) {
			die "ERROR: One of your pools has more than one sample with the same barcode designation.\n";
		}
		$barcodes_seen{$_} = 1;
	}
	return 1;
}

sub make_get_data_file {
	my ($spreadsheet_hash, $filepath) = @_;
	open (GETDATA, ">$filepath") || die "ERROR: Can't open $filepath!\n";
	my %paths_created;
	# For each row in the filtered spreadsheet hash...
	for my $p (@$spreadsheet_hash) {
		# If a machine, run and lane are defined...
		if ($p->{machine} && $p->{run} && $p->{lane}) {
			if ($p->{path}) {
				my $source_path = $p->{path};
				my $link_path = $p->{machine}.'_'.$p->{run}.'_'.$p->{lane}.'.seq';
				# If we haven't seen this file before (don't want to create the
				#	same link over and over)...
				if (!$paths_created{$link_path}) {
					if (-e $link_path) {
						print GETDATA "echo 'Skipping creation of symbolic ",
							"link to $source_path (a link to this file ",
							"already exists)...'\n";
					}
					else {
						print GETDATA "echo Creating symbolic link to ",
							"$source_path...\n";
						print GETDATA "ln -s $source_path $link_path\n";
						print GETDATA "echo Done.\n";
					}
					$paths_created{$link_path}=1;
				}
			}
			else {
				die "\nERROR: You have not specified a full file path for pool $p->{pool}!\n\n";
			}
		}
		else {
			die "\nERROR: the machine and/or run and/or lane # for pool $p->{pool} are not ",
			"specified. Check your spreadsheet!\n\n";
		}
	}
	close GETDATA;
	return 1;
}

sub squash_genomes {
	my $dir = shift;
	print "Starting genome squashing...\n\n";
	if (! -d $dir)	{	`mkdir $dir`;	}
	else {	print "Note: Directory '$dir' already exists.\n";	}
	`find ./$genomesdir -maxdepth 1 -type f -print | xargs $squash_exec_path $squashedgenomesdir`;
}

sub download_genomes {
	my ($species_list_ref, $source) = @_;
	if (! -d $genomesdir)	{	`mkdir $genomesdir`;	}
	else {	print "Warning: Directory 'genomes' already exists.\n";	}

	if ($source eq 'ncbi') {
		print "\nStarting genome downloads from NCBI...\n\n";
		download_files_from_genbank($species_list_ref, $genomesdir);
	}
	elsif ($source eq 'microbialomics') {
		print "\nStarting genome downloads from microbialomics...\n\n";
		download_files_from_microbialomics($species_list_ref, $genomesdir);
	}
}

sub download_files_from_genbank {
	my ($species, $outputdir) = @_;
	my %species_for_acc = genomecodes::declare_speciesname_for_genbankacc();
	for my $s (@$species) {
		my $filefriendlyname = $species_for_acc{$s};
		$filefriendlyname =~ s/\s+/_/;
		$filefriendlyname =~ s/\.//;
		my $filepath = "$outputdir\/$filefriendlyname";
		open (OUT, ">$filepath") || die "Can't open $filepath!\n";
		my $gb = new Bio::DB::GenBank (-retrievaltype => 'tempfile', -format => 'Fasta');
		# -retrievaltype => 'tempfile' avoids keeping all of the file in memory
		# -format => 'Fasta' avoids downloading additional information we don't need (e.g. features)
		my $seq = $gb->get_Seq_by_acc($s);
		print OUT ">$filefriendlyname\n";
		print OUT $seq->seq();
		close OUT;
		print "Download of $species_for_acc{$s} genome \($s\) from NCBI complete\n";
	}
}

sub get_species_list {
	my $spreadsheet_info_ref = shift;
	my %species_from_hash;	# Used to maintain a non-redundant list of species
	my %genomecodes;
	if ($ncbi) {
		%genomecodes = genomecodes::declare_genbankacc_for_genomeabbrev(); 
	}
	else {	
		%genomecodes = genomecodes::declare_inhouseacc_for_genomeabbrev();
		# Key = abbreviation (e.g. BACCAC), value = in-house accession #
	}
	foreach my $p (@$spreadsheet_info_ref) {
		if ($p->{genomes}) {
			my $allspecies = $p->{genomes};
			my @splitspecies = split(/,/,$allspecies);
			foreach (@splitspecies) {
				s/^\s+//;	# Eliminate any leading whitespace
				s/\s+$//;	# Eliminate any trailing whitespace
				if ($genomecodes{$_}) {	# If a species abbreviation is defined
					$species_from_hash{$genomecodes{$_}} = 1;
				}
				# Use accession # as key if no species name is defined
				else {
					$species_from_hash{$_} = 1;
				}
			}
		}
	}
	my @species_list = keys %species_from_hash;
	return \@species_list;
}

sub download_files_from_microbialomics {
	my ($species, $outputdir) = @_;		# Reference to an array of species ids, location to which to write genome files
	my %species_for_acc = genomecodes::declare_speciesname_for_inhouseacc();
	# Note that the lookup for species names will currently cause problems if there is no species name defined for an 
	#	accession # in genomecodes.pm; need to check if a species name is defined, and if so use it, otherwise just use
	#	the accession # instead (I don't think this should cause any problems)
	# Ultimately I need an efficient way of building a genomecodes.pm file from J's database so that it doesn't need to be
	#	manually updated as changes are made and new strains/species are introduced
	print "Will use database $db_name\.\.\.\n\n";
	for my $s (@$species) {
		my $filefriendlyname;
		if (defined $species_for_acc{$s}) {
			$filefriendlyname = $species_for_acc{$s};
		}
		else {
			print "Note: There is no species name defined for accession $s in genomecodes.pm. Will use $s as the species name instead.\n"; 
			$filefriendlyname = $s;
		}
		$filefriendlyname =~ s/\s+/_/;
		$filefriendlyname =~ s/\.//;
		# Check if genome has already been downloaded...
		if (-e "$outputdir\/$filefriendlyname") {
			if (defined $species_for_acc{$s}) {
				print "Note: It looks like the genome for $species_for_acc{$s} has already been downloaded.\n"; 
			}
			else {
				print "Note: It looks like the genome for $s has already been downloaded.\n";
			}
		}
		else {
			my ($genome_id, $genome_seq) = get_genome_data($s);
			my $filepath = "$outputdir\/$filefriendlyname";
			open(OUT, ">$filepath") || die "Can't open $filepath!\n";
			print OUT ">$filefriendlyname\n$genome_seq";
			close OUT;
			if (defined $species_for_acc{$s}) {
				print "Download of $species_for_acc{$s} genome \($s\) from microbialomics complete\n";
			}
			else {
				print "Download of $s genome from microbialomics complete\n";
			}
		}
	}
}

sub get_genome_data {
	my $constant_id=shift;
	my $query  = "SELECT genome_id, genome_sequence FROM genome g WHERE g.genome_constant_id=?";
	my $sth = $dbh->prepare($query);
	$sth->execute($constant_id);
	if ($sth->rows == 0) {
		die "\nERROR: There is no genome with genome_constant_id $constant_id",
			" in $db_name. Check your accession # or genome abbreviation to be",
			" sure it is correct.\n\n";
	}
	while (my @res = $sth->fetchrow_array()) {
		return @res;
	}
	return;
}

# table_to_hash(table file)
sub table_to_hash {
	my $in = shift;
	open (IN, $in) or die "Can't open project information file $in: $!\n";
	my $line = <IN>;	# First line is tab-separated headers
	chomp $line;
	my @headers = split /\t/, $line;
	for (@headers) {	
		$_ = lc($_);	# Convert all headers to lowercase
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

# filter_by_group(ref to all spreadsheet hashes, ref to hash of acceptable group codes)
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

sub download_table {
	my ($url, $outfile) = @_;
	print STDERR "\n" . fancy_title("Downloading contents of $url to $outfile");
	sleep(1 + 2*rand());	# Don't piss off google (ensures if you are batch-processing multiple 
							# different analyses you don't tick off Google's servers for
							# making a lot of requests for info in a short time (which would look like a DoS attack)
	my $download_info = `wget '$url' -O $outfile -t 3`;
	# URL passed to wget must be surrounded by single quotes or other parameters will be missed
	# wget usage is: wget <where to look> -O <where to send the info> -t <how many times to try>
}

sub fancy_title {
	my $title = shift;
	my $head_bar;
	for (0 .. length($title)+11) { $head_bar .= '-'; }
	return "$head_bar\n|     $title     |\n$head_bar\n";
}

sub usage {
	my $error = shift;
	print "$error\n" if $error;
	print 	"Usage: perl batch_coproseq.pl -g <analysis group(s)> -m <mapping file>\n",
			"\t-e Number of errors/mismatches to allow in alignment (0-2)\n",
			"\t-g Code(s) specifying the analysis groups (defined in Google Spreadsheet) you want to include in this analysis\n",
			"\t-i IGS table file (optional)\n",
			"\t-k Google spreadsheet key for specifying spreadsheets other than the default\n",
			"\t-l Length of read (after trimming barcode) to use in alignment\n",
			"\t-m Mapping file specifying the barcode associated with each sample in your analysis\n",
			"\t-ncbi Specifies that genome files should be downloaded from the NCBI web server, rather than microbialomics\n",
			"\t-p Google spreadsheet page/ID for specifying sheets other than the first in a multi-page spreadsheet; 1st page by default (i.e. -i 0)\n",
			"\n";
	exit(1);
}


# Takes any run number from the Google spreadsheet of 4 digits or less, 
# excluding any decimals specifying run version numbers, and converts it to the
# format used in the sequencing results file system (e.g. if run value is 55.1, 
# this number is converted to 0055.1)
sub run_to_directory {
	my $text=shift;
	my $dir_name_len = 4;	# Update if # of runs ever exceeds 9999
	my $dir_base = "/srv/seq/solexa/results/";
	my @pieces = split '\.', $text;				# Separate run # from run version # (e.g. 55 from .1)
	my $dir = shift @pieces;
	if (length($dir) <= 4) {
		# Append leading zeros until run # is 4 digits
		while (length($dir) < $dir_name_len) {
			$dir = "0" . $dir;
		}
		unshift @pieces, $dir;
		my $final_dir = join ".", @pieces;
		return $dir_base . $final_dir . "/";
	}
	else {
		# Denotes there was an error
		return "0";
	}
}

sub check_options {
	if (!$mapping_file_path) { 
		usage("\nERROR: You must pass the required parameter -m\n"); 
	}
	if (($mismatches_allowed < 0) || ($mismatches_allowed > 2)) {
		usage("\nERROR: The number of mismatches allowed is restricted to values from 0-2\n");
	}
	if (!$group) {
		usage("\nERROR: You must pass the required parameter -g\n");
	}
	if ($readsize > 32) {
		usage("\nERROR: The CoPro-seq workflow cannot currently align sequences greater than 32bp.  Please specify another read length.\n");
	}
	if ($GEO) {
		if (! -d $GEOdir)	{	`mkdir $GEOdir`;	}
	}
}

sub write_GEO_sh {
	my ($outfile, $commands) = @_;
	open(SH, $outfile) || die "Could not create output file at $outfile!\n";
	print SH $commands;
	close SH;
}

# make_modified_bc_file($p->{run}, $p->{lane}, \%$mapping_hash{$p->{pool}}, \@filteredsamples, $filepath);

# Need:
# The mapping hash for the current pool (to know which barcodes go with which samples)
#	Won't need to be the whole mapping hash...just the (%h{sample} = barcode sequence) sub-hash that's
#		normally part of %HoH{pool}{sample} = BC sequence
# The list of (filtered) samples from that pool to include
# The run and lane numbers (for creating the filenames of the barcode-split SCARFs in the .bc_new file)
# The .bc_new filepath

# Example output:
# AAAT	run100_lane1_AAAT.scarf
# TGGT	run100_lane1_TGGT.scarf

sub make_modified_bc_file {
	my ($run, $lane, $mapping_hash_specific_to_pool, $samples_array, $filepath) = @_;
	open(BC, ">$filepath") || die "Can't open modified .bc file $filepath!\n";
	my ($bcseq, $bcsplit_scarf_file);
	foreach(@$samples_array) { 
		$bcseq = $$mapping_hash_specific_to_pool{$_};
		$bcsplit_scarf_file = "run".$run."_lane".$lane."_".$bcseq."_".$_."\.scarf";
		print BC $bcseq."\t".$bcsplit_scarf_file."\n";
	}
	close BC;
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

sub prepare_references {
	my $species_list_ref = shift;
	my @species_names;
	if ($ncbi) {
		download_genomes($species_list_ref, 'ncbi');
		my %species_for_acc = genomecodes::declare_speciesname_for_genbankacc();
		foreach(@$species_list_ref) {
			push(@species_names, $species_for_acc{$_});
		}
	}
	else { 
		download_genomes($species_list_ref, 'microbialomics');
		my %species_for_acc = genomecodes::declare_speciesname_for_inhouseacc();
		# Key = accession #, value = species name
		foreach(@$species_list_ref) {	# $species_list_ref points to array of accession #'s, 
										# some of which may not be defined in species_for_acc
			if (defined $species_for_acc{$_}) {
				push(@species_names, $species_for_acc{$_});
			}
			else {
				push(@species_names, $_);	# Use accession # as species name if species name not defined in genomecodes.pm
			}
		}
	}
	print "\n";
	return @species_names;
}

sub write_GEO_readme {
	my $file = shift;
	open(README, $file) || die "Can't create README output file at $file!\n";
	print	README "Legend for .hitratios files\n",
					"===========================\n",
					"First line: <Sample name>\n",
					"Column 1: Reference genome name\n",
					"Column 2: \"Unique raw counts\", i.e., the number of sequencing reads from the\n",
					"\tsample's total dataset that could be aligned to a single position in the reference\n",
					"\tgenome in column 1, but not to any other position in any of the genomes included in\n",
					"\tthe .hitratios file\n",
					"Column 3: \"Normalized counts\", i.e., the number resulting from adjusting the\n",
					"\tunique raw counts value by its respective genome's length and uniqueness (for a k-mer\n",
					"\tsize equivalent to the read length) relative to all other reference genomes in\n",
					"\tthe .hitratios file\n",
					"Column 4: \"Normalized relative proportion/percentage\", i.e., the total\n",
					"\tnormalized counts attributed to a particular reference genome divided by the sum\n",
					"\tof all normalized counts attributable to all reference genomes in the .hitratios file\n";
	close README;
}

# Input: comma-delimited string specifying analysis groups (e.g., 'A,B')
# Return: reference to a group lookup hash (key = group, value = 1)
sub parse_groups {
	my $group_string = shift;
	my %groups;
	my @elements = split(/,/, $group_string);
	foreach (@elements) {
		# Eliminate leading and trailing whitespace that might be present
		# Otherwise, specifying -g 'A, B' would result in groups "A" and " B"
		s/^\s+//;
		s/\s+$//;
		$groups{$_} = 1;
	}
	return \%groups;
}

sub check_for_absent_groups {
	my ($spreadsheet_hash_ref, $groups_hash_ref) = @_;
	my (%groups_in_spreadsheet, %groups_missing_in_spreadsheet);
	# Create lookup hash of all groups found in analysis spreadsheet
	for my $p (@$project_data_hash) { 
		if (defined $p->{group}) {
			$groups_in_spreadsheet{$p->{group}} = 1;
		}
	}
	# Check to see if any of the groups specified by the user are not
	# found in the spreadsheet (i.e., are not present in the lookup hash)
	foreach (keys %$groups_hash_ref) {
		if (! defined $groups_in_spreadsheet{$_}) {
			$groups_missing_in_spreadsheet{$_} = 1;
		}
	}
	if (keys %groups_missing_in_spreadsheet > 0) {
		my $missing;
		foreach (keys %groups_missing_in_spreadsheet) {
			$missing .= $_."\n";
		}
		die	"ERROR: ".scalar(keys %groups_missing_in_spreadsheet)." group ",
			"names passed with '-g' could not be found in your analysis ",
			"spreadsheet:\n\n$missing\n";
	}
}

# APPENDIX A
#===============================================================================
# $filtered_data_hash is a reference to an array:
#	\@array;
# This referenced array contains references to hashes (one hash per sample)
#	@array = (\%sample_1, \%sample_2, \%sample_3, ..)
# Each hash contains info for one sample (key = column name; value = value)
# 	$hash{run} = 100
# 	$hash{lane} = 1
#	$hash{genomes} = "BACCAC,BACOVA,BACTHE"

# To cycle through this final data structure, use a loop:
#	for my $p (@$project_data_hash) { $p->{run} }
# Example above would access all values for the "run" column
