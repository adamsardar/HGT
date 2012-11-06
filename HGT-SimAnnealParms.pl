#!/usr/bin/env perl

=head1 NAME

HGT-parallel<.pl>

=head1 USAGE

HGT-parallel.pl

Example usage: 

HGT-parallel.pl --check n -itr 5000 -p 10 -i ./tree.nwk -o output --singlesim/--simulations

=head1 SYNOPSIS

A script for analysising the global extent of horizontal gene trasnfer by studying the phyletic pattern of domain architectures across a given tree.

There are a variety of options and modes to specify:

	-h
	
	Shows this help document
	
	--check (y/n) DEFAULT: n
	
	Perform internal checks as to whether a genome is 'include = y' within the SUPERFAMILY database
	
	-i  --tree / -ha --hash
	
	input newick tree file, which will then be quieried against the SUPERFAMILY databse. As an alternative, you may provide a data dump using -ha produced using -di or --singlesim
	
	--simualtions (options of specifying a number of processor cores to use using -p, which is strongly advised).
	
	This will perform full simulations under the model specified using -m (DEFAULT: poisson) and report the placement of 
	
		Additional paramters: -fnr --fals_negative_rate (Flase negative rate)
								 -p --processor_cores (Number of threads to create to perform simualtions)
								 -m --model (model to use in performing simualtions. Choose from 'Julian', 'poisson' or 'corrpoisson')
								 -dm --delmodel (model to use in assigning dleetion rates. Choose from 'Julian', 'Uniform', 'Geometric' and 'Power')
								 -s --store (minor speedup optimisation, it caches simualtion results so that they can be reused if identical paramters occur more than once. DEFAULT: FALSE)
	
	-o --output
	
	Output filename stub
	
	--singlesim
	
	Using the model conditions specified, this will perform a single simulations per domain architecture and then spit out a datastructure of the same form as -di.
	You can then load these back into the program using -ha
	
	-di --dump_input
	
	Dump intput. This will dump out the perl datastructure for use in inspecting or reloading data. Reload using -ha
	

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output and for dumping out datastructures
B<Math::Random> Used in Monte Carlo simulations steps
B<Parallel::ForkManager> Used for multi-threaded simualtions
B<DBI> Needed for database querying
B<Time::HiRes> Used to produce time estimates of event durations
B<List::Util> Sum function needed in producing summary statistics
B<File::Temp> Parallel implementation writes to temporary files in ram. This module ensures that they are cleaned up
B<Supfam::*> Supfam Toolkit provides a variety of custom tools needed. From simualtions to parsing newick trees.
=cut


use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Carp;

use DBI;

use Supfam::Utils;
use Supfam::hgt;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;

use File::Temp;
use Time::HiRes;

use Parallel::ForkManager;
use Math::Random qw(random_uniform_integer random_uniform random_exponential);
use List::Util qw(sum);#Used in generating summary statistics

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $TotalTic = Time::HiRes::time; #USed in timing the total runtime of script

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'HGTResults';
my $TreeFile;
my $SpeciesAlignmentFile;
my $maxProcs = 0;
my $test = 'n';
my $FalseNegativeRate = 0.00;
my $completes = 'n'; #A flag to specify whether or not to include architectures including _gap_ assignmenets
my $Iterations = 500;
my $model = 'poisson';
my $store = 0;
my $check = 'y'; #Perform a sanity check on the tree? This should be 'y' unless under extreme circumstances
my $dumpinput;
my $fullsims; #flag for performing full posterior quantile simulations
my $singlesim;
my $HGTpercentage = 0;
my $delmodel = 'Julian';

my $CommandOps = join("  ",@ARGV);

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$OutputFilename,
           "tree|i:s" => \$TreeFile,
           "hash|ha:s" => \$SpeciesAlignmentFile,
           "processor_cores|p:i" => \$maxProcs,
           "completes|comp:s" => \$completes,
           "no_iternations|itr:i" => \$Iterations,
           "fals_negative_rate|fnr:f" => \$FalseNegativeRate,
           "model|m:s" => \$model,
           "store|s!" => \$store,
           "check|c:s" => \$check,
           "dump_input|di!" => \$dumpinput,
           "simulations!" => \$fullsims,
           "singlesim!" => \$singlesim,
           "delmodel|dm:s" => \$delmodel,
           "HGTpercentage|ht:i" => \$HGTpercentage,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

die "Inappropriate model chosen; models avaialvle are Julian, poisson, corrpoisson, negbin and corrnegbin\n" unless ($model eq 'Julian' || $model eq 'poisson' || $model eq 'corrpoisson' || $model eq 'negbin' || $model eq 'corrnegbin');
#---------------------------------------

`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';
#Path to a piece of RAM that we can write to. This could be on hard disk, but on *nix systems /dev/shm is a shared RAM disk folder. We want a temporary folder in her that will be cleaned up on Prgram exit

# Make a temporary directory for output data 
my $RAWPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $HTMLPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $DELSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $SIMULATIONSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $SELFTERSTPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Produce a tree hash, either from SQL or a provided treefile
my ($root,$TreeCacheHash);
my $DomCombGenomeHash = {};

open RUNINFO, ">HGT_info.$OutputFilename";

if($fullsims){
	print STDERR "Runing a full set of simulations per trait and calculating the observed value placement ...\n";
	print STDERR "No of iterations per run is: $Iterations\n" if($fullsims);
	print STDERR "False Negative Rate:".$FalseNegativeRate."\n" if($fullsims);
	print STDERR "Simualtion Model used: $model\n";
	print STDERR "Cores used: $maxProcs\n" if ($maxProcs > 0);
}
if($singlesim){
	print STDERR "Runing a single simulation per trait and outputting the pseudo-observations ...\n";
	print STDERR "Deletion Model used: $delmodel";
	print STDERR " - HGT percentage: $HGTpercentage %" if($delmodel eq 'Julian');
	print STDERR "\n";
}
print STDERR "Command line invocation: $0 $CommandOps\n";
print STDERR "Complete Architectures?: $completes\n\n\n";


if($fullsims){
	print RUNINFO "Runing a full set of simulations per trait and calculating the observed value placement ...\n";
	print RUNINFO "No of iterations per run is: $Iterations\n" if($fullsims);
	print RUNINFO "False Negative Rate:".$FalseNegativeRate."\n" if($fullsims);
	print RUNINFO "Simualtion Model used: $model\n";
	print RUNINFO "Cores used: $maxProcs\n" if ($maxProcs > 0);
}
if($singlesim){
	print RUNINFO "Runing a single simulation per trait and outputting the pseudo-observations ...\n";
	print RUNINFO "Deletion Model used: $delmodel";
	print RUNINFO " - HGT percentage: $HGTpercentage %" if($delmodel eq 'Julian');
	print RUNINFO "\n";
}
print RUNINFO "Command line invocation: $0 $CommandOps\n";
print RUNINFO "Complete Architectures?: $completes\n\n\n";


if($TreeFile){
	
	open TREE, "<$TreeFile" or die $!.$?;
	my $TreeString = <TREE>;
	close TREE;
	
	($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

	my $dbh = dbConnect();
	my $sth;
	
	my @TreeGenomes = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
	my @TreeGenomesNodeIDs = @{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
	
	#---------Get a list of domain archs present in tree----------------------
	my $lensupraquery = join ("' or len_supra.genome='", @TreeGenomes); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run, but perl DBI only allows for a single value to occupy a wildcard

	my $tic = Time::HiRes::time;
	
	if($completes eq 'n'){
		
		$sth = $dbh->prepare("SELECT DISTINCT len_supra.genome,comb_index.comb 
							FROM len_supra JOIN comb_index ON len_supra.supra_id = comb_index.id 
							WHERE len_supra.ascomb_prot_number > 0 
							AND $lensupraquery 
							AND comb_index.id != 1;"); #comb_id =1 is '_gap_'
			
	}elsif($completes eq 'y'){
	
		$sth = $dbh->prepare("SELECT DISTINCT len_supra.genome,comb_index.comb 
								FROM len_supra JOIN comb_index ON len_supra.supra_id = comb_index.id 
								WHERE len_supra.ascomb_prot_number > 0 
								AND $lensupraquery AND comb_index.id != 1 
								AND comb_index.comb NOT LIKE '%_gap_%';"); 
								#select only architectures that are fully assigned (don't contain _gap_) comb_id =1 is '_gap_'
	}elsif($completes eq 'nc'){
	
		$sth = $dbh->prepare("SELECT DISTINCT len_supra.genome,comb_index.comb 
								FROM len_supra JOIN comb_index ON len_supra.supra_id = comb_index.id 
								WHERE len_supra.ascomb_prot_number > 0 
								AND $lensupraquery AND comb_index.id != 1 
								AND comb_index.comb NOT LIKE '_gap_%' 
								AND comb_index.comb NOT LIKE '%_gap_';");
								#select only architectures that are fully assigned (don't contain _gap_) comb_id =1 is '_gap_'
	}else{
		
		die "Inappropriate flag for whether or not to include architectures containing _gap_";
	}
	
	$sth->execute();
	
	while (my ($genomereturned,$combreturned) = $sth->fetchrow_array() ){
	
		die "$combreturned\n" if($combreturned =~ m/_gap_/ && $completes eq 'y'); # sanity check, will likely remove from final version of script
	
		$DomCombGenomeHash->{$combreturned} = {} unless (exists($DomCombGenomeHash->{$combreturned}));
		$DomCombGenomeHash->{$combreturned}{$genomereturned}++;
	}
	
	my $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Dom Arch hash:".($toc-$tic)."seconds\n";
	
	dbDisconnect($dbh);

}elsif($SpeciesAlignmentFile){
	
	print STDERR "Using species trait file input rather than read in values from SUPERFAMILY.\n";
	
	my $SpecieData = EasyUnDump($SpeciesAlignmentFile);
	($root,$TreeCacheHash,$DomCombGenomeHash) = @$SpecieData;
	
}else{
	
	die "no tree file provided as tree to calculate HGT upon\n";	
}

print STDERR "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";
print RUNINFO "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";

close RUNINFO;

open TREEINFO, ">HGT_tree.$OutputFilename";
my $NewickTree = ExtractNewickSubtree($TreeCacheHash, $root,1,0);
print TREEINFO $NewickTree."\n";
close TREEINFO;
# Dump tree into an output file


if($dumpinput){
	
	print STDERR "Creating a dump of input data\n";
	EasyDump('input_spcies_data.dat',[$root,$TreeCacheHash,$DomCombGenomeHash]);
}

#--Check Input Tree------------------------
my $dbh = dbConnect();
				
my @TreeGenomes = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
my @TreeGenomesNodeIDs = @{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree

my $genquery = join ("' or genome.genome='", @TreeGenomes); $genquery = "(genome.genome='$genquery')";# An ugly way to make the query run

my $sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();

while (my @query = $sth->fetchrow_array() ) {
		
	die "Genome: $query[0] is in the tree but should not be!\n"  unless ($query[0] eq "y" && $query[1] eq '' || $check eq 'n');
} #A brief bit of error checking - make sure that the tree doesn't contain any unwanted genoemes (e.g strains)

dbDisconnect($dbh);
#--------------------------
#Consider moving the above code segement into its own subroutine in Supfam::TreeFuncsNonBP


my $DomArchs = [];
@$DomArchs = keys(%$DomCombGenomeHash);
#These are all the unique domain architectures

print STDERR "Total No Of Dom Archs: ".@$DomArchs."\n";

#Single simultion for each domain architecture
croak "HGT percentage rate inaccurately set. Must be between 0 and 100%\n" unless($HGTpercentage >= 0 && $HGTpercentage <=100);

my $NoOfForks = $maxProcs;
$NoOfForks = 1 unless($maxProcs);

my $remainder = scalar(@$DomArchs)%$NoOfForks;
my $binsize = (scalar(@$DomArchs)-$remainder)/$NoOfForks;

my $ForkJobsHash = {};

for my $i (0 .. $NoOfForks-1){
	
	my @ForkJobList = @{$DomArchs}[$i*$binsize .. ($i+1)*$binsize-1];
	$ForkJobsHash->{$i}=\@ForkJobList;
}
#Create lists of jobs to be done by the relative forks

push(@{$ForkJobsHash->{0}},@{$DomArchs}[($NoOfForks)*$binsize .. ($NoOfForks)*$binsize+$remainder-1]) if ($remainder);

print STDERR "No Dom Archs in job batch is approx: ".$binsize."\n";
	
my $pm = new Parallel::ForkManager($maxProcs) if ($maxProcs);# Initialise

my $counter = 0;

foreach my $fork (0 .. $NoOfForks-1){
	
	my $ArchsListRef = $ForkJobsHash->{$fork};
		
	# Forks and returns the pid for the child:
	if ($maxProcs){$pm->start and next};

	open DELS, ">$DELSPATH/DelRates".$$.".dat" or die "Can't open file $DELSPATH/DelRates".$!;
	

	foreach my $DomArch (@$ArchsListRef){
		
		$counter++;
		my $newtic =  Time::HiRes::time;

		my ($CladeGenomes,$NodesObserved);
		
		my $NodeName2NodeID = {};
		map{$NodeName2NodeID->{$TreeCacheHash->{$_}{'node_id'}}= $_ }@TreeGenomesNodeIDs; #Generate a lookup table of leaf_name 2 node_id
		
		my $HashOfGenomesObserved = $DomCombGenomeHash->{$DomArch};
		@$NodesObserved = keys(%$HashOfGenomesObserved);
		
		my @NodeIDsObserved = @{$NodeName2NodeID}{@$NodesObserved};#Hash slice to extract the node ids of the genomes observed
		#Get the node IDs as the follwoing function doesn't work with the raw node tags
		
		
		
		my $NoGenomesObserved = scalar(@NodeIDsObserved);
		
		my $MRCA;
		my $deletion_rate = undef;
		my ($dels, $time);
		my $TotalBranchLength = undef;
		
		my ($lamba_best, $lambda_original,$BestLambda_PostQuant,$Original_PostQuant,$selftest,$diff);	
		
		my $print = 0;
		
		unless(scalar(@$NodesObserved) == 1){
			
			$MRCA = FindMRCA($TreeCacheHash,$root,\@NodeIDsObserved);#($TreeCacheHash,$root,$LeavesArrayRef)
			
			@$CladeGenomes = @{$TreeCacheHash->{$MRCA}{'Clade_Leaves'}}; # Get all leaf genomes in this clade	
			@$CladeGenomes = ($MRCA) if($TreeCacheHash->{$MRCA}{'is_Leaf'});
			
			my $CladeSize = scalar(@$CladeGenomes);

			($dels, $time) = DeletedJulian($MRCA,0,0,$HashOfGenomesObserved,$TreeCacheHash,$root,$DomArch); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree	
			
			next if($dels == 0);
			
			($lamba_best, $lambda_original,$BestLambda_PostQuant,$Original_PostQuant,$selftest) = DeletedSimAnneal($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$root,$DomArch,$CladeSize,100,$HGTpercentage,$model,$NoGenomesObserved);
			
			$diff = $lambda_original - $lamba_best;
			
			$print = 1;
			
			$deletion_rate = $dels/$time;
			
			$TotalBranchLength = $TreeCacheHash->{$MRCA}{'Total_branch_lengths'};
			
			
			
		}else{
			
			$MRCA = $NodeIDsObserved[0] ; #Most Recent Common Ancestor
		}
								
		print DELS "$DomArch\t$deletion_rate\t$dels:$time\t$TotalBranchLength\t$lamba_best\t$lambda_original\t$diff\t$BestLambda_PostQuant\t$Original_PostQuant\t$selftest\n" if($print);
		#print "$DomArch:$deletion_rate\n";
		
		
		if($verbose){
			my $newtoc =  Time::HiRes::time;
			my $times = $newtoc - $newtic;
			print STDERR $counter."'th DA done in ".$times." seconds and print flag was ".$print."\n";
		}
		
		
		last if ($counter > 40)
		
	}

	close DELS;
	
$pm->finish if ($maxProcs); # Terminates the child process


}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";

`cat $DELSPATH/* > ./SimAnnealDelRates.dat`;

my $TotalToc = Time::HiRes::time;
my $TotalTimeTaken = ($TotalToc - $TotalTic);
my $TotalTimeTakenHours = $TotalTimeTaken/(60*60);

open RUNINFOTIME, ">>HGT_info.$OutputFilename";	
print RUNINFOTIME $TotalTimeTaken." seconds\n";
print RUNINFOTIME $TotalTimeTakenHours." hours\n";
close RUNINFOTIME;
	
print STDERR $TotalTimeTaken." seconds\n";
print STDERR $TotalTimeTakenHours." hours\n";

#-------
__END__
