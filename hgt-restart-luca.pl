#! /usr/bin/perl -w
# Strict Pragmas

=head1 NAME

hgt-restart-luca<.pl>

=head1 USAGE

 hgt-restart-luca.pl [options -v,-d,-h] <ARGS>
 
 Example usage: 
 hgt-restart-luca.pl --check n -itr 5000 -p nocores -i treefile -o output

=head1 SYNOPSIS

A script for analysising the global extent of horizontal gene trasnfer by studying the phyletic pattern of domain architectures across a given tree. 

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------

use lib "/home/sardar/bin/perl-libs-custom/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__

use DBI;
use Supfam::Utils;
use Supfam::hgt;
use Supfam::SQLFunc;
use Supfam::TreeFuncs;

use File::Temp;
use Time::HiRes;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

use Parallel::ForkManager;

use List::Util qw(sum);#Used in generating summary statistics

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $TotalTic = Time::HiRes::time; #USed in timing the total runtime of script

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'HGTResults';
my $TreeFile;
my $maxProcs = 0;
my $test = 'n';
my $FalseNegativeRate = 0.00;
my $completes = 'n'; #A flag to specify whether or not to include architectures including _gap_ assignmenets
my $Iterations = 500;
my $model = 'poisson';
my $store = 0;
my $check = 'y'; #Perform a sanity check on the tree? This should be 'y' unless under extreme circumstances

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$OutputFilename,
           "tree|i:s" => \$TreeFile,
           "processor_cores|p:i" => \$maxProcs,
           "self_test|st:s" => \$test,
           "completes|comp:s" => \$completes,
           "no_iternations|itr:i" => \$Iterations,
           "fals_negative_rate|fnr:f" => \$FalseNegativeRate,
           "model|m:s" => \$model,
           "store|s:i" => \$store,
           "check|c:s" => \$check,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


die "Inappropriate model chosen\n" unless ($model eq 'Julian' || $model eq 'poisson' || $model eq 'corrpoisson');

#---------------------------------------

my $dbh = dbConnect();

`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';
#Path to a piece of RAM that we can write to. This could be on hard disk, but on *nix systems /dev/shm is a shared RAM disk folder. We want a temporary folder in her that will be cleaned up on Prgram exit

# Make a temporary directory for output data 
my $RAWPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $HTMLPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $DELSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $SIMULATIONSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $SELFTERSTPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $BIASPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my ($root,$TreeCacheHash,$tree) = BuildTreeCacheHash($TreeFile); # A massive limitation of this script is the use of BioPerl TreeIO which is SLOOOOW. This is a lookup hash to speed things up.

open RUNINFO, ">HGT_info.$OutputFilename";

print STDERR "No of iterations per run is: $Iterations\n";
print STDERR "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";
print STDERR "False Negative Rate:".$FalseNegativeRate."\n";
print STDERR "Model used: $model\n";
print STDERR "Treefile: $TreeFile \n";
print STDERR "Complete Architectures?: $completes\n";
print STDERR "Command line invocation: $0 ";
my $CommandLine = join(' ',@ARGV);
print STDERR $CommandLine."\n";

print RUNINFO "No of iterations per run is: $Iterations\n";
print RUNINFO "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";
print RUNINFO "False Negative Rate:".$FalseNegativeRate."\n";
print RUNINFO "Model used: $model\n";
print RUNINFO "Treefile: $TreeFile \n";
print RUNINFO "Complete Architectures?: $completes\n\n\n";
print RUNINFO "Command line invocation: $0 ";
print RUNINFO $CommandLine."\n\n\n";

close RUNINFO;

open TREEINFO, ">HGT_tree.$OutputFilename";
my $NewickTree = TreeHash2Newick($TreeCacheHash,$root);
print TREEINFO $NewickTree."\n";
close TREEINFO;
# Dump tree into an output file

#--Check Input Tree------------------------

my $NodeNameIDHash = {}; #This will be a hash to map from NodeName to the BioPerl NodeID $NodeNameIDHash->{NodeName}=NodeID
my @TreeGenomes = map{$_->id}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
map{$NodeNameIDHash->{$_->id}=$_}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; #Populate $NodeNameIDHash

my $genquery = join ("' or genome.genome='", @TreeGenomes); $genquery = "(genome.genome='$genquery')";# An ugly way to make the query run

my $sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();

while (my @query = $sth->fetchrow_array() ) {
		
	die "Genome: $query[0] is in the tree but should not be!\n"  unless ($query[0] eq "y" && $query[1] eq '' || $check eq 'n');
} #A brief bit of error checking - make sure that the tree doesn't contain any unwanted genoemes (e.g strains)
			
	#---------Get a list of domain archs present in tree----------------------

my $lensupraquery = join ("' or len_supra.genome='", @TreeGenomes); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run, but perl DBI only alows for a single value to occupy a wildcard

my $DomCombGenomeHash = {};

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
							AND comb_index.comb NOT LIKE '%_gap_%';"); #comb_id =1 is '_gap_'
}else{
	die "Inappropriate flag for whether or not to include architectures containing _gap_";
}

$sth->execute();

while (my ($genomereturned,$combreturned) = $sth->fetchrow_array() ){

	die "$combreturned\n" if($combreturned =~ /_gap_/ && $completes eq 'y'); # sanity check, will likely remove from final version of script

	$DomCombGenomeHash->{$combreturned} = {} unless (exists($DomCombGenomeHash->{$combreturned}));
	$DomCombGenomeHash->{$combreturned}{$genomereturned}++;
}

my $toc = Time::HiRes::time;
print STDERR "Time taken to build the Dom Arch hash:".($toc-$tic)."seconds\n";

my $DomArchs = [];
@$DomArchs = keys(%$DomCombGenomeHash);
#These are all the unique domain architectures

dbDisconnect($dbh) ; 

my $NoOfForks = $maxProcs;RUNINFO
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

#----------------------------------------------------


#Main-loop-------------------------------------------

print STDERR "Total No Of Dom Archs: ".@$DomArchs."\n";

my $pm = new Parallel::ForkManager($maxProcs) if ($maxProcs);# Initialise

foreach my $fork (0 .. $NoOfForks-1){
	
	my $ArchsListRef = $ForkJobsHash->{$fork};
			
	# Forks and returns the pid for the child:
	if ($maxProcs){$pm->start and next};
		
	my $CachedResults = {}; #Allow for caching of distributions after Random Model to speed things up
		
	open HTML, ">$HTMLPATH/$OutputFilename".$$.".html" or die "Can't open file $HTMLPATH/$OutputFilename".$!;
	open OUT, ">$RAWPATH/$OutputFilename".$$.".-RawData.colsv" or die "Can't open file $RAWPATH/$OutputFilename".$!;
	open DELS, ">$DELSPATH/DelRates".$$.".dat" or die "Can't open file $DELSPATH/DelRates".$!;
	open RAWSIM, ">$SIMULATIONSPATH/SimulationData".$$.".dat" or die "Can't open file $SIMULATIONSPATH/SimulationData".$!;		
	open SELFTEST, ">$SELFTERSTPATH/SelfTestData".$$.".dat" or die "Can't open file $SELFTERSTPATH/SelfTestData".$!;
	open BIAS, ">$BIASPATH/BiasEstimate".$$.".dat" or die "Can't open file $BIASPATH/BiasEstimate".$!;
			
	foreach my $DomArch (@$ArchsListRef){
		
	my ($CladeGenomes,$NodesObserved);
			
		my $HashOfGenomesObserved = $DomCombGenomeHash->{$DomArch};
		@$NodesObserved = keys(%$HashOfGenomesObserved);
		my @NodeIDsObserved = map{$NodeNameIDHash->{$_}}@$NodesObserved ; 
		#Get the node IDs as the follwoing function doesn't work with the raw node tags
	
		my $MRCA;
		my $deletion_rate;
		my ($dels, $time) = (0,0);
		
		unless(scalar(@NodeIDsObserved) == 1){
		
		$MRCA = $tree->get_lca(-nodes => \@NodeIDsObserved) ; #Find most Recent Common Ancestor	
		
		 if($model eq 'Julian' || $model eq 'poisson' || $model eq 'corrpoisson'){
		 	
				($dels, $time) = DeletedJulian($MRCA,0,0,$HashOfGenomesObserved,$TreeCacheHash,$root,$DomArch); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree	
				
			}else{
				die "Inappropriate model selected";
			}
			
		$deletion_rate = $dels/$time;
		
		}else{
			$deletion_rate = 0;	
			$MRCA = $NodeIDsObserved[0] ; #Most Recent Common Ancestor
		}
		
		@$CladeGenomes = @{$TreeCacheHash->{$MRCA}{'Clade_Leaves'}}; # Get all leaf genomes in this clade	
		
		print DELS "$DomArch:$deletion_rate:$dels\n" unless ($deletion_rate == 0);
		#print "$DomArch:$deletion_rate\n";
		
		my ($selftest,$distribution,$RawResults,$DeletionsNumberDistribution);
		
		if($deletion_rate > 0){
	
			unless($CachedResults->{"$deletion_rate:@$CladeGenomes"} && $store){
						
				if($model eq 'Julian'){
									
					($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelJulian($MRCA,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);
																					
				}elsif($model eq 'poisson'){
					
					($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelPoisson($MRCA,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);

				}elsif($model eq 'corrpoisson'){
					
					($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelCorrPoisson($MRCA,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);

				}else{
					die "No appropriate model selected";
				}
				
			$CachedResults->{"$deletion_rate:@$CladeGenomes"} = [$selftest,$distribution,$RawResults,$DeletionsNumberDistribution];		
			
			}else{
				($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = @{$CachedResults->{"$deletion_rate:@$CladeGenomes"}};
			}
			
			my $RawSimData = join(',',@$RawResults);
			print RAWSIM @$CladeGenomes.','.@NodeIDsObserved.':'.$DomArch.':'.$RawSimData."\n";
			#Print simulation data out to file so as to allow for testing of convergence
			
		}else{
			
			($selftest,$distribution) = ('NULL',{});
		}
		
		
#-------------- Output
		 
		my $NoGenomesObserved = scalar(@NodeIDsObserved);
		my $CladeSize = scalar(@$CladeGenomes);
				
		unless($selftest eq 'NULL'){ #Essentially, unless the deletion rate is zero
			
		#Gain an estimate of bias in the simulation (shoudl be 0)
		my @SimulatedNumberOfDeletions;
		map{@SimulatedNumberOfDeletions=(@SimulatedNumberOfDeletions,(($_)x($DeletionsNumberDistribution->{$_})))}keys(%$DeletionsNumberDistribution); #Produces a list of the simualted number of deletions
		
		my $DeletionNumberMean = sum(@SimulatedNumberOfDeletions)/scalar(@SimulatedNumberOfDeletions);
		my $BiasEstimate = $DeletionNumberMean-$dels; #This is a measure as to how bias the simulation is - the mean of the number of deletions in a simualtion shoudl be the number of deletions observed
		#
			
			
			my $CumulativeCount= 0;
			my @CumulativeDistribution; #(P(Nm<nr)) nr (number of genomes in reality) is the index, Nm is the random variable
			
			#Create cumulative distibution list
			for my $NumberOfGenomes (0 .. $CladeSize){
		
				$CumulativeCount += $distribution->{$NumberOfGenomes} if(exists($distribution->{$NumberOfGenomes})); # Cumulative count is a sum of the frequency of genome observation up to this point
				$CumulativeDistribution[$NumberOfGenomes]=$CumulativeCount;	
			}
			
			my $ProbLT = ($CumulativeDistribution[$NoGenomesObserved-1])/$Iterations;#Cumlative probability of less genomes in model simulations than in reality
			$selftest +=1 if ($selftest == 0); #Prevents selfcheck from being less than 0.
			my $SelfTestLT = ($CumulativeDistribution[$selftest-1])/$Iterations;#Cumlative probability of less genomes in model simulations than in self check
			
			#Self test is a measure of how reliable the simualtion is and whether we have achieved convergence - one random genome is chosen as a substitute for 'reality'.
			
				        print HTML "<a href=http://http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/maketree.cgi?genomes=";
				        print HTML join(',', @$NodesObserved);
				        print HTML ">$DomArch</a> Score: $ProbLT<BR>\n";
						
						print STDERR $DomArch."\n" unless($DomArch);
						$DomArch = 'NULL' unless($DomArch);
						
						print OUT "$DomArch:$ProbLT\n";
						print SELFTEST "$DomArch:$SelfTestLT\n";
						#The output value 'Score:' is the probability, givn the model, that there are more genomes in the simulation than in reality.
						print BIAS "$DomArch:$BiasEstimate\n";
		}
}

	close HTML;
	close OUT;
	close DELS;
	close RAWSIM;
	close SELFTEST;
	close BIAS;
	
$pm->finish if ($maxProcs); # Terminates the child process
}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";

`cat $HTMLPATH/* > ./$OutputFilename.html`;
`cat $DELSPATH/* > ./.DelRates.dat`;
`cat $RAWPATH/* > ./$OutputFilename-RawData.colsv`;
`cat $SIMULATIONSPATH/* > ./RawSimulationDists$Iterations-Itr$$.dat`;
`cat $SELFTERSTPATH/* > ./SelfTest-RawData.colsv`;
`cat $BIASPATH/* > ./BiasEstimates.colsv`;


`Hist.py -f "./$OutputFilename-RawData.colsv" -o $OutputFilename.png -t "Histogram of Cumulative p-Values" -x "P(Nm < nr)" -y "Frequency"	` ;
`Hist.py -f "./SelfTest-RawData.colsv" -o SelfTest.png -t "Histogram of Self-Test Cumulative p-Values" -x "P(Nm < nm)" -y "Frequency"	` ;
`Hist.py -f "./.DelRates.dat" -o ParDelRates.png -t "Histogram of Non-zero DeletionRates" -x "Deletion Rate" -y "Frequency"	-l Deletions`;
`Hist.py -f ./BiasEstimates.colsv -o Bias.png --xlab 'Simualtion' --ylab 'Bias' --title 'Estimate of per simulation bias in deletion probability distribution'`;

# Plot a couple of histograms for easy inspection of the data


open PLOT, ">./.ParPlotCommands.txt" or die $!;
print PLOT "Hist.py -f ./$OutputFilename-RawData.colsv -o $OutputFilename.png -t 'Histogram of Scores' -x 'Score' -y 'Frequency'\n\n\n";
print PLOT "Hist.py -f './.DelRates.dat' -o 'ParDelRates.png' -t 'Histogram of Non-zero DeletionRates' -x 'Deletion Rate' -y 'Frequency'	-l 'Deletions'\n\n\n";
print PLOT "Hist.py -f './SelfTest-RawData.colsv' -o SelfTest.png -t 'Histogram of Self-Test Cumulative p-Values' -x 'P(Nm < nm')' -y 'Frequency'\n\n\n";
print PLOT `Hist.py -f ./BiasEstimates.colsv -o Bias.png --xlab 'Simualtion' --ylab 'Bias' --title 'Estimate of per simulation bias in deletion probability distribution'`;
close PLOT;

my $TotalToc = Time::HiRes::time;
my $TotalTimeTaken = ($TotalToc - $TotalTic);
my $TotalTimeTakenHours = $TotalTimeTaken/(60*60);

open RUNINFOTIME, ">>HGT_info.$OutputFilename";
print RUNINFOTIME $TotalTimeTaken." seconds\n";
print RUNINFOTIME $TotalTimeTakenHours." hours\n";
close RUNINFOTIME;

#-------
__END__


