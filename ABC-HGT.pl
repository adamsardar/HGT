#! /usr/bin/perl -w
# Strict Pragmas

=head1 NAME

ABC-HGT<.pl>

=head1 USAGE

 ABC-HGT.pl [options -v,-d,-h] <ARGS>
 
 Example usage: 

ABC-HGT.pl -r ../Eukaryotes/DombCombSupfamTree/poisson/EukaryotesDomArchUnnormalised-RawData.colsv -i ../EukaryoteTreeDomaarc.tree -o Poisson ABCCheck -p 0 -itr 100 -e 1 --check n


=head1 SYNOPSIS

A means by which to test if a maximum likelihood estimate of 

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
use lib "$ENV{HOME}/bin/perl-libs-custom/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
B<Math::Random> Used in Monte Carlo simulations steps
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__

use DBI;

use Supfam::Utils;
use Supfam::hgt;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;

use File::Temp;
use Time::HiRes;

use Parallel::ForkManager;
use Math::Random qw(random_uniform_integer random_gamma);
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
my $epsilon = 1;#Tolerance used in ABC
my $check = 'y'; #Perform a sanity check on the tree? This should be 'y' unless under extreme circumstances
my $OutputDir = 'MLvsABC';
my ($rawdatafile) = glob("*-RawData.colsv*");

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
           "epsilon|e:s" => \$epsilon,
           "store|s:f" => \$epsilon,
           "check|c:s" => \$check,
           "rawdata|r:s" => \$rawdatafile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

die "Inappropriate model chosen\n" unless ($model eq 'Julian' || $model eq 'poisson' || $model eq 'corrpoisson');
#---------------------------------------

my $dbh = dbConnect();

mkdir("./$OutputDir");

`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';

my $SIMULATIONSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;

#----------------------------------------------------------------------------------------------------------------

#Produce a tree hash, either from SQL or a provided treefile
my ($root,$TreeCacheHash);

if($TreeFile){
	
	open TREE, "<$TreeFile" or die $!.$?;
	my $TreeString = <TREE>;
	close TREE;

	($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

}else{
	
	die "no tree file or SQL tree node or tree file provided as tree to calculate HGT upon\n";	
}

#Read in Scores from a previous simulation
my $DomArch2ScoresHash = {};
open SCORES, "<$rawdatafile" or die $?;

while (my $line = <SCORES>){
	
	chomp($line);
	my ($DomArch,$Score) = split(/:/,$line);
	$DomArch2ScoresHash->{$DomArch}=$Score;
}
close SCORES;



open RUNINFO, ">ABCHGT_info.$OutputFilename";

print STDERR "No of iterations per run is: $Iterations\n";
print STDERR "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";
print STDERR "False Negative Rate:".$FalseNegativeRate."\n";
print STDERR "Model used: $model\n";
print STDERR "Cores used: $maxProcs\n";
print STDERR "Treefile: $TreeFile \n";
print STDERR "Complete Architectures?: $completes\n";
print STDERR "Command line invocation: $0 \n";

print RUNINFO "No of iterations per run is: $Iterations\n";
print RUNINFO "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";
print RUNINFO "False Negative Rate:".$FalseNegativeRate."\n";
print RUNINFO "Model used: $model\n";
print RUNINFO "Cores used: $maxProcs\n";
print RUNINFO "Treefile: $TreeFile \n";
print RUNINFO "Complete Architectures?: $completes\n\n\n";
print RUNINFO "Command line invocation: $0 \n";

close RUNINFO;

open TREEINFO, ">ABCHGT_tree.$OutputFilename";
my $NewickTree = ExtractNewickSubtree($TreeCacheHash, $root,1,0);
print TREEINFO $NewickTree."\n";
close TREEINFO;
# Dump tree into an output file

#--Check Input Tree------------------------


my @TreeGenomes = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
my @TreeGenomesNodeIDs = @{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree


my $genquery = join ("' or genome.genome='", @TreeGenomes); $genquery = "(genome.genome='$genquery')";# An ugly way to make the query run

my $sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();

while (my @query = $sth->fetchrow_array() ) {
		
	die "Genome: $query[0] is in the tree but should not be!\n"  unless ($query[0] eq "y" && $query[1] eq '' || $check eq 'n');
} #A brief bit of error checking - make sure that the tree doesn't contain any unwanted genoemes (e.g strains)
			
	#---------Get a list of domain archs present in tree----------------------

my $lensupraquery = join ("' or len_supra.genome='", @TreeGenomes); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run, but perl DBI only allows for a single value to occupy a wildcard

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

	die "$combreturned\n" if($combreturned =~ m/_gap_/ && $completes eq 'y'); # sanity check, will likely remove from final version of script

	$DomCombGenomeHash->{$combreturned} = {} unless (exists($DomCombGenomeHash->{$combreturned}));
	$DomCombGenomeHash->{$combreturned}{$genomereturned}++;
}

my $toc = Time::HiRes::time;
print STDERR "Time taken to build the Dom Arch hash:".($toc-$tic)."seconds\n";

my $DomArchs = [];
@$DomArchs = keys(%$DomCombGenomeHash);
#These are all the unique domain architectures

dbDisconnect($dbh);

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

#----------------------------------------------------


#Main-loop-------------------------------------------

print STDERR "Total No Of Dom Archs: ".@$DomArchs."\n";

my $pm = new Parallel::ForkManager($maxProcs) if ($maxProcs);# Initialise

foreach my $fork (0 .. $NoOfForks-1){
	
	my $ArchsListRef = $ForkJobsHash->{$fork};
		
		
	# Forks and returns the pid for the child:
	if ($maxProcs){$pm->start and next};
		
	my $CachedResults = {}; #Allow for caching of distributions after Random Model to speed things up
	
	open OUT, ">$SIMULATIONSPATH/SimulationData".$$.".dat" or die "Can't open file $SIMULATIONSPATH/SimulationData".$!;	
	
	foreach my $DomArch (@$ArchsListRef){
	
	my ($CladeGenomes,$NodesObserved);
		
		my $NodeName2NodeID = {};
		map{$NodeName2NodeID->{$TreeCacheHash->{$_}{'node_id'}}= $_ }@TreeGenomesNodeIDs; #Generate a lookup table of leaf_name 2 node_id
		
		my $HashOfGenomesObserved = $DomCombGenomeHash->{$DomArch};
		@$NodesObserved = keys(%$HashOfGenomesObserved);
		
		my @NodeIDsObserved = @{$NodeName2NodeID}{@$NodesObserved};#Hash slice to extract the node ids of the genomes observed
		#Get the node IDs as the follwoing function doesn't work with the raw node tags
		
		my $NumberOfGenomes = scalar(@$NodesObserved);
		
		my $MRCA;
		my $deletion_rate;
		my ($dels, $time) = (0,0);
		
		unless(scalar(@$NodesObserved) == 1){
			
			$MRCA = FindMRCA($TreeCacheHash,$root,\@NodeIDsObserved);#($TreeCacheHash,$root,$LeavesArrayRef)

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
		@$CladeGenomes = ($MRCA) if($TreeCacheHash->{$MRCA}{'is_Leaf'});
		
		my $cladesize = scalar(@$CladeGenomes);
		
		my ($selftest,$distribution,$RawResults,$DeletionsNumberDistribution);
				
		my $PosteriorEstimatesOFDeltionRate = [];	#This are the crucial piece of data generated from this ABC process. The posterior estimates of the parameters
				
		if($deletion_rate > 0){
	
			
			my @PriorRates = random_gamma(2*$Iterations,$time,$dels); 
			#Supply a gamma distribution paramterised by alpha=total time deletions occur in and beta = total number deletions
			

			#While loop over prior array
			while (scalar(@$PosteriorEstimatesOFDeltionRate) < $Iterations){
			
				
				my $simdelrate = pop(@PriorRates);
				
				if(scalar(@PriorRates) == 0){
				
					@PriorRates = random_gamma($Iterations,$time,$dels);
				}
				#If the pool of prior values is runnign out, supply a load more into the same array space
				
				#Perfrom simulation according to model as chosen
					if($model eq 'Julian'){
										
						($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelJulian($MRCA,$FalseNegativeRate,1,$simdelrate,$TreeCacheHash);
																						
					}elsif($model eq 'poisson'){
						
						($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelPoisson($MRCA,$FalseNegativeRate,1,$simdelrate,$TreeCacheHash);
	
					}elsif($model eq 'corrpoisson'){
						
						($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelCorrPoisson($MRCA,$FalseNegativeRate,1,$simdelrate,$TreeCacheHash);
	
					}else{
						die "No appropriate model selected";
					}
					
				
				map{push(@$PosteriorEstimatesOFDeltionRate,$simdelrate) if(abs($NumberOfGenomes-$_) < $epsilon)}@$RawResults; #Add simulated deletion rate to the posterior if a simulation performed at that rate created a reasonably close value
				
			}	
			
			open TEMPPOSTERIOR, ">./.$$-TempPosterior.dat" or die $!.$?;
			print TEMPPOSTERIOR join("\n",@$PosteriorEstimatesOFDeltionRate);
			close TEMPPOSTERIOR;
			
			
			my $Score = 'Undefined';
			
			if(exists($DomArch2ScoresHash->{$DomArch})){
				
				$Score = $DomArch2ScoresHash->{$DomArch} ;
			}
			
			my $outputfilename;
			
			if(length($DomArch.".png") < 255){
			
				$outputfilename = $DomArch.".png";
				
			}else{
				
				my $FName = substr($DomArch,0,200);
				$outputfilename = $FName.".png";
			}
			
			
			`Hist.py -f './.$$-TempPosterior.dat' -o ./$OutputDir/$outputfilename -t '$DomArch ABC vs Max. LH' -x 'Deletion Rate' -y 'Frequency' -l 'Score:$Score\nClade Size:$cladesize Observed:$NumberOfGenomes\nEpsilon:$epsilon' --vline $deletion_rate --column 0`;
			
			my $AverageDifference = 0;
			map{$AverageDifference += ($deletion_rate-$_)}@$PosteriorEstimatesOFDeltionRate;
			
			print OUT $AverageDifference.":".$Score."\n";
			
			#TODO produce a summary of the distance between the ML value and all the ABC estimates
			#TODO plot this as a line plot
			
		}
		
		
	}


$pm->finish if ($maxProcs); # Terminates the child process

close OUT;

}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";

`cat $SIMULATIONSPATH/* > ./RawScoreVsDelRate.dat`;
`EZLines.py -f './RawScoreVsDelRate.dat' -o ./RawScoreVsDelRate.png -t 'ABC vs Max. LH Average' -l 'Epsilon:$epsilon' -x 0 -y 1 --xlab 'Average Difference Between ML and ABC Estimate' --ylab 'Posterior Quantile Score'`;
			

my $TotalToc = Time::HiRes::time;
my $TotalTimeTaken = ($TotalToc - $TotalTic);
my $TotalTimeTakenHours = $TotalTimeTaken/(60*60);

open RUNINFOTIME, ">>ABCHGT_info.$OutputFilename";
print RUNINFOTIME $TotalTimeTaken." seconds\n";
print RUNINFOTIME $TotalTimeTakenHours." hours\n";
close RUNINFOTIME;

print STDERR $TotalTimeTaken." seconds\n";
print STDERR $TotalTimeTakenHours." hours\n";

#-------
__END__


