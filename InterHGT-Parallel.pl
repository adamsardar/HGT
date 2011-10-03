#! /usr/bin/perl -w
# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
 if(-d "$ENV{HOME}/workspace/Oates/lib"){use lib "$ENV{HOME}/workspace/Oates/lib/"};
 if(-d "$ENV{HOME}/AmazonModules" && ! -d "$ENV{HOME}/workspace/Oates/lib"){use lib "$ENV{HOME}/AmazonModules/"};

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::Utils;
use Supfam::hgt;
use Supfam::SQLFunc;

use File::Temp;
use Time::HiRes;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

use Parallel::ForkManager;

use Sys::Hostname;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'IntraHGTResults';
my $TreeFile="./primates.tree";
my $maxProcs = 0;
my $test = 'n';
my $FalseNegativeRate = 0.05;
my $completes = 'y'; #A flag to specify whether or not to include architectures including _gap_ assignmenets
my $Iterations = 500;
my $model = 'Julian';
my $store = '0';
my $outgroup;


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
           "fals_negative_rate|fnr:i" => \$FalseNegativeRate,
           "model|m:s" => \$model,
           "store|s:i" => \$store,
           "outgroup|out=s" => \$outgroup,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------

my $host = hostname;

my $dbh;

if($host eq 'luca'){
	$dbh = dbConnect();
}else{
 $dbh = DBI->connect("DBI:mysql:superfamily;localhost"
                                        ,''
                                        , undef
                                        ,{RaiseError =>1}
                                    ) or die ;
#These two alternate DBI configurations make it easy to run the script on the amazon EC2 cloud
}

# Connect to database
`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';
#Path to a piece of RAM that we can write to. This could be on hard disk, but on *nix systems /dev/shm is a shared RAM disk folder. We want a temporary folder in her that will be cleaned up on Prgram exit

# Make a temporary directory for output data 
my $RAWPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $HTMLPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $DELSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $CONFIDENCEPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;


#my $TEMPDISTPATH= File::Temp->newdir( DIR => './' , CLEANUP => 1) or die $!;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my $tic = Time::HiRes::time; 

my $TreeCacheHash = {}; # A massive limitation of this script is the use of BioPerl TreeIO which is SLOOOOW. This is a lookup hash to speed things up.

my $input = new Bio::TreeIO(-file   => "$TreeFile",
                            -format => "newick") or die $!;
                            
my $tree = $input->next_tree;
my $root = $tree->get_root_node;

#Read in and initialise tree

## Initialise all the tree values that we might need in this hash to minimise calls to BioTree.



$TreeCacheHash->{$root}={};
my @RootDescendents = $root->get_all_Descendents;
$TreeCacheHash->{$root}{'all_Descendents'}=\@RootDescendents;
map{$TreeCacheHash->{$_}={}}@{$TreeCacheHash->{$root}{'all_Descendents'}};
#Bootstrap the BioPerl internal nodeids as keys to Cache Hash


foreach my $node ($root->get_all_Descendents,$root) {

	#each_descendant
	$TreeCacheHash->{$node}{'each_Descendent'} = [];
	push(@{$TreeCacheHash->{$node}{'each_Descendent'}},$node->each_Descendent);
	
	#alldescendents
	my @AllDescendents = $node->get_all_Descendents;
	$TreeCacheHash->{$node}{'all_Descendents'} =  \@AllDescendents;

	#Clade_leaves
	my @CladeLeaves = grep{$_->is_Leaf == 1}@{$TreeCacheHash->{$node}{'all_Descendents'}};
	$TreeCacheHash->{$node}{'Clade_Leaves'} = \@CladeLeaves;
	
	#Branch Length
	my $Branch = $node->branch_length;
	$TreeCacheHash->{$node}{'branch_length'} = $Branch ;
	$TreeCacheHash->{$node}{'branch_length'} = 0 if ($node = $root);
	
	#Total Branch Lengths in tree - I'd love a section in here: maybe in the future

}


#is_Leaf Entries
map{$TreeCacheHash->{$_}{'is_Leaf'} = 0}@{$TreeCacheHash->{$root}{'all_Descendents'}};
map{$TreeCacheHash->{$_}{'is_Leaf'} = 1}@{$TreeCacheHash->{$root}{'Clade_Leaves'}};

my $toc = Time::HiRes::time; 
print STDERR "Time taken to build the Tree Cache hash:".($toc-$tic)."seconds\n";

print STDERR "No of iterations per run is: $Iterations\n";

#--Check Input Tree------------------------

my $NodeNameIDHash = {}; #This will be a hash to map from NodeName to the BioPerl NodeID $NodeNameIDHash->{NodeName}=NodeID
my @TreeGenomes = map{$_->id}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
map{$NodeNameIDHash->{$_->id}=$_}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; #Populate $NodeNameIDHash

my $genquery = join ("' or genome.genome='", @TreeGenomes); $genquery = "(genome.genome='$genquery')";# An ugly way to make the query run

my $sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();

while (my @query = $sth->fetchrow_array() ) {
		
	die "Genome: $query[0] is in the tree but should not be!\n"  unless ($query[0] eq "y" && $query[1] eq '' );
} #A brief bit of error checking - make sure that the tree doesn't contain any unwanted genoems
			
	#---------Get a list of domain archs present in tree----------------------

my $lencombquery = join ("' or len_comb.genome='", @TreeGenomes); $lencombquery = "(len_comb.genome='$lencombquery')";# An ugly way to make the query run, but perl DBI only alows for a single value to occupy a wildcard

my $IngroupDomCombGenomeHash = {};
my $OutgroupDomCombGenomeHash = {};

		 $tic = Time::HiRes::time;

if($completes eq 'n'){
	$sth = $dbh->prepare("SELECT len_comb.genome,len_comb.comb,genome.domain FROM len_comb JOIN genome ON len_comb.genome=genome.genome WHERE ($lencombquery) AND len_comb.comb  != '_gap_';");
}elsif($completes eq 'y'){
	$sth = $dbh->prepare("SELECT len_comb.genome,len_comb.comb,genome.domain FROM len_comb JOIN genome ON len_comb.genome=genome.genome WHERE ($lencombquery) AND len_comb.comb  NOT LIKE '%_gap_%';");
}else{
	die "Inappropriate flag for whether or not to include architectures containign _gap_";
}

$sth->execute();

my $DomArchHash ={};

while (my ($genomereturned,$combreturned,$life_domain) = $sth->fetchrow_array() ){
	
	if($life_domain eq $outgroup){
	
			$OutgroupDomCombGenomeHash->{$combreturned} = {} unless (exists($OutgroupDomCombGenomeHash->{$combreturned}));
			$OutgroupDomCombGenomeHash->{$combreturned}{$genomereturned}++;
				
	}else{
		
			$IngroupDomCombGenomeHash->{$combreturned} = {} unless (exists($OutgroupDomCombGenomeHash->{$combreturned}));
			$IngroupDomCombGenomeHash->{$combreturned}{$genomereturned}++;
	}

	die "$combreturned\n" if($combreturned =~ /_gap_/ && $completes eq 'y'); # sanity check, will likely remove from final version of script
	$DomArchHash->{$combreturned}++;
}

die "Poor choice of outgroup or tree" unless (keys(%$OutgroupDomCombGenomeHash));

	 $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Dom Arch hashes:".($toc-$tic)."seconds\n";

my $DomArchs = [];
@$DomArchs = keys(%$DomArchHash);
#These are all the unique domain architectures

dbDisconnect($dbh) ; 

my $NoOfForks = $maxProcs;
$NoOfForks = 1 unless($maxProcs);

my $remainder = scalar(@$DomArchs)%$NoOfForks;
my $binsize = (scalar(@$DomArchs)-$remainder)/$NoOfForks;

my $ForkJobsHash = {};

for my $i (0 .. $NoOfForks-1){
	
	my @temp = @{$DomArchs}[$i*$binsize .. ($i+1)*$binsize-1];
	$ForkJobsHash->{$i}=\@temp;
}

push(@{$ForkJobsHash->{0}},@{$DomArchs}[($NoOfForks)*$binsize .. ($NoOfForks)*$binsize+$remainder-1]) if ($remainder);

#----------------------------------------------------

print STDERR "No Dom Archs in job batch is approx: ".$binsize."\n";

#Main-loop-------------------------------------------

print STDERR "Total No Of Dom Archs: ".@$DomArchs."\n";

my $pm = new Parallel::ForkManager($maxProcs) if ($maxProcs);# Initialise

foreach my $fork (0 .. $NoOfForks-1){
	
	my $ArchsListRef = $ForkJobsHash->{$fork};
			
	# Forks and returns the pid for the child:
	if ($maxProcs){$pm->start and next};
	
	#$$ is the PID of the process
		
	my $CachedResults = {}; #Allow for caching of distributions after Random Model to speed things up
		
	open HTML, ">$HTMLPATH/$OutputFilename".$$.".html" or die "Can't open file $HTMLPATH/$OutputFilename".$!;
	open OUT, ">$RAWPATH/$OutputFilename".$$.".-RawData.colsv" or die "Can't open file $RAWPATH/$OutputFilename".$!;
	open DELS, ">$DELSPATH/DelRates".$$.".dat" or die "Can't open file $DELSPATH/DelRates".$!;
	open CONFIDENCE, ">$CONFIDENCEPATH/StandardErrorEstimate".$$.".dat" or die "Can't open file $CONFIDENCEPATH/".$!;
					
	foreach my $DomArch (@$ArchsListRef){
		
		if($IngroupDomCombGenomeHash->{$DomArch} && $OutgroupDomCombGenomeHash->{$DomArch}){
			#i.e. if the dom arch has been seen in both groups
				
			my ($CladeGenomes,$NodesObserved);
					
				my $HashOfGenomesObserved = {};
				map{$HashOfGenomesObserved->{$_}++}(keys(%{$IngroupDomCombGenomeHash->{$DomArch}}),keys(%{$OutgroupDomCombGenomeHash->{$DomArch}}) );
				#No value should be more than 1, else there is a genome that both in the ingroup AND the outgroup! There is a check for this in the deleted function below -> see hgt.pm

				@$NodesObserved = keys(%$HashOfGenomesObserved);
				my @NodeIDsObserved = map{$NodeNameIDHash->{$_}}@$NodesObserved ; 
				#Get the node IDs as the follwoing function doesn't work with the raw node tags
			
				my $MRCA;
				my $deletion_rate;
				my ($dels, $time) = (0,0);
				
				unless(scalar(@NodeIDsObserved) == 1){
				
				$MRCA = $tree->get_lca(-nodes => \@NodeIDsObserved) ; #Find most Recent Common Ancestor	
				
				 if($model eq 'Julian'){
						
						($dels, $time) = DeletedJulian($MRCA,0,0,$HashOfGenomesObserved,$TreeCacheHash,$root); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree	
						
					}elsif($model eq 'sardar'){
						die "Significant work needed on this portion of the script";
						($dels, $time) = Deleted($tree,$MRCA,0,0,$NodesObserved,0); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree
						
						my $EvolutionaryTime = 0; # Evolutionary time from ancestor
						$EvolutionaryTime = $TreeCacheHash->{$MRCA}{'branch_length'} unless($MRCA eq $root);
						$time -= $EvolutionaryTime;
						#This little hack is to correct for the fact that the distance between the original clade parent and its ancestor has been included. (Frustrating, but this allows it to be a recursive function)
		
					}else{
						die "Inappropriate model selected";
					}
					
				$deletion_rate = $dels/$time;
				
				}else{
					$deletion_rate = 0;	
					$MRCA = $NodeIDsObserved[0] ; #Most Recent Common Ancestor
				}
				
				@$CladeGenomes = @{$TreeCacheHash->{$MRCA}{'Clade_Leaves'}}; # Get all leaf genomes in this clade	
				
				print DELS "$DomArch:$deletion_rate\n";
				#print "$DomArch:$deletion_rate\n";
				
				my $subtree = Bio::Tree::Tree->new(-root => $MRCA, -nodelete => 1);
				
				my ($selftest,$distribution);
				
				if($deletion_rate > 0){	
			
					unless($CachedResults->{"$deletion_rate:@$CladeGenomes"} && $store){
								
						if($model eq 'Julian'){
											
							($selftest,$distribution) = RandomModelJulian($subtree,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);
																							
						}elsif($model eq 'sardar'){
							
							die "Serious work needed on Sardar Model";
							($selftest,$distribution) = RandomModel('ERROR',$FalseNegativeRate,$Iterations,$deletion_rate,'ERROR');
						}else{
							die "No appropriate model selected";
						}
						
					$CachedResults->{"$deletion_rate:@$CladeGenomes"} = [$selftest,$distribution];		
					
					}else{
						($selftest,$distribution) = @{$CachedResults->{"$deletion_rate:@$CladeGenomes"}};
					}
					
				}else{
					
					($selftest,$distribution) = ('NULL',{});
				}
		
		#-------------- Output
				 
				my $ResultsHash={};
				my $Average = 0;
						
				unless($selftest eq 'NULL'){
						
				print CONFIDENCE $DomArch.":".int(100*0.5*sqrt(scalar(@$CladeGenomes)/$Iterations))."\n";
						
			#	EasyDump("./$TEMPDISTPATH/ParDistribution$DomArch.dump",[$DomArch,$distribution]); ######
					
					#---- Calculate discrete distance from median
					my $CumulativeCount = 0;
					
						for my $NumberOfGenomes (0 .. scalar(@$CladeGenomes)){
							
							if($NumberOfGenomes == scalar(@$NodesObserved)){
								
								my $GenomesSeenSoFar = $CumulativeCount; #This is a count of the number of model (simulation) genomes seen so far
								
									if(exists($distribution->{$NumberOfGenomes})){
										
											my $ModelFrequency = $distribution->{$NumberOfGenomes}; #ModelFrequency is the number of times that the number of actual genomes observed has been seen in simulation					
											$Average +=  ($GenomesSeenSoFar+(1/2)*($ModelFrequency)); # This is the sum of an arithmetic preogression starting at $GenomesSeeSoFar and ending at $GenomesSeeSoFar + ($ModelFrequency -1) then divided by (1/$ModelFrequency)
				
											for(1 .. $ModelFrequency){$ResultsHash->{$GenomesSeenSoFar++} = (1/$ModelFrequency);}
																
									}else{
										
											$Average += $CumulativeCount;
											$ResultsHash->{$CumulativeCount} += 0.5;
											$ResultsHash->{$CumulativeCount+1} += 0.5;
										
									}last;
							}
							
							 $CumulativeCount += $distribution->{$NumberOfGenomes} if(exists($distribution->{$NumberOfGenomes})); # Cumulative count is a sum of the frequency of genome observation up to this point
						}
					
						
						# ---- I have to honest and say that I don't relaly know what the above code does ... I think that this HIDEOUS piece of code is going to be used in calculating distance from the median. It sums the models seen up until the bin containing the genome of interest.
			
				        print HTML "<a href=http://http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/maketree.cgi?genomes=";
				        print HTML join(',', @$NodesObserved);
				        print HTML ">$DomArch</a> Score: $CumulativeCount<BR>\n";
				        
				        #------
				
					print OUT "$DomArch:$CumulativeCount\n";
					#The output value 'Score:' is the number of model run genomes seen up to the bin containing the observed genome number	
				}	
				
		
		}
	
}

	close HTML;
	close OUT;
	close DELS;
	close CONFIDENCE;
$pm->finish if ($maxProcs); # Terminates the child process
}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";

`cat $HTMLPATH/* > ./$OutputFilename.html`;
`cat $DELSPATH/* > ./.DelRates.dat`;
`cat $RAWPATH/* > ./$OutputFilename-RawData.colsv`;
`cat $CONFIDENCEPATH/* > ./$OutputFilename-STDERRESTIMATES.colsv`;

#`cat ./$TEMPDISTPATH/* > ./ParDists.dump`;

if(-e '~/bin/Hist.py'){
`Hist.py -f "./$OutputFilename-RawData.colsv" -o $OutputFilename.png -t "Histogram of Scores" -x "Score" -y "Frequency"	`;
`Hist.py -f "./.DelRates.dat" -o ParDelRates.png -t "Histogram of DeletionRates" -x "Deletion Rate" -y "Frequency"	-l Deletions`;
# Plot a couple of histograms for easy inspection of the data
}

open PLOT, ">./.ParPlotCommands.txt" or die $!;
print PLOT "Hist.py -f ./$OutputFilename-RawData.colsv -o $OutputFilename.png -t 'Histogram of Scores' -x 'Score' -y 'Frequency'\n\n\n";
print PLOT "Hist.py -f './.DelRates.dat' -o 'ParDelRates.png' -t 'Histogram of DeletionRates' -x 'Deletion Rate' -y 'Frequency'	-l 'Deletions'";
close PLOT;
#-------
__END__


=head1 NAME

I<.pl>

=head1 USAGE

 .pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to...

=head1 AUTHOR

B<Joe Bloggs> - I<Joe.Bloggs@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2010 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

