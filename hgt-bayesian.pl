#! /usr/bin/env perl

#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;
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

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $count =0;

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'HGTResults';
my $TreeFile;
my $maxProcs = 0;
my $test = 'n';
my $FalseNegativeRate = 0.05;
my $completes = 'y'; #A flag to specify whether or not to include architectures including _gap_ assignmenets
my $check = 'y'; #Perform a sanity check on the tree? This should be 'y' unless under extreme circumstances


#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$OutputFilename,
           "tree|i:s" => \$TreeFile,
           "processor_cores|p:i" => \$maxProcs,
           "completes|comp:s" => \$completes,
           "fals_negative_rate|fnr:f" => \$FalseNegativeRate,
           "check|c:s" => \$check,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------

#my $dbh = dbConnect();

my $dbh = DBI->connect("DBI:mysql:superfamily;localhost"
                                        ,''
                                        , undef
                                        ,{RaiseError =>1}
                                    ) or die ;
#These two alternate DBI configurations make it easy to run the script on the amazon EC2 cloud
# Connect to database

`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';
#Path to a piece of RAM that we can write to. This could be on hard disk, but on *nix systems /dev/shm is a shared RAM disk folder. We want a temporary folder in her that will be cleaned up on Prgram exit

# Make a temporary directory for output data 
my $RAWPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;
my $DELSPATH= File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;

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

print STDERR "Number of genomes in tree: ".scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})."\n";

print STDERR "False Negative Rate:".$FalseNegativeRate."\n";
#print STDERR "Standard Error of mean estimators will be approx: ".int(100*0.5*sqrt(scalar(@{$TreeCacheHash->{$root}{'Clade_Leaves'}})/$Iterations))."% \n";
#Standard error (statistics) is defined as 0.5*sqrt(no_genomes/no_samples(iterations)).

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

my $lencombquery = join ("' or len_comb.genome='", st recent figures suggest that we get through 47m gallons of the gloopy stu@TreeGenomes); $lencombquery = "(len_comb.genome='$lencombquery')";# An ugly way to make the query run, but perl DBI only alows for a single value to occupy a wildcard

my $DomCombGenomeHash = {};

if($completes eq 'n'){
	$sth = $dbh->prepare("SELECT len_comb.genome,comb_index.comb FROM len_comb JOIN comb_index ON len_comb.comb_id=comb_index.id WHERE $lencombquery AND comb_index.comb != '_gap_';");
}elsif($completes eq 'y'){
	$sth = $dbh->prepare("SELECT len_comb.genome,comb_index.comb FROM len_comb JOIN comb_index ON len_comb.comb_id=comb_index.id WHERE $lencombquery AND comb_index.comb NOT LIKE '%_gap_%';");
}else{
	die "Inappropriate flag for whether or not to include architectures containign _gap_";
}

$sth->execute();

	$tic = Time::HiRes::time;

while (my ($genomereturned,$combreturned) = $sth->fetchrow_array() ){

	die "$combreturned\n" if($combreturned =~ /_gap_/ && $completes eq 'y'); # sanity check, will likely remove from final version of script

	$DomCombGenomeHash->{$combreturned} = {} unless (exists($DomCombGenomeHash->{$combreturned}));
	$DomCombGenomeHash->{$combreturned}{$genomereturned}++;
}

	 $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Dom Arch hash:".($toc-$tic)."seconds\n";

my $DomArchs = [];
@$DomArchs = keys(%$DomCombGenomeHash);
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
		
	open OUT, ">$RAWPATH/$OutputFilename".$$.".-RawData.colsv" or die "Can't open file $RAWPATH/$OutputFilename".$!;
	open DELS, ">$DELSPATH/DelRates".$$.".dat" or die "Can't open file $DELSPATH/DelRates".$!;

	foreach my $DomArch (@$ArchsListRef){
		
		my ($CladeGenomes,$NodesObserved);
			
		my $HashOfGenomesObserved = $DomCombGenomeHash->{$DomArch};
		@$NodesObserved = keys(%$HashOfGenomesObserved);
		my @NodeIDsObserved = map{$NodeNameIDHash->{$_}}@$NodesObserved ; 
		#Get the node IDs as the follwoing function doesn't work with the raw node tags
	
		my $MRCA; # Most recent common ancestor
		my $deletion_rate;
		my ($dels, $time);
		
		unless(scalar(@NodeIDsObserved) == 1){
		
			$MRCA = $tree->get_lca(-nodes => \@NodeIDsObserved); #Find most Recent Common Ancestor	
			($dels, $time) = DeletedJulian($MRCA,0,0,$HashOfGenomesObserved,$TreeCacheHash,$root); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree		
			$deletion_rate = $dels/$time;
			
		}else{
			
			$deletion_rate = 0;	
			$MRCA = $NodeIDsObserved[0] ; #Most Recent Common Ancestor
		}
		
		@$CladeGenomes = @{$TreeCacheHash->{$MRCA}{'Clade_Leaves'}}; # Get all leaf genomes in this clade (as NODEIDs)
						
		print DELS "$DomArch:$deletion_rate\n";
		#print "$DomArch:$deletion_rate\n";

		my $DomArchLH; #This will be the likelihood of the domain archtecture phylogeteic spread given the oberserved data.
		
		if($deletion_rate > 0){
									
			$DomArchLH = RandomModel($CladeGenomes,\@NodeIDsObserved,$TreeCacheHash,$deletion_rate,$MRCA);

		}else{
			
			$DomArchLH = ('NULL');
		}
		
		#no warnings 'numeric';
		
		
		## ----- OUTPUT --- ##
		
		my $LogLH = log($DomArchLH) unless($DomArchLH eq 'NULL' || $DomArchLH == 0);
	#	print STDERR $count++." $DomArch : $DomArchLH\n "if ($DomArchLH == 0 && $DomArchLH ne 'NULL');
		
		print OUT "$DomArch:$LogLH\n" unless($DomArchLH eq 'NULL' || $DomArchLH == 0);
	
	}	

	close OUT;
	close DELS;
	
$pm->finish if ($maxProcs); # Terminates the child process
}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";

`cat $DELSPATH/* > ./.DelRates.dat`;
`cat $RAWPATH/* > ./$OutputFilename-RawLH.colsv`;


`Hist.py -f "./$OutputFilename-RawLH.colsv" -o $OutputFilename.png -t "Histogram of LogLH Values" -x "Score" -y "Frequency"	`;
#`Hist.py -f "./.DelRates.dat" -o ParDelRates.png -t "Histogram of DeletionRates" -x "Deletion Rate" -y "Frequency"	-l Deletions`;
# Plot a couple of histograms for easy inspection of the data

open PLOT, ">./.ParPlotCommands.txt" or die $!;
print PLOT "Hist.py -f ./$OutputFilename-RawLH.colsv -o $OutputFilename.png -t 'Histogram of LogLH Values' -x 'Score' -y 'Frequency'\n\n\n";
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

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

