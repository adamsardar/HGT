#! /usr/bin/perl -w

=head1 NAME

HGTConvergnce<.pl>

=head1 USAGE

HGTConvergnce.pl [options -v,-d,-h] -f INPUTFILE -o OUTPUTFILE

=head1 SYNOPSIS

A script to test for convergence of asymetric differences in densities of the resulting distributions
from HGT scripts.

INPUTFILE is the (possibly concatenated) results from hgt-restart.pl. 

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use List::Util qw[min max];
use Term::ProgressBar;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $inputfile;
my $outputfile = 'convergence.png';

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "file|f:s" => \$inputfile,
        	"output|o:s" => \$outputfile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my $RawResultsHash = {};
# This will contain all results of simulation runs: DomArch => @[data]
my $Iterations = [];
#A list of the number of iterations that each dom arch ahs been run over
my $TreeDataHash = {};
#Details regarding the tree used to generate this data

open INPUT, "<$inputfile" or die "$? $!";

#Rack open the output file of a simulation and store in a hash

while (my $line = <INPUT>){
	
	my ($TreeDetails,$DomArch,$Results) = split(':',$line);
	my ($CladeSize,$No_Genomes) = split(',',$TreeDetails);
	
	
	unless (exists($TreeDataHash->{$DomArch})){
		$TreeDataHash->{$DomArch}={};
		$TreeDataHash->{$DomArch}{'CladeSize'} = $CladeSize;
		$TreeDataHash->{$DomArch}{'Number_Of_Genomes_Observed'} = $No_Genomes;
	}else{
		die "Different data appended to file, there's an error somehwere: a different size/type of tree was used to generate this" if ($TreeDataHash->{$DomArch}{'CladeSize'} != $CladeSize || $TreeDataHash->{$DomArch}{'Number_Of_Genomes_Observed'} != $No_Genomes);
		#Some data validity checking
	}
	
	$RawResultsHash->{$DomArch}=[] unless(exists($RawResultsHash->{$DomArch}));
	push(@{$RawResultsHash->{$DomArch}},split(',',$Results));	
}

close INPUT;

my $CladeSizeArray = [];
my $CrossEntropyHash = {}; #This will be a hash of all the individual cross entropy terms expectation terms

my $DistProbArray = [];
my $ModelProbArray = [];

foreach my $DA (keys(%$TreeDataHash)){
	
	#Find two values - $Distmle, which is the maximum likelihood estimator of the underlying distribution and $Modellh which is the model (simulated) likelihood
	my $Distmle;
	
	my ($CladeSize,$No_Genomes) = ($TreeDataHash->{$DA}{'CladeSize'},$TreeDataHash->{$DA}{'Number_Of_Genomes_Observed'});
	push(@$CladeSizeArray,$CladeSize);
	
	$Distmle = $No_Genomes/$CladeSize;
	
	next if($CladeSize < 40);
	
	my ($ModelDistributions,$ModelProbabilities) = ({},{});
	my $Iterations = scalar(@{$RawResultsHash->{$DA}});

	map{$ModelDistributions->{$_}++}@{$RawResultsHash->{$DA}}; #A raw count of frequency of a particular outcome of a simulation
	map{$ModelProbabilities->{$_}=($ModelDistributions->{$_}/$Iterations)}(keys(%$ModelDistributions)); #Normalise the above hash so as to form probabilities

	my $CrossEnt = 0;
	
	foreach my $ModelOutcome (keys(%$ModelProbabilities)){ #$ModelOutcome is the number of genomes
		
		my $ModelProb = $ModelProbabilities->{$ModelOutcome};
		my $ModelGenomeProb = $ModelOutcome/$CladeSize;
		
		my $CrossEntContribution;
		$CrossEntContribution = $Distmle*(log(1/$ModelGenomeProb)/log(2)) + (1-$Distmle)*log(1/(1-$ModelGenomeProb))/log(2) unless( $ModelGenomeProb == 0||$ModelGenomeProb==1);
		$CrossEntContribution = $Distmle*(log(1/$ModelGenomeProb)/log(2)) if($ModelGenomeProb ==1);
		$CrossEntContribution = (1-$Distmle)*log(1/(1-$ModelGenomeProb))/log(2) if($ModelGenomeProb==0);
		
		$CrossEnt+=$CrossEntContribution*$ModelProb;
	}
	
	#E[H(p,M)] = SUMj -Rj(Pi*lb{Mj})+(1-Pi)*lb{(1-Mj)}), where Rj is the outcome probability, Mj is the liklihood of a random genome possesing domain Arch Ai and Pi is the MLE of that same event
		
	$CrossEntropyHash->{$DA}=$CrossEnt;
}

#Output file of each dom arch KL dist

open OUT, "> CrossEntropyColsv.dat";
		
foreach my $DomainArchitecture (keys(%$CrossEntropyHash)){

	my $Cross = $CrossEntropyHash->{$DomainArchitecture};
	print OUT "$DomainArchitecture:$Cross\n";
}

close OUT;

#Output hist of CrossEnts (labelled with total entropy at top)

`Hist.py -f "CrossEntropyColsv.dat" -o CrossEntropy.png -t "Histogram of Expected Cross Entropies " -x "Cross Entropy (bits)" -y "Freq."	-b 100`;

#Output hist of clade sizes

open CLADES, "> CladeSizeColsv.dat";
		
foreach my $CladeSize (@$CladeSizeArray){

	print CLADES "Clade:$CladeSize\n";
}

close CLADES;


#open DISTS, "> DistsColsv.dat";
#		
#foreach my $Dist (@$DistProbArray){
#
#	print DISTS "Dist:$Dist\n";
#}
#
#close DISTS;

#open MODEL, "> ModelColsv.dat";
#		
#foreach my $Model (@$ModelProbArray){
#
#	print MODEL "Model:$Model\n";
#}
#
#close MODEL;

`Hist.py -f "CladeSizeColsv.dat" -o CladeSizeColsv.png -t "Histogram of Clade Size" -x "Clade Size" -y "Freq."`;
#`Hist.py -f "DistsColsv.dat" -o Dists.png -t "Histogram of MLE of underlying distribution" -x "P(DA in genome)" -y "Freq."`;
#`Hist.py -f "ModelColsv.dat" -o ModelColsv.png -t "Histogram of Model Probability Dist" -x "P(DA in genome)" -y "Freq."`;


__END__

