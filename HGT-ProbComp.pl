#! /usr/bin/env perl

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

#####

my $ScoresLT = {}; #Less than
my $ScoresLTET = {}; #Less Than Or Equal To
my $ScoresExpectedLT = {}; #Expected scores Less Than 
my $ScoresExpectedLTET = {}; #Expected scores Less Than Or Equal To


foreach my $DomArch (keys(%$RawResultsHash)) {
	
	my $CladeSize = $TreeDataHash->{$DomArch}{'CladeSize'};
	my $NoGenomesObserved = $TreeDataHash->{$DomArch}{'Number_Of_Genomes_Observed'};
	
	my $distribution = {};
	my @Simulations = @{$RawResultsHash->{$DomArch}};
	
	map{$distribution->{$_}++}@Simulations;
	
	my $CumulativeCount= 0;
	my @CumulativeDistribution; #(P(Nm<nr)) nr (number of genomes in reality) is the index, Nm is the random variable
	
	#Create cumulative distibution list
	for my $NumberOfGenomes (0 .. $CladeSize){

		$CumulativeCount += $distribution->{$NumberOfGenomes} if(exists($distribution->{$NumberOfGenomes})); # Cumulative count is a sum of the frequency of genome observation up to this point
		$CumulativeDistribution[$NumberOfGenomes]=$CumulativeCount;	
	}
	
	my $ExpectedValueLT=0;
	my $ExpectedValueLTET=0;
	
	for my $NumberOfGenomes (0 .. $CladeSize){

		$ExpectedValueLT+=(($distribution->{$NumberOfGenomes})/@Simulations)*($CumulativeDistribution[$NumberOfGenomes-1]/@Simulations) if(exists($distribution->{$NumberOfGenomes}));
		$ExpectedValueLTET+=(($distribution->{$NumberOfGenomes})/@Simulations)*($CumulativeDistribution[$NumberOfGenomes]/@Simulations) if(exists($distribution->{$NumberOfGenomes}));
	}		
	
	$ScoresLT->{$DomArch} = ($CumulativeDistribution[$NoGenomesObserved-1])/@Simulations;#Cumlative probability of less genomes in model simulations than in reality
	$ScoresLTET->{$DomArch} = 1-($CumulativeDistribution[$NoGenomesObserved])/@Simulations;	#Cumlative probability of less than or equal genomes in model simulations than in reality
	$ScoresExpectedLT->{$DomArch}=$ExpectedValueLT;
	$ScoresExpectedLTET->{$DomArch}=$ExpectedValueLTET;

}

open OUTLT, "> .LessThan$outputfile.Colsv";
open OUTLTET, "> .GreaterThan$outputfile.Colsv";
open EXPLT, "> .ExpectedLessThan$outputfile.Colsv";
open EXPLTET, "> .ExpectedGreaterThan$outputfile.Colsv";

foreach my $FinalDomArch (keys(%$ScoresLT)){

		my $scoreLT = $ScoresLT->{$FinalDomArch};
		my $scoreLTET = $ScoresLTET->{$FinalDomArch};
		
		my $expLT = $ScoresExpectedLT->{$FinalDomArch};
		my $expLTET = $ScoresExpectedLTET->{$FinalDomArch};
		
		print OUTLT "$FinalDomArch:$scoreLT\n";
		print OUTLTET "$FinalDomArch:$scoreLTET\n";
		print EXPLT "$FinalDomArch:$expLT\n";
		print EXPLTET "$FinalDomArch:$expLTET\n";
}

close OUTLT;
close OUTLTET;
close EXPLT;
close EXPLTET;

`Hist.py -f ".LessThan$outputfile.Colsv" -o LessThan.png -t "Histogram of P(No. model genomes < No observed genomes)" -x "P(Nm < nr)" -y "Frequency"	`;
`Hist.py -f ".GreaterThan$outputfile.Colsv" -o GreaterThan.png -t "Histogram of P(No. model genomes > No observed genomes)" -x "P(Nm > nr)" -y "Frequency"	`;

`Hist.py -f ".ExpectedLessThan$outputfile.Colsv" -o ExpLessThan.png -t "Histogram of E[P(No. model genomes < No observed genomes)]" -x "E[P(Nm < Nr)]" -y "Frequency"	`;
`Hist.py -f ".ExpectedGreaterThan$outputfile.Colsv" -o ExpGreaterThan.png -t "Histogram of E[P(No. model genomes <= No observed genomes)]" -x "E[P(Nm <= Nr)]" -y "Frequency"`;



__END__

