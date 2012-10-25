#! /usr/bin/env perl

=head1 NAME

RepeatHGTAnalysis<.pl>

=head1 USAGE

RepeatHGTAnalysis.pl [options -v,-d,-h] <ARGS> -t --tree treefile -r --rawdata RawDataFile -d --delrates DelRatesFile -s --simulations SimulationsDataFile -o --output output directory

=head1 SYNOPSIS

A script to perform a repeat of HGT analysis, using stored simualtion values from a pervious study. Consider also calling PostHGTAnalysis so that we can fully drill down and study the results.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2012 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "/home/sardar/bin/perl-libs-custom";

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
use DBI;
use Supfam::Utils;
use Supfam::hgt;
use Math::Random qw(random_uniform_integer);

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my ($treefile) = glob("HGT_tree.*");
my ($rawdatafile) = glob("*-RawData.colsv*");
my ($DelRatesFile) = glob(".DelRates.dat");
my ($SimulationsDataFile) = glob("./RawSimulationDists*");

my $OutputDir = 'RepeatAnalysis';

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t:s" => \$treefile,
           "rawdata|r:s" => \$rawdatafile,
           "output|o:s" => \$OutputDir,
           "delrates|dr:s" => \$DelRatesFile,
           "simulations|s:s" => \$SimulationsDataFile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Create an output directory if it doesn't already exist
mkdir("./$OutputDir");

open RAWSIMS, "<$SimulationsDataFile" or die $!.$?;

my $RawSimsHash = {};
#Slurp in all the results in raw sim data file. Note that you can just concantenate together two data files to get this.
while (my $line = <RAWSIMS>){
	
	chomp($line);
	my ($GenomesCounts,$DomArch,$SimulatedGenomes) = split(/:/,$line);
	
	my $InfoHash = {};
	$RawSimsHash->{$DomArch} = $InfoHash;
	
	$InfoHash->{'Simulations'}=[] unless(exists($InfoHash->{'Simulations'}));
	push(@{$InfoHash->{'Simulations'}},split(/,/,$SimulatedGenomes));
	
	my ($cladesize,$observedgenomes) = split(',',$GenomesCounts);
	
	if(exists($InfoHash->{'ObservedNumberOfGenomes'})){
		
		die "Multiple entries in raw results file for $DomArch disagree on clade size and/or number of observed genomes!\n" unless ($InfoHash->{'ObservedNumberOfGenomes'}[0] == $observedgenomes && $InfoHash->{'ObservedNumberOfGenomes'}[1] eq $cladesize);
		
	}else{
	
		$InfoHash->{'ObservedNumberOfGenomes'}=[$observedgenomes,$cladesize];
	}
}

close RAWSIMS;

#Do a sanity check over input. All values should have same number of simulations. Die if not

#Choose a random DA and measure how many iterations of simualtion there were

my @DomainArchitectures = keys(%$RawSimsHash);
my $Iterations = scalar(@{$RawSimsHash->{$DomainArchitectures[rand(scalar(@DomainArchitectures))]}{'Simulations'}});

#Calculate a score and self-test score for each simulation set. Push onto @scores

my @Scores;

open SCORES, ">./$OutputDir/RepeatScores-RawData.colsv" or die $!.$?;
open SELFTEST, ">./$OutputDir/SelfTest-RawData.colsv" or die $!.$?;

foreach my $DomainArchitecture (@DomainArchitectures){
	
	my ($ObservedNumberOfGenomes,$CladeSize)=@{$RawSimsHash->{$DomainArchitecture}{'ObservedNumberOfGenomes'}};
	
	my $CurrentDAIterations = scalar(@{$RawSimsHash->{$DomainArchitecture}{'Simulations'}});
	die "Differing number of iterations of simulation for DA $DomainArchitecture. Soemthing is very wrong - have you concatenated two very different results files together?\n " unless ($Iterations == $CurrentDAIterations);
	
	my @RawResults = @{$RawSimsHash->{$DomainArchitecture}{'Simulations'}};
	
	my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@RawResults)-1));		
	my $SelftestValue = $RawResults[$selftest_index]; # A single uniform random simulation value
	
	my $Distribution = {};
	map{$Distribution->{$_}++}@RawResults;
	
	my $PosteriorQuantile = calculatePosteriorQuantile($ObservedNumberOfGenomes,$Distribution,$Iterations+1,$CladeSize);
	#Iterations plus one is a consequence of the fact that our observed value is being treated as though it were a 'simulated' value
	
	#$SelftestValue +=1 if ($SelftestValue == 0);
	#This is an error for previous scipts, left in here so that we can switch it on the easily demonstrate the error source
	
	$Distribution->{$SelftestValue}--;
	my $SelfTestQuantile = calculatePosteriorQuantile($SelftestValue,$Distribution,$Iterations,$CladeSize);
	
	#Print scores to files
	
	print SCORES $DomainArchitecture.':'.$PosteriorQuantile."\n";
	print SELFTEST $DomainArchitecture.':'.$SelfTestQuantile."\n";
	
	push(@Scores,$PosteriorQuantile);
}

close SCORES;
close SELFTEST;

`Hist.py -f "./$OutputDir/RepeatScores-RawData.colsv" -o "./$OutputDir/RepeatScores.png" -t "Histogram of Cumulative p-Values" -x "P(Nm < nr)" -y "Frequency"` ;
`Hist.py -f "./$OutputDir/SelfTest-RawData.colsv" -o ./$OutputDir/SelfTest.png -t "Histogram of Self-Test Cumulative p-Values" -x "P(Nm < nm)" -y "Frequency"` ;


print STDERR "PostHGTAnalysis.pl -t $treefile --rawdata ./$OutputDir/RepeatScores-RawData.colsv --delrates $DelRatesFile --simulations $SimulationsDataFile -o $OutputDir/InDepthAnalysis/ -w 3000\n" if($debug); 

`PostHGTAnalysis.pl -t $treefile --rawdata "./$OutputDir/RepeatScores-RawData.colsv" --delrates $DelRatesFile --simulations $SimulationsDataFile -o "./$OutputDir/InDepthAnalysis" -w 3000`;

#Calculate overall left-right symmetry
#Calculate left-right symmetry in 0.1-0.9

my $NumberLHSscores = grep{$_ < 0.5}@Scores;
my $NumberRHSscores = grep{$_ > 0.5}@Scores;

my $Asymmetry = 100*($NumberRHSscores - $NumberLHSscores)/scalar(@Scores);

$NumberLHSscores = grep{$_ < 0.5 && $_ > 0.1}@Scores;
$NumberRHSscores = grep{$_ > 0.5 && $_ < 0.9}@Scores;

my $EightyPercentAsymmetry =  100*($NumberRHSscores - $NumberLHSscores)/($NumberLHSscores + $NumberRHSscores);

open ASYM, ">./$OutputDir/Asymmetry.txt" or die $!.$?;
print ASYM "Total Asymmetry: ".$Asymmetry."%\n";
print ASYM "Mid-Eighty Percent".$EightyPercentAsymmetry."%";
close ASYM;



__END__
