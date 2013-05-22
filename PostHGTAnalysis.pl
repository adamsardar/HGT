#!/usr/bin/env perl

=head1 NAME

PostHGTAnalysisI<.pl>

=head1 USAGE

 PostHGTAnalysis.pl [options -v,-d,-h] <ARGS> -t --tree treefile -r --rawdata RawDataFile -d --delrates DelRatesFile -s --simulations SimulationsDataFile -e --Examples #OfExamples -w --width WidthOfTreePlots 


=head1 SYNOPSIS

A script to perform some useful summary analysis after a large HGT run. This will produce tree diagrams with overlays of domain architectures and study the left-right symmetry.

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
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Carp;
use Supfam::Utils;


# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my ($treefile) = glob("HGT_tree.*");

my @RawFiles = glob("*-RawData.colsv*");
@RawFiles = grep{$_ !~ m/SelfTest/i}@RawFiles;
my ($rawdatafile) = $RawFiles[0];
print STDERR "Using raw data file: ".$rawdatafile."\n";
my ($DelRatesFile) = glob(".DelRates.dat");
my ($SimulationsDataFile) = glob("./RawSimulationDists*");
my $NumberOfExamples = 20;
my $Width = 1500;

my $OutputDir = 'IndepthStudyOfDAs';
my $SummaryDir = 'SummaryStatistics';

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t:s" => \$treefile,
           "rawdata|r:s" => \$rawdatafile,
           "delrates|dr:s" => \$DelRatesFile,
           "simulations|s:s" => \$SimulationsDataFile,
           "Examples|e:i" => \$NumberOfExamples,
           "width|w:i" => \$Width,
           "outputdir|o:s" => \$OutputDir,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Create an output directory if it doesn't already exist
mkdir("./$SummaryDir");
mkdir("./$OutputDir");

#Read in Del Rates
my $DomArch2DelsHash = {};

open DELS, "<$DelRatesFile" or croak "No Deletions File:\n\n".$!.$?;

while (my $line = <DELS>){
	
	chomp($line);
	my ($DomArch,$DelRate) = split(/:/,$line);
	$DomArch2DelsHash->{$DomArch}=$DelRate;
}
close DELS;

#Read in Scores
my $DomArch2ScoresHash = {};
open SCORES, "<$rawdatafile" or croak "No scores file".$!.$?;

while (my $line = <SCORES>){
	
	chomp($line);
	my ($DomArch,$Score) = split(/:/,$line);
	$DomArch2ScoresHash->{$DomArch}=$Score;
}
close SCORES;

my @DomArchs = keys(%$DomArch2ScoresHash);
my @StudiedDomArchs = map{$DomArchs[POSIX::ceil(rand(scalar(@DomArchs)))]}(1 .. $NumberOfExamples);
my $StudyDAHash = {};
map{$StudyDAHash->{$_}=undef}@StudiedDomArchs;

my $RawSimsHash = {};
map{$RawSimsHash->{$_}={}}@DomArchs;

open RAWSIMS, "<$SimulationsDataFile" or croak "No Simulations file".$!.$?;

while (my $line = <RAWSIMS>){
	
	chomp($line);
	my ($GenomesCounts,$DomArch,$SimulatedGenomes) = split(/:/,$line);
	
	next unless(exists($RawSimsHash->{$DomArch}));
	
	my $InfoHash = $RawSimsHash->{$DomArch};
	
	$InfoHash->{'Simulations'}=[] unless(exists($InfoHash->{'Simulations'}));
	push(@{$InfoHash->{'Simulations'}},split(/,/,$SimulatedGenomes));
	
	my ($cladesize,$observedgenomes) = split(',',$GenomesCounts);
	$InfoHash->{'ObservedNumberOfGenomes'}=[$observedgenomes,$cladesize];
}
close RAWSIMS;

open SUMMARY, ">./$SummaryDir/Summary.txt" or die $!."\t".$?;
print SUMMARY "DomArch\tCladeSize\tCladeOccupance\tCladeFraction\tDelRate\tPQ\n";

my $PredictivePosteriorHash = {};
@{$PredictivePosteriorHash}{(0.01,0.02,0.05,0.1,0.25,0.5)}=([],[],[],[],[],[]);

my $OneTailPredictivePosteriorHash = {};
@{$OneTailPredictivePosteriorHash}{(0.01,0.02,0.05,0.1,0.25,0.5)}=([],[],[],[],[],[]);

open CLADESIZE, ">.CladeSizes.dat" or die $!.$?;

# for number_of_examples random selection of results
foreach my $DomArch (@DomArchs){

	my $DelRate;
	$DelRate = 'unkown' unless(exists($DomArch2DelsHash->{$DomArch}));
	$DelRate = $DomArch2DelsHash->{$DomArch} if(exists($DomArch2DelsHash->{$DomArch}));
	
	my $Score = $DomArch2ScoresHash->{$DomArch};
	
	foreach my $PosteriorPreditiveQuartile (keys(%$PredictivePosteriorHash)){
		
		if($Score > 0.5){
			
			push(@{$PredictivePosteriorHash->{$PosteriorPreditiveQuartile}},$DomArch) if(abs($Score - 1) <= $PosteriorPreditiveQuartile/2);
			push(@{$OneTailPredictivePosteriorHash->{$PosteriorPreditiveQuartile}},$DomArch) if(abs($Score - 1) <= $PosteriorPreditiveQuartile/2);
		}else{
			
			push(@{$PredictivePosteriorHash->{$PosteriorPreditiveQuartile}},$DomArch) if(($Score) <= $PosteriorPreditiveQuartile/2);
		}
		#If the score liesoutside the major density of the posterior at some range, push it onto the list
		
	}
		
	my ($NumberOfGenomes,$cladesize) = @{$RawSimsHash->{$DomArch}{'ObservedNumberOfGenomes'}};
	
	
	if (exists($StudyDAHash->{$DomArch})){
	#Create css overlay
	`prepareTraitTreeOverlay.pl -da $DomArch -ss 1 -t $treefile`;
	
	#Pump out the tree into output dir, labelling the file with dom arch, score and del rate
	my $FileOut = "./$OutputDir/".$DomArch.".DelRate.".$DelRate.".Score.".$Score.".svg";
	
	`nw_display -sr -S -w $Width -c treedisplayoptions.css ./$treefile > $FileOut`;
	
	#Create css overlay
	`prepareTraitTreeOverlay.pl -da $DomArch -ss 0 -t $treefile`;
	
	#Pump out the tree into output dir, labelling the file with dom arch, score and del rate
	$FileOut = "./$OutputDir/".$DomArch.".DelRate.".$DelRate.".Score.".$Score.".WithSupras.svg";
	
	`nw_display -sr -S -w $Width -c treedisplayoptions.css ./$treefile > $FileOut`;
	
	open FH, ">.TempHistFile.dat";
	print FH join("\n",@{$RawSimsHash->{$DomArch}{'Simulations'}});
	close FH;
	
	`Hist.py -f './.TempHistFile.dat' -o ./$OutputDir/$DomArch.HistPlacement.png -t '$DomArch \nSimulations' -x 'Number Of Genomes' -y 'Frequency' -l 'Score:$Score\nDel Rate:$DelRate\nClade Size:$cladesize Observed:$NumberOfGenomes' --vline $NumberOfGenomes --column 0`;
	
	}

	print CLADESIZE $DomArch.":".$cladesize."\n";
	my $cladeFraction = $NumberOfGenomes/$cladesize;
	print SUMMARY $DomArch."\t".$cladesize."\t".$NumberOfGenomes."\t".$cladeFraction."\t".$DelRate."\t".$Score."\n";
	
}

close CLADESIZE;
close SUMMARY;

open PREDITIVESUMMARY, ">./$SummaryDir/QuartileDetails.dat" or die $!."\t".$?;
print PREDITIVESUMMARY "Quantile\tTwoTailedDensity\tTopHalf(OneTailed)Density%\tInverseTwoTail\n";	

foreach my $PosteriorQuartile (sort(keys(%$PredictivePosteriorHash))){
	
	my $ListOfDAs = $PredictivePosteriorHash->{$PosteriorQuartile};
	
	open PREDITIVEPOSTERIOR, ">./$SummaryDir/PerArchAcceptReject".$PosteriorQuartile.".dat" or die $!."\t".$?;
	print PREDITIVEPOSTERIOR join("\n",@$ListOfDAs);
	close PREDITIVEPOSTERIOR;
	
	my $oneTailedDAs = $OneTailPredictivePosteriorHash->{$PosteriorQuartile};
	
	my $TwoTailedpercent = 100*scalar(@$ListOfDAs)/scalar(@DomArchs);
	my $OneTailedpercent = 100*scalar(@$oneTailedDAs)/scalar(@DomArchs);
	my $quantile = 100*$PosteriorQuartile;
	my $inversedensity = 100 - $TwoTailedpercent;
	print PREDITIVESUMMARY $quantile."%\t".$TwoTailedpercent."%\t".$OneTailedpercent."\t".$inversedensity."%\n";	
}

close PREDITIVESUMMARY;


unless(-e "./CladeSizesHist.png"){
	`Hist.py -f './.CladeSizes.dat' -o ./CladeSizesHist.png -t 'Histogram Of Cladesizes' -x 'Clade Size' -y 'Frequency' -l 'Number of Dom Archs' --column 1`;
}

my @Scores = values(%$DomArch2ScoresHash);

my $NumberLHSscores = grep{$_ < 0.5}@Scores;
my $NumberRHSscores = grep{$_ > 0.5}@Scores;

my $Asymmetry = 100*($NumberRHSscores - $NumberLHSscores)/scalar(@Scores);

$NumberLHSscores = grep{$_ < 0.5 && $_ > 0.1}@Scores;
$NumberRHSscores = grep{$_ > 0.5 && $_ < 0.9}@Scores;

my $EightyPercentAsymmetry =  100*($NumberRHSscores - $NumberLHSscores)/($NumberLHSscores + $NumberRHSscores);

open ASYM, ">Asymmetry.txt" or die $!.$?;
print ASYM "Total Asymmetry: ".$Asymmetry."%\n";
print ASYM "Mid-Eighty Percent".$EightyPercentAsymmetry."%";
close ASYM;

#Calculate overall left-right symmetry
#Calculate left-right symmetry in 0.1-0.9

__END__