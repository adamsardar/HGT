#!/usr/bin/env perl

=head1 NAME

PostHGTAnalysisI<.pl>

=head1 USAGE

 PostHGTAnalysis.pl [options -v,-d,-h] <ARGS> (-p1 --prior1 PRIORLIBRARY1_FROM_ABC_HGT | -p2 --prior2 PRIORLIBRARY2_FROM_ABC_HGT | -s1 --summary1 SUMMARY_FILE_FROM_POSTHGTCOMPARE | -s1 --summary2 SUMMARY_FILE_FROM_POSTHGTCOMPARE)
 
=head1 SYNOPSIS

This script will compare a bunch of sumamries of data for nice presentation in thesis/papers

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
use Supfam::Utils;

use Carp;
use Carp::Assert;
use Carp::Assert::More;
use Statistics::Basic qw(:all);
use Devel::Size qw(size total_size);
use List::Util;
use List::MoreUtils qw(minmax);

use Math::SimpleHisto::XS;
#For binning data


# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my $Outfile = "ParamComparisons.dat";

my $PriorLibrary1;
my $PriorLibrary2;
my $SummaryFile1;
my $SummaryFile2;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "outfile|o:s" => \$Outfile,
           "prior1|p1:s" => \$PriorLibrary1,
           "prior2|p2:s" => \$PriorLibrary2,
           "summary1|s1:s" => \$SummaryFile1,
           "summary2|s2:s" => \$SummaryFile2,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Read in Summary files

my $SummaryFiles = {};
my $Summaryfileindex = 0;

my $UniqDomArches = {};

foreach my $file ($SummaryFile1,$SummaryFile2){
	
	next unless($file);
	$Summaryfileindex++;
	
	$SummaryFiles->{$Summaryfileindex}={};
	
	open SUM,"<$file" or die $!."`t".$?."\n";

	my $HeadersLine = <SUM>;

	while (my $line = <SUM>){
		
		chomp($line);
		my @LineValues =split(/\s+/,$line);
		
		my $DomArch = shift(@LineValues);
		#Again, this is the domain architecture entry
		
		$SummaryFiles->{$Summaryfileindex}{$DomArch}={};
		$SummaryFiles->{$Summaryfileindex}{$DomArch}{'CladeSizeML'}= shift(@LineValues);
		$SummaryFiles->{$Summaryfileindex}{$DomArch}{'CladeOccupanceML'}= shift(@LineValues);
		$SummaryFiles->{$Summaryfileindex}{$DomArch}{'CladeFractionML'}= shift(@LineValues);
		$SummaryFiles->{$Summaryfileindex}{$DomArch}{'DelRate'}= shift(@LineValues);
		$SummaryFiles->{$Summaryfileindex}{$DomArch}{'PQML'}= shift(@LineValues);

		$UniqDomArches->{$DomArch}=undef;

	}
	#Populate a hash of structure $HashRef->{"SummaryFile".$Summaryfileindex}{DomainArchitecture} = {delrat => 'val', 'cladefraction = 'val'  ...}
	
	close SUM;

}


my $PriorFiles = {};
my $Priorfileindex = 0;

foreach my $PriorDir ($PriorLibrary1,$PriorLibrary2){
	
	next unless($PriorDir);
	my @PriorFiles = glob("$PriorLibrary1/*Posterior.dat");
	
	$Priorfileindex++;
	
	foreach my $PriorFile (@PriorFiles){
		
		$PriorFile =~ m/([\w,]+)-Posterior\.dat/g;
		
		my $DA = $1;

		open PRIOR,"<$PriorFile" or die $!."\t".$?."\n";
		
		my $PriorRates = []; 
		
		while (my $line = <PRIOR>){
			
			chomp($line);
			push(@$PriorRates,$line);
		}
		
		$PriorFiles->{$Priorfileindex}{$DA}{'PriorRates'}=$PriorRates;
		$PriorFiles->{$Priorfileindex}{$DA}{'MeanPriorRate'} = mean(@$PriorRates);
		$PriorFiles->{$Priorfileindex}{$DA}{'MedianPriorRate'} = median(@$PriorRates);
		
		if ($Priorfileindex > 1){
			
			if(exists($PriorFiles->{1}{$DA})){
					#Stop perl kicking is we haven't seen the DA before
					
				my ($minprior,$maxprior) = minmax((@{$PriorFiles->{1}{$DA}{'PriorRates'}},@{$PriorFiles->{2}{$DA}{'PriorRates'}}));
					
				 my $HistPrior1 = Math::SimpleHisto::XS->new(
		    			min => $minprior, max => $maxprior, nbins => 101
		 		 );
		 		 
		 		my $binlowerboundries = $HistPrior1->bin_lower_boundaries();
		 		my $Nbins = scalar(@$binlowerboundries);
		 		croak "Mismactch in binsizes\n" unless ($Nbins-1 == 100);
		 		
		 		 my $HistPrior2 = Math::SimpleHisto::XS->new(
		    			bins => $binlowerboundries
		 		 );
				
				$HistPrior1->fill($PriorFiles->{1}{$DA}{'PriorRates'});
				$HistPrior2->fill($PriorFiles->{2}{$DA}{'PriorRates'});
					
				my $BinMidPoints = $HistPrior1->bin_centers();
				
				my $KL12 = 0;
				my $KL21 = 0;

				my $Prior1Size = scalar(@{$PriorFiles->{1}{$DA}{'PriorRates'}});
				my $Prior2Size = scalar(@{$PriorFiles->{2}{$DA}{'PriorRates'}});
				
				foreach my $bin (0 .. ($Nbins-2)){
					
					my $bin1size = $HistPrior1->binsize($bin);
					my $bin2size = $HistPrior2->binsize($bin);
					
					$KL12+=(log($bin1size*$Prior2Size/$bin2size*$Prior1Size)/log(2))*($bin1size/$Prior1Size);
					$KL21+=(log($bin2size*$Prior1Size/$bin1size*$Prior2Size)/log(2))*($bin2size/$Prior2Size);
				}
				
				$PriorFiles->{1}{$DA}{'KL12'}=$KL12;
				$PriorFiles->{2}{$DA}{'KL21'}=$KL21;
				#Calculate and then assaign the KL distance between the two distibutions
				
			}
			
			delete($PriorFiles->{1}{$DA}{'PriorRates'}) if(exists($PriorFiles->{1}{$DA}{'PriorRates'}));
			delete($PriorFiles->{2}{$DA}{'PriorRates'}) if(exists($PriorFiles->{2}{$DA}{'PriorRates'}));
			#No point keeping densities that we don't need
		}
		
		$UniqDomArches->{$DA}=undef;
		
		close PRIOR;
	}
}

if ($verbose){
	
	my $SumSize = total_size($SummaryFiles)/(1024*1024);
	my $PriorSize = total_size($PriorFiles)/(1024*1024);
	
	print STDERR "Size of SummaryFilesHash is ".$SumSize." Mbytes\n";
	print STDERR "Size of PriorLibraryHash is ".$PriorSize." Mbytes\n";
}

open PARAMCOMPARE,">$Outfile" or die $!."\t".$?;

my @DomArches2Study = keys(%$UniqDomArches);

print PARAMCOMPARE "DomArch	MeanP1	MedianP1	KLP1	MeanP2	MedianP2	KLP2	CladeSizeML1       CladeOccupanceML1  CladeFractionML1   DelRate1	PQML1	CladeSizeML2       CladeOccupanceML2  CladeFractionML2   DelRateML2	PQML2\n" ;

foreach my $DomArch (@DomArches2Study){

	print PARAMCOMPARE $DomArch."\t";
	
	for my $i (1,2) {

		if(exists($PriorFiles->{$i}{$DomArch})){
			
			print PARAMCOMPARE $PriorFiles->{$i}{$DomArch}{'MeanPriorRate'}."\t";
			print PARAMCOMPARE $PriorFiles->{$i}{$DomArch}{'MedianPriorRate'}."\t";
			
			if(exists($PriorFiles->{$i}{$DomArch}{'KL12'})){
				print PARAMCOMPARE $PriorFiles->{$i}{$DomArch}{'KL12'}."\t";
			}else{
				print PARAMCOMPARE "NA\t";
			}
			
		}else{
			
			print PARAMCOMPARE "NA\tNA\tNA\t";
		}
		
	}
	
	for my $j (1,2) {

		if(exists($SummaryFiles->{$j}{$DomArch})){
			
			print PARAMCOMPARE $SummaryFiles->{$j}{$DomArch}{'CladeSizeML'}."\t";
			print PARAMCOMPARE $SummaryFiles->{$j}{$DomArch}{'CladeOccupanceML'}."\t";
			print PARAMCOMPARE $SummaryFiles->{$j}{$DomArch}{'CladeFractionML'}."\t";
			print PARAMCOMPARE $SummaryFiles->{$j}{$DomArch}{'DelRate'}."\t";
			print PARAMCOMPARE $SummaryFiles->{$j}{$DomArch}{'PQML'};

		}else{
			
			print PARAMCOMPARE "NA\tNA\tNA\tNA\tNA";
		}
		
		if($j == 1){
			
			print PARAMCOMPARE "\t";
		}else{
			
			print PARAMCOMPARE "\n";
		}
		
	}
}


close PARAMCOMPARE;




__END__