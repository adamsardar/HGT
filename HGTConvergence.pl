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

open INPUT, "<$inputfile";

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

map{push(@{$Iterations},scalar(@{$RawResultsHash->{$_}}))}keys(%$RawResultsHash);
	
my $ITR = join("\n",@{$Iterations});

my $LR = [];
#This will be the left minus right density of the resulting histogram of scores

my $progress = Term::ProgressBar->new({name => 'Calculating Convergence Time',count => min(@{$Iterations}) , ETA => 'linear',});
  $progress->minor(0);
  my $next_update = 0;
  my $ProgressCount = 0;
#Initialise a progress bar

for my $index (1 .. min(@{$Iterations})){

$ProgressCount++;
	
	my $Scores = {};

	foreach my $DomArch (keys(%$RawResultsHash)) {
		
		my $CladeSize = $TreeDataHash->{$DomArch}{'CladeSize'};
		my $NoGenomes = $TreeDataHash->{$DomArch}{'Number_Of_Genomes_Observed'};
		
		my $distribution = {};
		my @Simulations = @{$RawResultsHash->{$DomArch}};
		my $c =0;
		
		map{$distribution->{$_}++}splice(@Simulations, int(rand(scalar(@Simulations))), 1) while scalar($c++) < $index;
		#Draw items from the list of simulation results, delete them from the original list and keep track of them in the hash %$distributions
		
		my $ResultsHash={};
		my $Average = 0;
		my $CumulativeCount = 0;
			
		for my $NumberOfGenomes (0 .. $CladeSize){
			
			if($NumberOfGenomes > $NoGenomes){
				
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

		$Scores->{$DomArch} = $CumulativeCount;
		
	}
	#This dull piece of code above works out the discrete distance from the median for each sampled distribution

	my $divider = (($index+1)/2);
	
	my $left = grep{$_ < $divider}(values(%$Scores));
	my $right = grep{$_ > $divider}(values(%$Scores));
	#Calculate the density the left and right of a histogram of scores 
	
	push(@$LR,($right-$left));
	#Store the result
	
	$next_update = $progress->update($ProgressCount) if $ProgressCount >= $next_update;

	if($index == min(@{$Iterations})){
		
		open OUT, "> .AggregatedDataColsv.raw";
		
		foreach my $FinalDomArch (keys(%$Scores)){

				my $score = $Scores->{$FinalDomArch};
				print OUT "$FinalDomArch:$score\n";
		}
		
		close OUT;

		open ASYM, '>./Asymmetry.dat' or die "$? $!";
		my $assymetry = 100*(($right-$left)/scalar(keys(%$Scores)));
		print ASYM $assymetry.'%';
		close ASYM;

		`Hist.py -f ".AggregatedDataColsv.raw" -o Aggregate.png -t "Histogram of Scores" -x "Score" -y "Frequency"	`;
	}
}
	
$progress->update($ProgressCount++);	
	
my $output = {};
my $count = 0;

map{$output->{$count++}=$_}@$LR;
	
open OUT, ">./.convergence.dat";
map{print OUT "$_".','."$output->{$_-1}"."\n"}(1 .. min(@{$Iterations}));
close OUT;

`EZLines.py -f ./.convergence.dat -o $outputfile -x 0 -y 1 --xlab 'Iterations' --ylab 'Asymetric Density' --title 'Convergence of Asymetric Density Estimate\n With Increasing Number of Iterations' -c b`;


__END__

