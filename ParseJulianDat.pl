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

my $DataHash = {};
my $count = 0;

open INPUT, "<$inputfile";

while (my $line = <INPUT>){
	
	my $distribution = {};

	my @GenomeSimData = split(',',$line);
	map{my @splits = split(' ', $_); $distribution->{$splits[0]}= $splits[1];}@GenomeSimData;

	$DataHash->{$count++}=$distribution;
}

close INPUT;

my $Scores ={};


foreach my $DomArch (keys(%$DataHash)) {
		
		my $distribution = $DataHash{$DomArch};
			
		my $ResultsHash={};
		my $Average = 0;
		my $CumulativeCount = 0;
			
		for my $NumberOfGenomes (0 .. $CladeSize){
			
			if($NumberOfGenomes == $NoGenomes){
				
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

open OUT, "> $inputfile-Colsv.raw";

foreach my $FinalDomArch (keys(%$Scores)){

		my $score = $Scores->{$FinalDomArch};
		print OUT "$FinalDomArch:$score\n";
}

close OUT;
`Hist.py -f "$inputfile-Colsv.raw" -o $inputfile.png -t "Histogram of Scores" -x "Score" -y "Frequency"	` if(-e '~/bin/Hist.py');


__END__

