#! /usr/bin/perl -w

=head1 NAME

HGTSingleDetailedSimulation<.pl>

=head1 USAGE

HGTSingleDetailedSimulation.pl -t| --tree <treefile> -dom|--domainarch <comb_id> 

=head1 SYNOPSIS

A script to take a dump of a hash from the Deletion Bias scripts and plot some diagnostic information

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

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


# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $dumpfile;

my $OutputDir = 'IndepthStudyOfDeletionBias';

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "dump|f=s" => \$dumpfile,
           "outdir|o:s" => \$OutputDir,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Create an output directory if it doesn't already exist
mkdir("./$OutputDir");

#Read in dumped data
my $DumpedHashRef = EasyUnDump($dumpfile);

#my $DumpedHashRef = {'SingleVals' =>  $SingleValues, 'ObservedDists' => $InterDeletionDistances,'scores' => $scores, 'selftest' => $selftest, 'Simulateddistances' => $SimulatedInterDeletionDistances};
			

my $ObservedDist = $DumpedHashRef->{'ObservedDists'};
my $SelfTestDists = $DumpedHashRef->{'selftest'};
my $Simulateddistances = $DumpedHashRef->{'Simulateddistances'};

open TEMPFILE, "> ./temp.dat" or die $!;
print TEMPFILE join("\n",@$Simulateddistances);
close TEMPFILE;

`Hist.py -f "./temp.dat" -o ./$OutputDir/SimulatedDists.png -t "Histogram of Raw Simulated InterDeletion Distances" -x "Phylogenetic Distance" -y "Frequency" -u 0` ;

open TEMPFILE, "> ./temp.dat" or die $!;
print TEMPFILE join("\n",@$ObservedDist);
close TEMPFILE;

`Hist.py -f "./temp.dat" -o ./$OutputDir/ObservedDists.png -t "Histogram of Raw Observed InterDeletion Distances" -x "Phylogenetic Distance" -y "Frequency" -u 0` ;


__END__