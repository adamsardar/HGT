#! /usr/bin/env perl

=head1 NAME

HGTAlgorithmBenchmarkI<.pl>

=head1 USAGE

 HGTAlgorithmBenchmark.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to compare implementations of HGT algorithms

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom";


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
use Supfam::TreeFuncsNonBP;
use Math::Random qw(random_uniform random_uniform_integer);

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'HGTBenchmark';
my $TreeFile;
my $FalseNegativeRate = 0.00;
my $completes = 'n'; #A flag to specify whether or not to include architectures including _gap_ assignmenets
my $Iterations = 1000;
my $store = 0;
my $check = 'y'; #Perform a sanity check on the tree? This should be 'y' unless under extreme circumstances

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$OutputFilename,
           "tree|i:s" => \$TreeFile,
           "no_iternations|itr:i" => \$Iterations,
           "false_negative_rate|fnr:f" => \$FalseNegativeRate,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Produce a tree hash, either from SQL or a provided treefile
my ($root,$TreeCacheHash);

if($TreeFile){
	
	open TREE, "<$TreeFile" or die $!.$?;
	my $TreeString = <TREE>;
	close TREE;

	($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

}else{
	
	die "no tree file or SQL tree node or tree file provided as tree to calculate HGT upon\n";	
}


my @FullRootClade = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};

my $CladeSize = scalar(@FullRootClade);

my ($RandomNumberOfCladeGenomes) =  random_uniform_integer(1,int($CladeSize/4),int(3*$CladeSize/4));
fisher_yates_shuffle(\@FullRootClade);

my @ObservedGenomes = @FullRootClade[(0 .. $RandomNumberOfCladeGenomes)];

my $HashOfGenomesObserved = {};
$HashOfGenomesObserved->{'ness'}={};
@{$HashOfGenomesObserved->{'ness'}}{@ObservedGenomes}=((1) x scalar(@ObservedGenomes));

my ($dels, $time) = DeletedJulian($root,0,0,$HashOfGenomesObserved,$TreeCacheHash,$root,'ness');

my $deletion_rate = $dels/$time;

my $tic = Time::HiRes::time;

my ($SelftestValueA,$distributionA,$RawResultsA,$DeletionsNumberDistributionA) = RandomModelCorrPoisson($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);

#my $distributionA = {};

#map{$distributionA->{$_}++}random_uniform_integer($Iterations,0,$CladeSize);

my $toc = Time::HiRes::time;

print "Done normal Poisson in ".($toc-$tic)."\n";

$tic = Time::HiRes::time;

#my $distributionB = {};
#map{$distributionB->{$_}++}random_uniform_integer($Iterations,0,$CladeSize);
my ($SelftestValueB,$distributionB,$RawResultsB,$DeletionsNumberDistributionB) = RandomModelCorrPoissonOptimised($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);
#my ($SelftestValueB,$distributionB,$RawResultsB,$DeletionsNumberDistributionB) = RandomModelCorrPoisson($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash);
$toc = Time::HiRes::time;

print "Done optimised Poisson ".($toc-$tic)."\n";

#Add pseudocount

#map{$distributionB->{$_} += 0.1; $distributionA->{$_} += 0.1;}(0 .. $CladeSize);

my $KLscore = KLdistance($distributionA,$distributionB);

print "KL score = ".$KLscore."\n";

__END__