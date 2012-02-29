#! /usr/bin/perl -w

=head1 NAME

ABC-HGT<.pl>

=head1 USAGE

 ABC-HGT.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script for use in conjunction with ABCtoolkit. ABCtoolkit takes care of all the fancy approximate bayesian computation gubbins and this script provides simulations based on a tree using the standard poisson model.

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
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
use DBI;

use Supfam::Utils;
use Supfam::hgt;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my $DeletionRate;
my $TreeFile;
my $model = 'poisson';
my $MRCA;


#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s"  => \$TreeFile,
           "del|d=f" => \$DeletionRate,
           "MRCA|m=s" => \$MRCA,
           "model|m:s" => \$model,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my ($root,$TreeCacheHash);

if($TreeFile){
	
	open TREE, "<$TreeFile" or die $!.$?;
	my $TreeString = <TREE>;
	close TREE;

	($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

}else{
	
	die "no tree file or SQL tree node or tree file provided as tree to calculate HGT upon\n";	
}


####

unless($model eq 'Julian' || $model eq 'poisson' || $model eq 'corrpoisson'){
			 	
			 die "Inappropriate model selected";
}

####



my ($selftest,$distribution,$RawResults,$DeletionsNumberDistribution);		


if($model eq 'Julian'){
									
			($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelJulian($MRCA,0,1,$DeletionRate,$TreeCacheHash);
																					
}elsif($model eq 'poisson'){
					
			($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelPoisson($MRCA,0,1,$DeletionRate,$TreeCacheHash);

}elsif($model eq 'corrpoisson'){
					
			($selftest,$distribution,$RawResults,$DeletionsNumberDistribution) = RandomModelCorrPoisson($MRCA,0,1,$DeletionRate,$TreeCacheHash);

}else{
			die "Inappropriate model selected";
}

my ($result) = keys(%$distribution); #Only a single result should be present

print $result."\n";


__END__