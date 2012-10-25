#!/usr/bin/env perl

=head1 NAME

CondenseOutput<.pl>

=head1 USAGE

CondenseOutput.pl [options -h] -rs --rawsims <raw_simulations_file.dat> -rd --rawdata <raw_data_file.colsv>

=head1 SYNOPSIS

Quick script to condense a bunch of HGT sim data files together. The defaults should be fine to rock with, but you can overule them

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2012 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Carp> More sophisticated error messaging.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Carp;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my @RAW_SIM_FILES = glob("RawSimulation*");
croak "Must supply a sim file!\n" unless(scalar(@RAW_SIM_FILES) >0);
my $RawSimFile = $RAW_SIM_FILES[0];

my @RAW_DAT_FILES = glob("*-RawData.colsv");
croak "Must supply a data file!\n" unless(scalar(@RAW_DAT_FILES) >0);
my $RawDatFile = $RAW_DAT_FILES[0];

#Set command line flags and parameters.
GetOptions(
           "help|h!" => \$help,
           "rawsims|rs:s" => \$RawSimFile,
           "rawdata|rd:s" => \$RawDatFile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

print STDERR "Using raw data file ".$RawDatFile."\n";
print STDERR "Using raw sim file ".$RawSimFile."\n";


my $RawSims = {};


open RAWDAT, "< $RawDatFile" or croak $RawDatFile.$?.$!; 


while(my $line = <RAWDAT>){
	
	chomp($line);
	my @data = split(':',$line);
	$RawSims->{$data[0]}=$data[1];
}
#$RawSims is now a hash of 'Domain arch' => Score

close RAWDAT;


open RAWSIMS, "< $RawSimFile" or croak $RawSimFile.$?.$!; 

print "DomainArchitecture\t CladeSize\tCladeOccupancy\tCladeFraction\tScore\n";

my $CladeFraction;

while(my $line = <RAWSIMS>){
	
	chomp($line);
	my ($cladeinfo,$domarch,undef) = split(':',$line);
	
	if(exists($RawSims->{$domarch})){
		
		my ($CladeSize,$CladeOccupancy) = split(',',$cladeinfo);
		
		$CladeFraction = $CladeOccupancy/$CladeSize;
		my $score = $RawSims->{$domarch};
		print $domarch."\t".$CladeSize."\t".$CladeOccupancy."\t".$CladeFraction."\t".$score."\n";
		
	}else{
	
		carp "No entry in domarch score hash for ".$domarch."!\n";	
	}

}

close RAWSIMS;


__END__
