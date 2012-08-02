#! /usr/bin/env perl

=head1 NAME

CheckTreesI<.pl>

=head1 USAGE

 .pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A old script of Julians to check input trees. This might well get deleted form the repo.

=head1 AUTHOR

B<Joe Bloggs> - I<Joe.Bloggs@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2010 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/workspace/Oates/lib/";


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
use Supfam::SQLFunc;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $InputFile;


#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "file|f=s" => \$InputFile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Get other command line arguments that weren't optional flags.
my ($file)= @ARGV;

#Print out some help if it was asked for or if no arguments were given.
#pod2usage(-exitstatus => 0, -verbose => 2) if not $file or $help;

# Sub definitions
#----------------------------------------------------------------------------------------------------------------

sub ReadTree{
#ARGS: ReadTree($treefile)
#read-in-tree-----------------------------------------
my $leaf='';
my @tree;
my $next=-1;
my $gen;
my $i;
my $flag=1;
my (%nodeup,%nodedown);
my %distances;
my $treefile=$_[0];
my $node;
my ($one,$two,$three);
my @middles;
my $middle;
my $end;

open TREE,("$treefile");
while (<TREE>){
  if (/\S/){
$leaf=$leaf.$_;
chomp $leaf;
}
}
close TREE;
$leaf =~ s/\;//g;
@tree=split /,/,$leaf;
until ($flag == 0){
$flag=0;
$next=-1;
foreach $i (0 .. scalar(@tree)-1){
$leaf = $tree[$i];
unless ($leaf eq ':'){
unless ($next == -1){
  if ($leaf =~ /^([\w:]+):(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+)\)(\S*)$/){
$distances{$1}=$2;
$end="$1$3";
$node=$gen;
$one=$1;
    foreach $middle (@middles){
$middle =~ /^(\S+):(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+),(\d+)$/;
#Match either a number like 2.34995959 or 0.23E-94
$node=$node.':'.$1;
$tree[$3]=':';
    }
$tree[$i]=':';
$node="$node:$end";
$tree[$next]=$node;
$node =~ s/:(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+)\)//g;$node =~ s/:(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+)$//g;$node =~ s/\)//g;$node =~ s/\(//g;
$gen =~ s/\)//g;$gen =~ s/\(//g;
$one =~ s/\)//g;$one =~ s/\(//g;
$nodeup{$gen}=$node;if (exists($nodedown{$node})){unless ($nodedown{$node}=~/,/){$nodedown{$node}=$nodedown{$node}.",$gen";}}else{$nodedown{$node}=$gen;}
$nodeup{$one}=$node;if (exists($nodedown{$node})){unless ($nodedown{$node}=~/,/){$nodedown{$node}=$nodedown{$node}.",$one";}}else{$nodedown{$node}=$one;}
    foreach $middle (@middles){
$middle =~ /^([\w\:]+):(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+),(\d+)$/;$one=$1;$two=$2;
$distances{$one}=$two;if (exists($nodedown{$node})){unless ($nodedown{$node}=~/,/){$nodedown{$node}=$nodedown{$node}.",$one";}}else{$nodedown{$node}=$one;}
$nodeup{$one}=$node;
    }
$flag=1;
$next=-1;
  }
elsif ($leaf =~ /\)/ or $leaf =~ /\(/){
$next=-1;
  }
else{
push @middles,"$leaf,$i";
}
}
if($leaf =~ /^(\S*)\(([\w:]+):(-?\d+\.?\d*|-?\d+\.?\d*E{1}-{1}\d+)$/){
$distances{$2}=$3;
@middles=();
$next=$i;
$gen="$1$2";
}
}
}
}
$distances{$node}=0;
return(\%nodeup,\%distances,\%nodedown);
#-----------------------------------------------------
}

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Read in trees:

#Read-tree Julian
my @temp=&ReadTree($InputFile);

my %parent=%{$temp[0]};
my %distances=%{$temp[1]};
my %childs=%{$temp[2]};
#---------

#Read-tree BioPerl
my $input = new Bio::TreeIO(-file   => "$InputFile",
                            -format => "newick");
my $tree = $input->next_tree;
#Read in and initialise tree

my $root = $tree->get_root_node;
#---------

print "Children\n";
while (my ($key, $value) = each(%childs)){
	
	print "\t".$key." =>".$value."\n";
}

print "Parents\n";
while (my ($key, $value) = each(%parent)){
	
	print "\t".$key." =>".$value."\n";
}
print "Distances\n";
while (my ($key, $value) = each(%distances)){
	
	print "\t".$key." =>".$value."\n";
}

print "\n\n\n";

my @BioPNodes = ($root->get_all_Descendents);

print scalar(grep{$_->is_Leaf}@BioPNodes)."\n";

foreach my $node (@BioPNodes){
	
	my @ChildNodes = $node->each_Descendent;
	
	my $AllDescendatnsString = join(',',(map{$_->id}grep{$_->is_Leaf}($node->get_all_Descendents,$node)));
	my $CladeString = '';
	
	foreach  my $child (@ChildNodes){
		
		my @clade = map{$_->id}grep{$_->is_Leaf}($child->get_all_Descendents,$child);
		my $CladeName = join(',',@clade);

		
		$CladeString .= ":".$CladeName;
	}
	
	print "\t".$AllDescendatnsString." => ".$CladeString."\n" unless($node->is_Leaf);
}

#	my @clade = map{$_->id}$root->get_all_Descendents;
#	my $CladeName = join(',',@clade);
#	
#	print "\t".$CladeName." => 0\n";

#Attempt to populate the CONFIG variable from a local XML file.
#my $CONFIG = XMLin("config.xml", ForceArray => 0, KeyAttr => [ ])
#    or warn "Can't open local XML config file.";

#Database parameters get these from a local config.xml, or give defaults
#my $host=(not ref $CONFIG->{'database'}->{'host'})?$CONFIG->{'database'}->{'host'}:'hostname.co.uk';
#my $user=(not ref $CONFIG->{'database'}->{'user'})?$CONFIG->{'database'}->{'user'}:'db_username';
#my $password=(not ref $CONFIG->{'database'}->{'password'})?$CONFIG->{'database'}->{'password'}:'db_password';

__END__

