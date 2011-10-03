#!/usr/bin/perl -w

use lib "$ENV{HOME}/workspace/Oates/lib/";

use strict;
use DBI;
use Supfam::SQLFunc;

no warnings "uninitialized";

#This program is supposed to: Calculate horizontal transfer
#Julian Gough 21.3.08

#SQL--
#my @temp=JgoughSQL->Password();
#my $dsn= "DBI:mysql:superfamily:supfam";
#my ($dbh,$sth);
#my $user_name=$temp[0];
#my $password=$temp[1];
#$dbh = DBI->connect ($dsn,$user_name,$password,{RaiseError => 1});

my $dbh = dbConnect;
my $sth;

#Variables------------------
my ( $usage, $treefile );
my (
    $excl,      $i,      $ii,       $jj,         $arch,
    $deletions, $clade,  $example,  $avegenomes, $genquery,
    $thisone,   $gennum, $selftest, $html
);
my (
    %exclude, %parent,  %genomes,      %distances, %childs,
    %stored,  %results, %distribution, @temp
);
my ( @archs, @gens );
my $j            = 0;
my $iarch        = 0;
my $iterations   = 100;
my $ratio        = 0;
my $falsenegrate = 0.45;
my $completes    = 'y';
my $average      = 0;
my $test         = 'n';

#---------------------------

#ARGUEMENTS-------------------------------------------
$usage = "hgt.pl <*treefile> <exclusions>\n";
die "Usage: $usage" unless ( @ARGV == 1 or @ARGV == 2 );
$treefile = $ARGV[0];
if   ( $treefile =~ /(\S+)\.\S+/ ) { $html = "$1.html"; }
else                               { $html = "$treefile.html"; }
open HTML, (">$html");
if   ( defined( $ARGV[1] ) ) { $excl = $ARGV[1]; }
else                         { $excl = 'none'; }

#-----------------------------------------------------

#exclusions (architectures)
unless ( $excl eq 'none' ) {
    open EX, ("$excl");
    while (<EX>) {
        if (/(\S+)/) {
            $exclude{$1} = 1;
        }
    }
    close EX;
}

#--------------------------

#Read-tree
@temp      = &ReadTree($treefile);
%parent    = %{$temp[0]};
%distances = %{$temp[1]};
%childs    = %{$temp[2]};

#---------

#Check-genomes-list
foreach $i ( keys(%parent) ) {
    if ( length($i) == 2 ) {
        push @gens, $i;
    }
}
$genquery = join "' or genome.genome='", @gens;
$genquery = "(genome.genome='$genquery')";
$sth =
  $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();
while ( @temp = $sth->fetchrow_array() ) {
    unless ( $temp[0] eq 'y' and $temp[1] eq '' ) {
    	
die "Genome: $temp[2] is in the tree but should not be!$temp[0] and $temp[1] \n";
    }
}

#------------------

# Check all the genomes in the input file (tree) are as anticipated

#Get-list-of-architectures---------------------------
$sth = $dbh->prepare(
"SELECT DISTINCT len_comb.comb FROM len_comb,genome WHERE len_comb.genome=genome.genome AND $genquery AND len_comb.comb != '_gap_';"
);
$sth->execute();
while ( @temp = $sth->fetchrow_array() ) {
    unless ( $completes eq 'y' and $temp[0] =~ /_gap_/ ) {
        unless ( exists( $exclude{$arch} ) ) {
            push(@archs, $temp[0]);
        }
    }
}

#----------------------------------------------------

#Main-loop-------------------------------------------
foreach $arch ( 0 .. scalar(@archs) - 1 ) {

    #$arch=$archs[(scalar(@archs)-1-$arch)];
    $arch = $archs[$arch];
    $iarch++;

    #get-genomes---
    %genomes = ();
    $example = 'none';
    $sth     = $dbh->prepare(
"SELECT DISTINCT len_comb.genome FROM len_comb,genome WHERE len_comb.genome=genome.genome AND $genquery AND len_comb.comb = '$arch';"
    );
    $sth->execute();
    while ( @temp = $sth->fetchrow_array() ) {
        if ( exists( $parent{ $temp[0] } ) ) {
            $genomes{ $temp[0] } = 1;
            $example = $temp[0];
        }
    }
    if ( $example eq 'none' ) {
        print STDERR "No genomes for this architecture: $arch\n";
        die;
    }

    #--------------
    #produce distributions
    $clade = &Clade( $example, \%parent, \%genomes );
    $deletions = &Deleted( $clade, \%parent, \%genomes, \%childs, \%distances );
    if ( $deletions > 0 ) {
        $thisone = ( length($clade) + 1 ) / 3;
        $thisone = $thisone . ":$deletions";
        unless ( exists( $stored{$thisone} ) ) {
            @temp =
              &RandomModel( $clade, $deletions, \%parent, \%childs, \%distances,
                $iterations, $falsenegrate );
            $selftest = pop(@temp);
            $stored{$thisone} = join ',', @temp;
        }

        #---------------------
        #work out addition to results for this architecture
        @temp         = split (/,/, $stored{$thisone});
        %distribution = @temp;
        $ii           = 0;
        
        
        if   ( $test eq 'y' ) { $gennum = $selftest; }
        else                  { $gennum = scalar( keys(%genomes) ); }
        for $i ( 0 .. ( length($clade) + 1 ) / 3 ) {
            if ( $i == $gennum ) {
                if ( exists( $distribution{$i} ) ) {
                    $jj = $ii + 1;
                    for ( 1 .. $distribution{$i} ) {
                      
                            $results{$jj} += 1 / $distribution{$i};
                            $average += $jj / $distribution{$i};
                            $jj++;
                    }
                }
                else {    
                        $results{$ii} += 0.5;
                        $results{ ( $ii + 1 ) } += 0.5;
                  	    $average += $ii;
                }
                
                last;
            }
            if ( exists( $distribution{$i} ) ) {
                $ii += $distribution{$i};
            }
        }
### This HIDEOUS piece of code is going to be used in calculating distance from the median. It sums the models seen up until the bin containing the genome of interest.

        #--------------------------------------------------
        #output
        $j++;
        print STDERR "Average: ",
          ( $average / $j - $iterations / 2 ) / $iterations,
          "      ... done $iarch of ", scalar(@archs), "\n";
        for $i ( 1 .. $iterations )the co {
            unless ( $i == 1 ) { print ","; }
            if ( exists( $results{$i} ) ) {
                print "$i $results{$i}";
            }
            else {
                print "$i 0";
            }
        }
        print "\n";
        print HTML
"<a href=http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/maketree.cgi?genomes=";
        print HTML join ',', keys(%genomes);
        print HTML ">$arch</a> Score: $ii<BR>\n";
$ii;
        #------
    }
}
close HTML;
dbDisconnect;

#----------------------------------------------------

sub Deviation {
    my $i    = $_[0];
    my $iter = $_[1];
    my %dist;

    return ( \%dist );
}

sub Clade {
    my $clade       = $_[0];
    my %parent      = %{ $_[1] };
    my %ArchGenomes = %{ $_[2] };

    my $flag = 1;
    my @CladeGenomes;
    my %children = ();

    until ( $flag == 0 ) {

        $flag = 0;
        @CladeGenomes = split /:/, $clade;

        foreach my $gen (@CladeGenomes) {

            $children{$gen} = 1;
        }

        foreach my $genome ( keys(%ArchGenomes) ) {

            unless ( exists( $children{$genome} ) ) {

                $flag = 1;
                die "Problem with the tree \n" if ( $parent{$clade} eq '' );
                $clade = $parent{$clade};
                last;    #Breaks foreach loop
            }
        }
    }
    return ($clade);

#Clade returns a ':' seperated list of all the genomes in the minimal clade of the tree that contains wihtin all of the genomes possesing the domain architecture of interest
}

sub Deleted {
    my $clade = $_[0];
    my ( $i, $flag );
    $i = $_[1];
    my %parent = %$i;
    $i = $_[2];
    my %genomes = %$i;
    $i = $_[3];
    my %childs = %$i;
    $i = $_[4];
    my %distances = %$i;
    my ( $dels, $time );
    $time = 0;
    $dels = 0;
    my @nodes;
    my @temp;

    if   ( length($clade) == 2 ) { @nodes = ($clade); }
    else                         { @nodes = split /,/, $childs{$clade}; }
    until ( scalar(@nodes) == 0 ) {
        $clade = pop(@nodes);

        #check to see if the clade is empty
        $flag = 0;
        @temp = split /:/, $clade;
        foreach $i (@temp) {
            if ( exists( $genomes{$i} ) ) {
                $flag = 1;
                last;
            }
        }

        #----------------------------------
        if ( $flag == 0 ) {
            $time += $distances{$clade} / 2;
            $dels++;
        }
        else {
            $time = $time + $distances{$clade};
            unless ( length($clade) ==  2 ) {
                @temp = split /,/, $childs{$clade};
                if ( scalar(@temp) != 2 ) {
                    print STDERR "More than one child to a node\n";
                    die;
                }
                foreach $i (@temp) {
                    push @nodes, $i;
                }
            }
        }
    }

    return ( ( $dels / $time ) );
}

sub RandomModel {
    my $clade        = $_[0];
    my $deletions    = $_[1];
    my $iterations   = $_[5];
    my $falsenegrate = $_[6];
    my $i;
    $i = $_[2];
    my %parent = %$i;
    $i = $_[3];
    my %childs = %$i;
    $i = $_[4];
    my %distances = %$i;
    my $cladetime = 0;
    my ( @temp, @nodes, @gens );
    my ( $delstodo, $iter, $gen, $extradels, $delpoint, $falses, $time,
        $selftest );
    my ( %genomes, %modelgens, %distribution );
    my $totalgens = 0;

    #work out time in clade
    if   ( length($clade) == 2 ) { @nodes = ($clade); }
    else                         { @nodes = split (/,/, $childs{$clade}); }
    until ( scalar(@nodes) == 0 ) {
        $clade = pop(@nodes);
        unless ( length($clade) == 2 ) {
            @temp = split (/,/, $childs{$clade});
            foreach $i (@temp) {
                push @nodes, $i;
            }
        }
        else {
            $genomes{$clade} = 1;
            push @gens, $clade;
        }
        $cladetime = $cladetime + $distances{$clade};
    }

    #----------------------
    $deletions = $deletions * $cladetime;

    #do iterations to get averages
    $selftest = int( rand($iterations) ) + 1;
    for $iter ( 1 .. $iterations ) {
        $falses    = 0;
        $delstodo  = int($deletions);
        $extradels = $deletions - int($deletions);
        if ( $extradels > rand(1) ) {
            $delstodo++;
        }

        #random model
        %modelgens = %genomes;
        until ( $delstodo == 0 ) {
            if ( scalar( keys(%modelgens) ) == 0 ) { last; }
            if ( $falsenegrate > rand(1) ) {
                $falses = int( rand( scalar(@gens) ) );
                $falses = $gens[$falses];
                delete( $modelgens{$falses} );
                $delstodo--;
                next;
            }
            $delpoint = rand($cladetime);
            $clade    = $_[0];
            $time     = 0;

            #tree
            @nodes = split (/,/, $childs{$clade});
            until ( scalar(@nodes) == 0 ) {
                $clade = pop(@nodes);
                unless ( length($clade) == 2 ) {
                    @temp = split (/,/, $childs{$clade});
                    foreach $i (@temp) {
                        push @nodes, $i;
                    }
                }
                $time = $time + $distances{$clade};
                if ( $time > $delpoint ) {
                    @temp = split (/:/, $clade);
                    foreach $gen (@temp) {
                        delete( $modelgens{$gen} );
                    }
                    last;
                }
            }

            #----
            $delstodo--;
        }
        $totalgens = scalar( keys(%modelgens) );
        if ( $iter == $selftest ) { $selftest = scalar( keys(%modelgens) ); }
        if ( exists( $distribution{$totalgens} ) ) {
            $distribution{$totalgens}++;
        }
        else { $distribution{$totalgens} = 1; }

        #------------
    }

    #-----------------------------

    @temp = %distribution;

    return ( @temp, $selftest );
}

sub ReadTree {

    #ARGS: ReadTree($treefile)
    #read-in-tree-----------------------------------------
    my $leaf = '';
    my @tree;
    my $next = -1;
    my $gen;
    my $i;
    my $flag = 1;
    my ( %nodeup, %nodedown );
    my %distances;
    my $treefile = $_[0];
    my $node;
    my ( $one, $two, $three );
    my @middles;
    my $middle;
    my $end;

    open TREE, ("$treefile");
    while (<TREE>) {
        if (/\S/) {
            $leaf = $leaf . $_;
            chomp $leaf;
        }
    }
    close TREE;
    $leaf =~ s/\;//g;
    @tree = split (/,/, $leaf);
    until ( $flag == 0 ) {
        $flag = 0;
        $next = -1;
        foreach $i ( 0 .. scalar(@tree) - 1 ) {
            $leaf = $tree[$i];
            unless ( $leaf eq ':' ) {
                unless ( $next == -1 ) {
                    if ( $leaf =~ /^([\w:]+):(-?\d+\.?\d*)\)(\S*)$/ ) {
                        $distances{$1} = $2;
                        $end           = "$1$3";
                        $node          = $gen;
                        $one           = $1;
                        foreach $middle (@middles) {
                            $middle =~ /^(\S+):(-?\d+\.?\d*),(\d+)$/;
                            $node = $node . ':' . $1;
                            $tree[$3] = ':';
                        }
                        $tree[$i]    = ':';
                        $node        = "$node:$end";
                        $tree[$next] = $node;
                        $node =~ s/:-?\d+\.?\d*\)//g;
                        $node =~ s/:-?\d+\.?\d*$//g;
                        $node =~ s/\)//g;
                        $node =~ s/\(//g;
                        $gen  =~ s/\)//g;
                        $gen  =~ s/\(//g;
                        $one  =~ s/\)//g;
                        $one  =~ s/\(//g;
                        $nodeup{$gen} = $node;

                        if ( exists( $nodedown{$node} ) ) {
                            unless ( $nodedown{$node} =~ /,/ ) {
                                $nodedown{$node} = $nodedown{$node} . ",$gen";
                            }
                        }
                        else { $nodedown{$node} = $gen; }
                        $nodeup{$one} = $node;
                        if ( exists( $nodedown{$node} ) ) {
                            unless ( $nodedown{$node} =~ /,/ ) {
                                $nodedown{$node} = $nodedown{$node} . ",$one";
                            }
                        }
                        else { $nodedown{$node} = $one; }
                        foreach $middle (@middles) {
                            $middle =~ /^([\w\:]+):(-?\d+\.?\d*),(\d+)$/;
                            $one             = $1;
                            $two             = $2;
                            $distances{$one} = $two;
                            if ( exists( $nodedown{$node} ) ) {
                                unless ( $nodedown{$node} =~ /,/ ) {
                                    $nodedown{$node} =
                                      $nodedown{$node} . ",$one";
                                }
                            }
                            else { $nodedown{$node} = $one; }
                            $nodeup{$one} = $node;
                        }
                        $flag = 1;
                        $next = -1;
                    }
                    elsif ( $leaf =~ /\)/ or $leaf =~ /\(/ ) {
                        $next = -1;
                    }
                    else {
                        push @middles, "$leaf,$i";
                    }
                }
                if ( $leaf =~ /^(\S*)\(([\w:]+):(-?\d+\.?\d*)$/ ) {
                    $distances{$2} = $3;
                    @middles       = ();
                    $next          = $i;
                    $gen           = "$1$2";
                }
            }
        }
    }
    $distances{$node} = 0;
    return ( \%nodeup, \%distances, \%nodedown );

    #-----------------------------------------------------
}

