#! /usr/bin/env perl

die "DO NOT use this script. It is old, not written by me and no longer maintained. It is kept here purely for legacy reasons!\n\n";

use lib "$ENV{HOME}/workspace/Oates/lib/";

use strict;
use warnings;
use DBI;
use Supfam::SQLFunc;

use Supfam::Utils;

use Time::HiRes;
use File::Temp;
#This program is supposed to: Calculate horizontal transfer
#Julian Gough 21.3.08

#SQL--
my $dbh = dbConnect;
my $sth;
#-----

no warnings 'uninitialized';

my $TEMPDISTPATH= File::Temp->newdir( DIR => './' , CLEANUP => 1);


#Variables------------------
my @temp;
my ($usage,$treefile);
my ($excl,$i,$ii,$jj,$arch,$deletions,$clade,$example,$avegenomes,$genquery,$thisone,$gennum,$selftest,$html);
my (%exclude,%parent,%genomes,%distances,%childs,%stored,%results,%distribution);
my (@archs,@gens);
my $j=0;my $iarch=0;
my $iterations=100;
my $ratio=0;
my $falsenegrate=0.00;
my $completes='y';
my $average=0;
my $test='n';
my $raw;
my $store=0; 
my $outputname;
my $check = 'n';
#---------------------------

#ARGUEMENTS-------------------------------------------
$usage="hgt.pl <*treefile> <exclusions>\n";
die "Usage: $usage" unless (@ARGV == 1 or @ARGV == 2);
$treefile=$ARGV[0];
if ($treefile =~ /(\S+)\.\S+/){$html="$1.html"; $raw = "$1-RAW.colsv"; $outputname=$1;}else{$html="$treefile.html";}
open HTML,(">$html");
if (defined($ARGV[1])){$excl = $ARGV[1];}else{$excl='none';}
	open DELS, ">./.JulDelRates.dat" or die $!;
	open RAW,  ">./$raw" or die $!;
#-----------------------------------------------------

#exclusions (architectures)
unless ($excl eq 'none'){
open EX,("$excl") or die $!;
while (<EX>){
if (/(\S+)/){
$exclude{$1}=1;
}
}
close EX;
}
#--------------------------

#Read-tree
@temp=&ReadTree($treefile);

%parent=%{$temp[0]};
%distances=%{$temp[1]};
%childs=%{$temp[2]};
#---------

#Check-genomes-list
foreach $i (keys(%parent)){
  if (length($i) == 2){
push @gens,$i;
  }
}
$genquery = join "' or genome.genome='",@gens;
$genquery="(genome.genome='$genquery')";

$sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();
while (@temp=$sth->fetchrow_array()){
  unless ($temp[0] eq 'y' and $temp[1] eq '' or $check eq 'n'){
print STDERR "Genome: $temp[2] is in the tree but should not be!\n";die;
  }
}
#------------------

$genquery = join "' or protein.genome='",@gens;
$genquery="(protein.genome='$genquery')";

my $lengenquery = join "' or len_comb.genome='",@gens;
$lengenquery="(len_comb.genome='$lengenquery')";

#Get-list-of-architectures---------------------------
$sth = $dbh->prepare("SELECT DISTINCT len_comb.comb FROM len_comb WHERE $lengenquery AND len_comb.comb != '_gap_';");
$sth->execute();
while (@temp=$sth->fetchrow_array()){
unless ($completes eq 'y' and $temp[0] =~ /_gap_/){
unless (exists($exclude{$temp[0]})){
	## ORIGINAL: unless (exists($exclude{$arch)){ - this can't be right, as $arch has yet to be set
push @archs,$temp[0];
}
}
}
#-----------------------------------------------
		
#Main-loop--------------------------------------
foreach $arch (0 .. scalar(@archs)-1){
	
	#my $tic = Time::HiRes::time; ########################
	
#$arch=$archs[(scalar(@archs)-1-$arch)];
$arch=$archs[$arch];

$iarch++;

#get-genomes---
%genomes=();$example='none';
$sth = $dbh->prepare("SELECT DISTINCT len_comb.genome FROM len_comb WHERE $lengenquery AND len_comb.comb = '$arch';");
$sth->execute();
while (@temp=$sth->fetchrow_array()){

if (exists($parent{$temp[0]})){
$genomes{$temp[0]}=1;
$example=$temp[0];
}
}

if ($example eq 'none'){
die "No genomes for this architecture: $arch\n";
}

#--------------
#produce distributions
$clade=&Clade($example,\%parent,\%genomes);
$deletions=&Deleted($clade,\%parent,\%genomes,\%childs,\%distances);


print DELS "$arch:$deletions\n";
#print "$arch:$deletions\n";

if ($deletions > 0){

#print "\n".$arch." => ";

$thisone=(length($clade)+1)/3;$thisone=$thisone.":$deletions";

my $NoGenomesInCalde;

  unless (exists($stored{$thisone}) && $store){
@temp=&RandomModel($clade,$deletions,\%parent,\%childs,\%distances,$iterations,$falsenegrate);


$selftest=pop(@temp);
$NoGenomesInCalde=pop(@temp);
my $DistHashPointer = pop(@temp);
$stored{$thisone}=$DistHashPointer;
}

%distribution=%{$stored{$thisone}};

$ii=0;
if ($test eq 'y'){$gennum=$selftest;}else{$gennum=scalar(keys(%genomes));}
for $i (0 .. (length($clade)+1)/3){
if ($i == $gennum){
  if (exists($distribution{$i})){
$jj=$ii+1;
for (1 .. $distribution{$i}){
  if (exists($results{$jj})){
$results{$jj}=$results{$jj}+1/$distribution{$i};
$average=$average+$jj/$distribution{$i};
}
else{
$results{$jj}=1/$distribution{$i};
$average=$average+$jj/$distribution{$i};
}
$jj++;
}
}
else{
  if (exists($results{$ii})){
$results{$ii}=$results{$ii}+0.5;
}
else{
$results{$ii}=0.5;
}
  if (exists($results{($ii+1)})){
$results{($ii+1)}=$results{($ii+1)}+0.5;
}
else{
$results{($ii+1)}=0.5;
}
$average=$average+$ii;
}
last;
}
  if (exists($distribution{$i})){
$ii=$ii+$distribution{$i};
}
}
#--------------------------------------------------
#output
$j++;print STDERR "Average: ",($average/$j-$iterations/2)/$iterations,"      ... done $iarch of ",scalar(@archs),"(".$iterations." iteration runs)\n";
#for $i (1 .. $iterations){
#  unless ($i == 1){print ",";}
#  if (exists($results{$i})){
#print "$i $results{$i}";
#  }
#else{
#print "$i 0";
#}
#}
#Commented out so as to simplify output
#print "\n";
print HTML "<a href=http://http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/maketree.cgi?genomes=";
print HTML join ',',keys(%genomes);
print HTML ">$arch</a> Score: $ii<BR>\n";

print RAW "$arch".':'."$ii\n"
#------
}

#my $toc = Time::HiRes::time; ########################
#print STDERR "Julian Loop takes:".($toc-$tic)."seconds\n"; ######################

}

close HTML;
close DELS;
close RAW;

`Hist.py -f "./$raw" -o $outputname.png -t "Histogram of Scores" -x "Score" -y "Frequency"	`;
`Hist.py -f "./.JulDelRates.dat" -o JulDelRates.png -t "Histogram of DeletionRates" -x "Deletion Rate" -y "Frequency"	-l Deletions`;
# Plot a couple of histograms for easy inspection of the data

open PLOT, ">./.JulPlotCommands.txt" or die $!;
print PLOT "Hist.py -f ./$raw -o $outputname.png -t 'Histogram of Scores' -x 'Score' -y 'Frequency'\n\n\n";
print PLOT "Hist.py -f './.JulDelRates.dat' -o 'JulDelRates.png' -t 'Histogram of DeletionRates' -x 'Deletion Rate' -y 'Frequency'	-l 'Deletions'";
close PLOT;

dbDisconnect($dbh) ;
#----------------------------------------------------

sub Deviation{
my $i=$_[0];
my $iter=$_[1];
my %dist;

return (\%dist);
}

sub Clade{
my $clade=$_[0];
my $flag=1;
my @gens;
my $i;
my %children=();
my %parent=%{$_[1]};
my %genomes=%{$_[2]};

until ($flag == 0){
$flag=0;
@gens=split /:/,$clade;
  foreach $i (@gens){
$children{$i}=1;
  }
foreach $i (keys(%genomes)){
  unless (exists($children{$i})){
$flag=1;

if ($parent{$clade} eq ''){die "Problem with the tree -$parent{$clade}- -$clade-\n";}
$clade=$parent{$clade};
last;
  }
}
}
return ($clade);
}

sub Deleted{
my $clade=$_[0];#MRCA
my ($i,$flag);
my %parent=%{$_[1]};#NodeUpHash
my %genomes=%{$_[2]};#Genomes with observed dom arch
my %childs=%{$_[3]};#NodeDownHash
my %distances=%{$_[4]};#Evolutionary distance from node to parent hash
my ($dels,$time);
$time=0;$dels=0;
my @nodes;
my @temp;

if (length($clade) == 2){@nodes=($clade);}else{@nodes=split /,/,$childs{$clade};}
until (scalar(@nodes) == 0){
$clade=pop(@nodes);
#check to see if the clade is empty
$flag=0;
@temp=split /:/,$clade;
foreach $i (@temp){
if (exists($genomes{$i})){
$flag=1;last;
}
}
#----------------------------------
if ($flag == 0){
$time=$time+$distances{$clade}/2;
$dels++;
}
else{
$time=$time+$distances{$clade};
unless (length($clade) == 2){
@temp=split /,/,$childs{$clade};
if (scalar(@temp) != 2){print STDERR "More than one child to a node\n";die;}
foreach $i (@temp){
push @nodes,$i;
}
}
}
}

return(($dels/$time));
}

sub RandomModel{
my $clade=$_[0];
my $CladeBack = $clade;
my $deletions=$_[1];
my $iterations=$_[5];
my $falsenegrate=$_[6];
my $i;$i=$_[2];my %parent=%$i;$i=$_[3];my %childs=%$i;$i=$_[4];my %distances=%$i;
my $cladetime=0;
my (@temp,@nodes,@gens);
my ($delstodo,$iter,$gen,$extradels,$delpoint,$time,$selftest);
my (%genomes,%distribution);
my $totalgens=0;

#work out time in clade
if (length($clade) == 2){@nodes=($clade);}else{@nodes=split /,/,$childs{$clade};}


until (scalar(@nodes) == 0){
$clade=pop(@nodes);
unless (length($clade) == 2){

map{push(@nodes,$_)}(split (',',$childs{$clade}));
}else{
$genomes{$clade}=1;
push @gens,$clade;
}
$cladetime=$cladetime+$distances{$clade};
}

#print "CladeTime =".$cladetime;

my $NoGenomes = keys(%genomes);


#----------------------
my $deletionsexpected=$deletions*$cladetime;

#do iterations to get averages
$selftest=int(rand($iterations))+1;
	


for $iter (1 .. $iterations){

	$delstodo=int($deletionsexpected);
	$delstodo++ if (($deletionsexpected-$delstodo) > rand(1));
	#print "Blah: $delstodo\n";
	
	
	#random model
	my %modelgens=%genomes;
	
		
my $count =0;
	while($delstodo){
		$count++;
				@nodes=();
				last unless (scalar(keys(%modelgens)));
		  	  
		 if ($falsenegrate > rand(1)){
			my $falses= $gens[int(rand(scalar(@gens)))];
			delete($modelgens{$falses});
			$delstodo--;
			next;
		  }
		  
		  
		$delpoint=rand($cladetime);
		$clade=$CladeBack;$time=0;

		@nodes=split(',',$childs{$clade});
		
		while ($clade=pop(@nodes)){
						
			map{push(@nodes,$_)}(split (',',$childs{$clade})) unless (length($clade) == 2);
			$time+=$distances{$clade};
			
			if ($time > $delpoint){
				
			map{delete($modelgens{$_})}(split (':',$clade));
			last;
			}
		}
		#----
		$delstodo--;
	}
	$totalgens=scalar(keys(%modelgens));
	if ($iter == $selftest){$selftest=scalar(keys(%modelgens));}
	$distribution{$totalgens}++;
		 
}

return (\%distribution,$NoGenomes,$selftest);
}

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


