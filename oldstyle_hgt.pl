#!/usr/bin/perl -w
use lib "$ENV{HOME}/scripts/gough/lib";
use strict;
use DBI;
use JgoughSQL;
use Jgough;
use JgoughPhylo;

#This program is supposed to: compare orthology and paralogy using tree information
#Julian Gough 013.11.02

#SQL
my @temp=JgoughSQL->Password();
my $dsn= "DBI:mysql:superfamily:localhost";
my ($dbh,$sth);
my $user_name=$temp[0];
my $password=$temp[1];
#-----
my $prop;
#my @props=('0.33','0.40','0.50','0.60','0.66','0.75','0.80','0.83','1.00');
my @props=('0.33','0.50','0.75','0.83','1.00');
my $diffiterations=100;  #makes a big difference to speed
my $treefile;
my $usage;
my %genomes;
my %combs;
my $genquery;
my $comb;
my @noofgenomes;
my $i;
my $combos='';
my $suffix='';
my $set;
my $genome;
my $sf;
my $type;
my $html='y';
my $completes='n';
my %nodes;
my %groups=();
my %groups2=();
my @temp2;
my $flag;
my ($n,$x,$j,$above,$this);
my %grpszs;
my $high=0;
my $bin;
my @bins;
my @vals;
my @order;
my %binvals;
my $sum=0;
my %combs2;
my $i2;
my $random;
my @ranpool;
my $last;
my $rani;
my %bingrps;
my %bingrps2;
my @grps;
my @grps2;
my $sum2;
my @temp3;
my $ij;
my $totalseqs=0;
my $ff;
my $point;
my $node;
my $failsafe;
my $deleted;
my %dltd;
my $gnm;
my @gnms;
my %gotnodes;
my $diff;
my %diffdist;
my $nd;
my $tot=0;
my @difffromrand;
my @difffromrandline;
my $difff;
#sets--------------
my @gens;
#------------------

#ARGUEMENTS-------------------------------------------
$usage="orthoparalogy.pl <treefile> <output name> <only count combos y/n> <comb/supfam> <only count completes y/n>\n";
die "Usage: $usage" unless (@ARGV >= 2);
print STDERR "output.1 is the number of combinations occurring in X genomes, X=1..n\n";
print STDERR "output.2 is the number of combinations occurring in X groups, X=1..n\n";
print STDERR "output.3 is another program: convergent.pl\n";
print STDERR "output.4 is the same as (2), but random\n";
print STDERR "output.5 is a cumulative distribution plot of deletion rates in architectures\n";
print STDERR "output.6 is the difference from random\n";
print STDERR "output.7 is no. of distinct groups away from random\n";
print STDERR "output.8 is the ordered architectures and their occurrence\n";
$treefile=$ARGV[0];
$set=$ARGV[1];
$combos=$ARGV[2];
$type=$ARGV[3];
if (defined($ARGV[4])){
if ($ARGV[4] =~ /[Yy]/){
$completes='y';
$suffix=$suffix.'.p';
}
}
unless (defined($type)){
$type='comb';
}
unless (defined($combos)){
$combos='y';
$suffix=$suffix.'.c';
}
elsif ($combos =~ /[Yy]/){
$suffix=$suffix.'.c';
}
unless ($type =~ /[Cc]/){
$suffix=$suffix.'.s';
}
$dbh = DBI->connect ($dsn,$user_name,$password,{RaiseError => 1});
#-----------------------------------------------------

@temp=JgoughPhylo->ReadTree($treefile);
$i=$temp[0];
%nodes=%$i;

#Check-genomes-list
@temp=();
foreach $i (keys(%nodes)){
  if (length($i) == 2){
push @temp,$i;
  }
}
$genquery = join "' or genome='",@temp;
$genquery="(genome='$genquery')";
$sth = $dbh->prepare("SELECT include,password,genome FROM genome WHERE $genquery;");
$sth->execute();
while (@temp=$sth->fetchrow_array()){
$genomes{$temp[2]}=1;
  unless ($temp[0] eq 'y' and $temp[1] eq ''){
print STDERR "Genome: $temp[2] is in the tree but should not be!\n";die;
  }
}
#------------------

#Read-in-combs----------------------------------------
if ($type =~ /[Cc]/){
$sth = $dbh->prepare("SELECT genome,comb FROM comb WHERE $genquery;");
}
else{
$sth = $dbh->prepare("SELECT genome,sf FROM ass WHERE evalue <= 0.02 AND $genquery;");
}
$sth->execute();
while (@temp=$sth->fetchrow_array()){
unless ($temp[1] =~ /_gap_/ and $completes eq 'y'){
  if ($type =~ /[Ss]/ or $combos =~ /[Nn]/ or ($temp[1] =~ /\d+,\S+,/ or $temp[1] =~ /,\S+,\d+/ or $temp[1] =~ /\d+,\d+/)){
$totalseqs++;
$ff=0;
if (exists($combs{$temp[1]})){
  if (exists(${$combs{$temp[1]}}{$temp[0]})){
$ff=1;
  }
}
  if ($ff == 1){
${$combs{$temp[1]}}{$temp[0]}++;
  }
else{
${$combs{$temp[1]}}{$temp[0]}=1;
}
}
}
}
#-----------------------------------------------------

#loop#through#every#combination#######################################
$prop=1;
print STDERR scalar(keys(%combs))," combs to do\n";
foreach $comb (keys(%combs)){
$diffiterations=100;

@temp=keys(%{$combs{$comb}});
$sum++;
if ($sum/100 == int($sum/100)){
print STDERR "$sum combs of ",scalar(keys(%combs))," done now\n";
}
    if (scalar(@temp) > 1){

#Get-required-deletion-rates----
$this=$temp[0];
$failsafe=0;
$i=0;
until ($i == scalar(@temp)){
$failsafe++;if ($failsafe > 100){die "'until' loop not escaped\n";}
$i=0;
$j=0;
foreach $genome (split /:/,$this){
$j++;
if (exists($combs{$comb})){
  if (exists(${$combs{$comb}}{$genome})){
$i++;
  }
}
}
$last=$this;
$this=$nodes{$this};
}
$groups{$comb}=&Grouper($comb,$prop,\%combs);
#-------------------------------

#Average-the-diffs-------------
$difff=0;
for $ij (1 .. $diffiterations){
&MakeRand;
%{$combs2{$comb}}=();
foreach $gnm (split /:/,$last){
  unless (exists($dltd{$gnm})){
${$combs2{$comb}}{$gnm}=1;
  }
}
$groups2{$comb}=&Grouper($comb,$prop,\%combs2);
@temp=split /,/,$groups{$comb};
@temp2=split /,/,$groups2{$comb};
$diff=scalar(@temp2)-scalar(@temp);
$tot++;
if (exists($diffdist{$diff})){
$diffdist{$diff}++;
}
else{
$diffdist{$diff}=1;
}
$difff=$difff+$diff;
}
push @difffromrand,($difff/$diffiterations);
$i2=join ',',keys(%{$combs{$comb}});
push @difffromrandline,"$comb\t$i2\t$groups{$comb}\n";
#------------------------------

#Binning-----
$bin=(1-$i/$j);
if (exists($binvals{$bin})){
$binvals{$bin}++;
@temp = split /,/,$groups{$comb};
$bingrps{$bin}=$bingrps{$bin}.','.(scalar(@temp));
@temp = split /,/,$groups2{$comb};
$bingrps2{$bin}=$bingrps2{$bin}.','.(scalar(@temp));
}
else{
$binvals{$bin}=1;
@temp = split /,/,$groups{$comb};
$bingrps{$bin}=scalar(@temp);
@temp = split /,/,$groups2{$comb};
$bingrps2{$bin}=scalar(@temp);
}
#------------

}
 }

######################################################################

#print-out-deletion-rates--
@bins=();
@vals=();
open OUT,(">$set.5$suffix.dat");
  foreach $i (keys(%binvals)){
push @bins,$i;
push @vals,$binvals{$i};
  }
@order=Jgough->OrderArray2(@bins);
$j=0;
foreach $i (@order){
$j=$j+$vals[$i];
print OUT "$bins[$i] ",$j/$sum,"\n";
}
print OUT "\n";
close OUT;
#--------------------------

#print-out-grpszs
@bins=();
open OUT,(">$set.6$suffix.dat");
  foreach $i (keys(%binvals)){
push @bins,$i;
if (defined($bingrps{$i})){
push @grps,$bingrps{$i};
}
else{
push @grps,0;
}
if (defined($bingrps2{$i})){
push @grps2,$bingrps2{$i};
}
else{
push @grps2,0;
}
  }
@order=Jgough->OrderArray2(@bins);
foreach $i (@order){
@temp=split /,/,$grps[$i];
$sum=0;
foreach $j (@temp){
$sum=$sum+$j;
}
@temp=split /,/,$grps2[$i];
$sum2=0;
foreach $j (@temp){
$sum2=$sum2+$j;
}
print OUT "$bins[$i] ",($sum/$sum2),"\n";
}
close OUT;
#----------------

#print-out-difference-dist--
@bins=();
@vals=();
open OUT,(">$set.7$suffix.dat");
  foreach $diff (keys(%diffdist)){
push @bins,$diff;
push @vals,$diffdist{$diff};
  }
@order=Jgough->OrderArray2(@bins);
foreach $i (@order){
print OUT "$bins[$i] ",$vals[$i]/$tot,"\n";
}
print OUT "\n";
close OUT;
#---------------------------

#print-out-ordered-list-----
open OUT,(">$set.8$suffix.dat");
@order=Jgough->OrderArray2(@difffromrand);
foreach $i (@order){
print OUT $difffromrand[$i],"\t",$difffromrandline[$i];
}
print OUT "\n";
close OUT;
#---------------------------

#other-outputs-----------------
#&Group(2,\%combs,@props);
#&Group(4,\%combs2,@props);
#&NoOfGenomes;
#------------------------------










sub NoOfGenomes{
#Number-of-genomes-in-each-comb-----------------------
  foreach $comb (keys(%combs)){
    if ($html eq 'y'){
&Html;
    }
    if (defined($noofgenomes[scalar(keys(%{$combs{$comb}}))])){
$noofgenomes[scalar(keys(%{$combs{$comb}}))]++;
    }
else{
$noofgenomes[scalar(keys(%{$combs{$comb}}))]=1;
}
  }
open OUT,(">$set.1$suffix.dat");
for $i (0 .. scalar(@noofgenomes)-1){
  if (defined($noofgenomes[$i])){
print OUT  "$noofgenomes[$i]\n";
}
else{
print OUT "0\n";
}
}
close OUT;
#-----------------------------------------------------
return;
}

sub Html{
    if (scalar(keys(%{$combs{$comb}})) > 40){
print "<BR><B>Genome occurrance</B>: ",scalar(keys(%{$combs{$comb}})),"<BR>";
      foreach $genome (keys(%{$combs{$comb}})){
print ", <a href=cgi-bin/gen_list.cgi?genome=$genome>$genome</a>: ${$combs{$comb}}{$genome}";
      }
print "<BR>\nCombination:<BR>";
foreach $sf (split /,/,$comb){
  if ($sf eq '_gap_'){
print "$sf<BR>";
  }
else{
$sth = $dbh->prepare("SELECT description FROM des WHERE id=$sf;");
$sth->execute();
@temp=$sth->fetchrow_array();
print "::<a href=scop.mrc-lmb.cam.ac.uk/scop-1.63/search.cgi?sunid=$sf>$temp[0]</a><BR>";
}
}
}
return;
}



sub Group{
my %groups;
my @props=@_;
my $number=shift(@props);
my $pointer=shift(@props);
my %combshash=%$pointer;
  if (-e "$set.$number.dat"){
system ("rm $set.$number.dat");
  }
  foreach $prop (@props){
$high=0;
#GROUPS-----------------------------------------------
  foreach $comb (keys(%combshash)){
$groups{$comb}=&Grouper($comb,$prop,\%combshash);
}
unless (scalar(keys(%groups)) == scalar(keys(%combshash))){
die "Mismatch of combs: ",scalar(keys(%groups))," ",scalar(keys(%combshash)),"\n";
  }
#-----------------------------------------------------
&Groups($number,\%groups);
}
return(\%groups);
}

sub Groups{
my %grpszs;
my $num=$_[0];
my $pointer=$_[1];
my %groups=%$pointer;
#Print-out-groups-------------------------------------
  foreach $comb (keys(%groups)){
@temp = split /,/,$groups{$comb};
$i=@temp;
if ($i > $high){
$high=$i;
}
if (exists($grpszs{$i})){
$grpszs{$i}++;
}
else{
$grpszs{$i}=1;
}
  }
unless ($num == 0){
open OUT,(">>$set.$num.dat");
for $i (1 .. $high){
  if (exists($grpszs{$i})){
print OUT "$grpszs{$i}\n";
}
else{
print OUT "0\n";
}
}
print OUT "\n";
close OUT;
}
#-----------------------------------------------------
return;
}


sub Grouper{
my $comb=$_[0];
my $prop=$_[1];
my $point=$_[2];
my %combs=%$point;
my %done=();
my $i;
my ($x,$n,$this,$j);
my $groups;
my @temp;
my $flag;

    foreach $genome (keys(%{$combs{$comb}})){
      unless (exists($done{$genome})){

$this=$genome;
$flag=1;
until ($flag == 0){
$flag=0;
$above=$nodes{$this};
@temp=split /:/,$above;
$n=0;
$x=0;
foreach $i (@temp){
$n++;
if (exists($combs{$comb})){
  if (exists(${$combs{$comb}}{$i})){
$x++;
  }
}
}
      if ($x/$n >= $prop){
$flag=1;
$this=$above;
      }
unless ($flag == 1 and exists($nodes{$this})){
$j=0;
  foreach $i (split /:/,$this){
$j++;
$done{$i}=1;
  }
if (defined($groups)){
$groups="$groups,$j";
}
else{
$groups=$j;
}
$flag=0;
}
}

      }
    }
return($groups);
}

sub MakeRand{
#make-rand------------------------------------------
#make-pool-of-all-nodes
@ranpool=split /:/,$last;
my @temp=@ranpool;
%gotnodes=();
foreach $node (@temp){
$failsafe=0;
until ($node eq $last){
$failsafe++;if ($failsafe > 10000){die "'until' loop not escaped\n";}
$node=$nodes{$node};
unless (exists($gotnodes{$node})){
push @ranpool, $node;
$gotnodes{$node}=1;
}
}
}
#----------------------
#try-nodes-at-random---
$failsafe=0;
$deleted=0;
%dltd=();
until ($deleted == ($j-$i) or scalar(@ranpool) < 1){
$failsafe++;if ($failsafe > 10000){die "'until' loop not escaped\n";}
$random=int(rand(scalar(@ranpool)));
$node=$ranpool[$random];
splice(@ranpool,$random,1);
unless (exists($dltd{$node})){
@gnms=split /:/,$node;
if (scalar(@gnms) <= ($j-$i-$deleted)){
foreach $gnm (@gnms){
$dltd{$gnm}=1;
$deleted++;
$nd=$gnm;
until ($nd eq $last){
$nd=$nodes{$nd};
$dltd{$nd}=1;
}
}
}
}
}
#----------------------
#---------------------------------------------------
}
