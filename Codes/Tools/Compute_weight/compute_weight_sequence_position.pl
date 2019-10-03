#!/usr/bin/perl -w
#########
##
########

use strict;



if ($#ARGV < 0)
{
    helper();
}

sub helper 
{
	print STDERR "Title  : COMPUTE WEIGHTS FOR SEQUENCES AND POSITIONS FROM MULTIPLE ALIGNEMENT (from HENIKOFF & HENIKOFF 1994)\n";
    print STDERR "Author : JCG (jean-christophe.gelly\@univ-paris-diderot.fr)\n";
    print STDERR "Name   : $0\n";
	print STDERR "Usage  : $0 <file alignment fasta>\n";
	print STDERR "Example: perl $0 globin.fasta\n";
    print STDERR "Output : Most useful information is the \"Normalized_weight\" of the sequence. For an alignment the sum of all weight equal 1.0\n";

	exit 1;
}

#############

my $file=shift;

open(F,"$file") or die "Cannot open $file !\n";
my @tab_file=<F>;
close F;

############
print STDERR "#POSITION BASED SEQUENCE WEIGHT FROM HENIKOFF & HENIKOFF 1994\n";
#############
#Parse file
my $current_name;
my $current_seq;

my @tab_name;
my @tab_seq;

my $number_seq=-1;
for (my $i=0;$i <= $#tab_file ; $i++)
{
    my $line=$tab_file[$i];
    next if ($line =~/^\s+/);
    chomp $line;

    if ($line =~/^>(.+)$/)
    {
        $number_seq++;
        $tab_name[$number_seq]=$1;
        $tab_seq[$number_seq]="";
    }
    else
    {
        $tab_seq[$number_seq].=$line;
    }
}
$number_seq++;

if ($number_seq==0)
{
    print STDERR "No sequence in file!\n";
}

#############
#Parse position
my @tab_aa=split('',"ACDEFGHIKLMNPQRSTVWY-_");
#print "@tab_aa\n";
my @tab_seq_tab_position;

my $number_pos=-1;

for (my $i=0 ; $i < $number_seq ; $i++)
{
 #   print "$i ($number_seq): $tab_seq[$i]\n";
    my @tab_position=split('',$tab_seq[$i]);
    $number_pos=-1;
    for (my $j=0 ; $j <= $#tab_position ; $j++)
    {
        $number_pos++;
#       print "$j:$tab_position[$j] ";
        ${$tab_seq_tab_position[$i]}[$j]=$tab_position[$j];
    }
#  print "\n";
}
$number_pos++;
###############
# Compute value foreach position in sequence
my @tab_sequence_tab_position_computation;
for (my $j=0 ; $j < $number_pos ; $j++)
{
    my $number_different_aa=0;
    my %hash_compute_aa;
    #Initialisation
    foreach my $aa (@tab_aa)
    {
        $hash_compute_aa{"$aa"}=0;
 #       print "\$aa:$aa \$hash_compute_aa{\$aa}:$hash_compute_aa{$aa}\n";
    }
    #Comptabilisation
    for (my $i=0 ; $i < $number_seq ; $i++)
    {
  #      print "SEQ:$i position:$j ->\"${$tab_seq_tab_position[$i]}[$j]\"\n";
        my $aa=${$tab_seq_tab_position[$i]}[$j];
        if (exists $hash_compute_aa{$aa})
        {
            if ($hash_compute_aa{"$aa"}==0)
            {
                $number_different_aa++;
            }
        }
        else
        {

            print STDERR "WARNING ! UNKNOWN Amino Acid \"$aa\" replace by \"_\"\n";
            $aa="_";
        }
        $hash_compute_aa{$aa}++;
    }
    #Computation
    my $total_value_position=0;
    for (my $i=0 ; $i < $number_seq ; $i++)
    {
        my $aa=${$tab_seq_tab_position[$i]}[$j];
        my $value=1/($number_different_aa*$hash_compute_aa{$aa});
#       printf("$aa %2s / (%2s * %2s ) \n","1",$number_different_aa,$hash_compute_aa{$aa});
        ${$tab_sequence_tab_position_computation[$i]}[$j]=$value;
        $total_value_position+=$value;
    }
    #   print "TOTAL:$total_value_position\n";
}

##############
#COMPUTE VALUE FOR EACH SEQUENCE
my @tab_total_value_seq;
my @tab_total_normalized_value_seq;
for (my $i=0 ; $i < $number_seq ; $i++)
{
    my $total_value_seq=0;
    for (my $j=0 ; $j < $number_pos ; $j++)
    {
        $total_value_seq+=${$tab_sequence_tab_position_computation[$i]}[$j];
    }
    $tab_total_value_seq[$i]=$total_value_seq;
    $tab_total_normalized_value_seq[$i]=$total_value_seq/$number_pos;
    print ">$tab_name[$i]\n";
    printf("Normalized_weight: %5.3f\n",$tab_total_normalized_value_seq[$i]);
    for (my $j=0 ; $j < $number_pos ; $j++)
    {
        printf("%-5s ",${$tab_seq_tab_position[$i]}[$j]);
    }
    print "\n";
    for (my $j=0 ; $j < $number_pos ; $j++)
    {
        printf("%5.3f ",${$tab_sequence_tab_position_computation[$i]}[$j]);
    }
    print "\n";

}

################
#COMPUTE NORMALIZED VALUE FOR EACH SEQUENCE



