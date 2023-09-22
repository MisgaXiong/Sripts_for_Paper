#!/usr/bin/perl -w

=head1 Description
	Count the number of reads in the fastq file
Usage
	perl fqstat.pl -FQ1 <FQ R1> -FQ2 <FQ R2>
Parameters
	-FQ1	[str]	Input the name of forward fastq file
	-FQ1	[str]	Input the name of reverse fastq file
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.13 16:25 0.0.1
=cut

use strict;
use warnings;
use Getopt::Long;

my ($R1, $R2, $help) = ("", "", "");

GetOptions(
	'FQ1=s' => \$R1,
	'FQ2=s' => \$R2,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$R1) or (!$R2) or ($help));

my ($total_reads, $R1_reads, $R2_reads) = (0, 0, 0);
my @R1_names = split/\./,$R1;
if($R1_names[$#R1_names] =~ m/gz/){
	open(FORWARD, "gzip -dc $R1|");
}
else{
	open(FORWARD, "<", $R1);
}
while(<FORWARD>){
	$R1_reads += 1;
}
close FORWARD;
$R1_reads /= 4;

my @R2_names = split/\./,$R2;
if($R2_names[$#R2_names] =~ m/gz/){
	open(FORWARD, "gzip -dc $R2|");
}
else{
	open(FORWARD, "<", $R2);
}
while(<FORWARD>){
	$R2_reads += 1;
}
close FORWARD;
$R2_reads /= 4;

$total_reads = $R1_reads + $R2_reads;

my @out;
push @out, join("\t", ("Number of total reads", $total_reads));
push @out, join("\t", ("R1 reads", $R1_reads));
push @out, join("\t", ("R2 reads", $R2_reads));
print join("\n", @out)."\n";
