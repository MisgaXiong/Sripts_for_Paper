#!/usr/bin/perl -w

use warnings;
use strict;

my (%lin2ctg, %ctg2acc, %acc2taxid, %ctg2info);
open(FILE, "<", $ARGV[0]);
my $lin = 1;
while(<FILE>){
	chomp;
	my @line = split/\t/;
	my @subline = split/\|/,$line[0];
	$subline[1] =~ s/bp//;
	if($line[3]/$subline[1] >= 0.175){
		$ctg2acc{$line[0]}{$lin} = $line[1];
		@{$ctg2info{$line[0]}{$lin}} = ($line[3], $line[2], $line[$#line]);
		$acc2taxid{$line[1]} = 1;
		$lin2ctg{$lin} = $line[0];
		$lin += 1;
	}
}
close FILE;
open(FILE, "<", $ARGV[1]);
my @out;
push @out, join("\t", ("Accession ID", "Taxonomy ID", "Query length", "Aligned length", "Identity", "BitScore"))."\n";
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if(!$acc2taxid{$line[0]}){
		next;
	}
	else{
		$acc2taxid{$line[0]} = $line[2];
	}
}
for my $lin (sort{$a <=> $b} keys %lin2ctg){
	my $ctgs = $lin2ctg{$lin};
	my @subline = split/_/,$ctgs;
	push @out, join("\t", ($lin2ctg{$lin}, $acc2taxid{$ctg2acc{$lin2ctg{$lin}}{$lin}}, $subline[3], @{$ctg2info{$lin2ctg{$lin}}{$lin}}))."\n";
}
my $name = $ARGV[0];
$name =~ s/_blastn//;
open(OUT, ">", $name."_acc2taxid.map");
print OUT @out;
close OUT;
