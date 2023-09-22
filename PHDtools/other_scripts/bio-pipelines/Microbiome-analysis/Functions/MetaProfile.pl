#!/usr/bin/perl -w

use warnings;

use strict;

my %species;

open(FILE, "<", $ARGV[0]); #打开taxid2name或species.ant的映射文件;

while(<FILE>){

	chomp;
	
	if($_ !~ m/Accession ID/){
	
		my @line = split/\t/,$_;
		
		@{$species{$line[0]}{$line[1]}} = ($line[2], $line[3]); #species{contigs名称}{物种名称} = (contigs长度, 比对长度);
	
	}

}

close FILE;

my %species_abundance;

my %species_alignlen;

my %species_contiglen;

open(FILE, "<", $ARGV[1]); #打开RPM文件

my @unidentify;

while(<FILE>){

	chomp;
	
	if($_ !~ m/ContigsID/){
	
		my @line = split/\t/,$_;
		
		my @species_name = keys %{$species{$line[0]}};
		
		if(scalar(@species_name) >= 1){
		
			push @{$species_abundance{$species_name[0]}}, $line[1];
		
			push @{$species_contiglen{$species_name[0]}}, @{$species{$line[0]}{$species_name[0]}}[0];
		
			push @{$species_alignlen{$species_name[0]}}, @{$species{$line[0]}{$species_name[0]}}[1];
		
		}
		else{
		
			push @unidentify, $line[0]."\n";
		
		}
	
	}

}

my @out;

push @out, join("\t", ("Species", "Number of contigs", "Mean length of contigs", "Median length of contigs","Max length of contigs",
"Min length of contigs", "Mean length of align", "Median length of align", "Max length of align", "Min length of align",
"Mean abundance", "Median abundance", "Max abundance", "Min abundance"));

push @out, "\n";

for my $key (keys %species_contiglen){

	my @contigs_len = MFour(@{$species_contiglen{$key}});
	
	my @align_len = MFour(@{$species_alignlen{$key}});
	
	my @abundance = MFour(@{$species_abundance{$key}});
	
	push @out, join("\t", ($key, @contigs_len, @align_len[1..4], @abundance[1..4]));
	
	push @out, "\n";

}

my $name = $ARGV[1];

$name =~ s/_RPM\.tsv//;

open(OUT, ">", $name."_metagenomic.profile");

print OUT @out;

close OUT;

open(OUT, ">", $name."_contigs_unidentify.tsv");

print OUT @unidentify;

close OUT;

sub MFour {

	my @arr = @_;
	
	my $number = 0;
	
	my $min = $arr[0];
	
	my $max = 0;
	
	my $total = 0;
	
	for my $i (sort{$a<=>$b} @arr){
	
		$number += 1;
		
		$total += $i;
		
		if($min > $i){
		
			$min = $i;
		
		}
		
		if($max < $i){
		
			$max = $i;
		
		}
	
	}
	
	my $mean = $total/$number;
	
	my $median;
	
	my @arr_sort = sort{$a<=>$b} @arr;
	
	if($number % 2 == 0){
	
		$median = ($arr_sort[$number/2-1] + $arr_sort[$number/2])/2;
	
	}
	elsif($number % 2 != 0){
	
		$median = $arr_sort[($number-1)/2];
	
	}
	
	return($number, $mean, $median, $max, $min);

}