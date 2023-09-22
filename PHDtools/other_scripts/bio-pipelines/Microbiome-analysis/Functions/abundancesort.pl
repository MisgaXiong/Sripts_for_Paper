#!/usr/bin/perl -w

use warnings;

use strict;

my @header;

my %abundance_table;

my %abundance_map;

my (@out, @vertout);

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;

	my @line = split/\t/,$_;
	
	if($line[0] eq "Species"){
	
		push @header, join("\t", (@line, "RPM", "RPKM"));
		
		push @header, "\n";
	
	}
	else{
	
		my $RPM = $line[1]*$line[10];
		
		my $RPKM = $RPM*(10**3)/($line[1]*$line[2]);
		
		push @{$abundance_table{$line[0]}}, join("\t", (@line, $RPM, $RPKM));
		
		push @{$abundance_table{$line[0]}}, "\n";
		
		$abundance_map{$line[0]} = $RPM;
	
	}

}

close FILE;

open(FILE, "<", $ARGV[1]);

my %vertplant;

for my $i (<FILE>){

	chomp($i);
	
	my @line = split/\t/,$i;
	
	my @subline = split/\-/,$line[$#line];
	
	foreach(@subline){
		
		if($_ eq "7742:clade" or $_ eq "35493:phylum"){
			
			$vertplant{$line[1]} = 1;
			
		}
		
	}

}

close FILE;

push @out, @header;

push @vertout, @header;

for my $species (sort{$abundance_map{$b}<=>$abundance_map{$a}} keys %abundance_map){

	if($vertplant{$species}){
	
		push @vertout, @{$abundance_table{$species}};
	
	}
	else{
	
		push @out, @{$abundance_table{$species}};
	
	}
	
}

my $name = $ARGV[0];

$name =~ s/\.profile//;

open(OUT, ">", $name."_sort.tsv");

print OUT @out;

close OUT;

$name =~ s/microbiome/vertplant/;

open(OUT, ">", $name."_sort.tsv");

print OUT @vertout;

close OUT;
