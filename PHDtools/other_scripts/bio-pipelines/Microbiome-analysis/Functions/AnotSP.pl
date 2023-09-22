#!/usr/bin/perl -w

use warnings;

use strict;

my %TaxoTree;

my %TaxoLeaf;

open(FILE, "<", $ARGV[0]);

for my $i (<FILE>){

	chomp($i);
	
	my @line = split/\t/,$i;
	
	my @subline = split/\-/,$line[$#line];
	
	if(grep { $_ eq "2:superkingdom" } @subline){
	
		if($subline[0] eq "species"){
		
			$TaxoTree{$line[0]} = $line[1];
		
		}
		else{
		
			for my $j (@subline[1..$#subline]){
			
				if($j =~ m/\:species/){
				
					my @target = split/\:/,$j;
					
					$TaxoLeaf{$line[1]} = $target[0];
					
					last;
				
				}
			
			}
		
		}
	
	}

}

close FILE;

my @out;

push @out, join("\t", ("Accession ID", "Species name", "Query length", "Aligned length"));

push @out, "\n";

open(FILE, "<", $ARGV[1]);

while(<FILE>){

	chomp;
	
	if($_ !~ m/Accession ID/){
	
		my @line = split/\t/,$_;
		
		if($TaxoLeaf{$line[1]}){
		
			$line[1] = $TaxoTree{$TaxoLeaf{$line[1]}};
			
			push @out, join("\t", @line);
			
			push @out, "\n";
		
		}
		else{
		
			push @out, join("\t", @line);
			
			push @out, "\n";
		
		}
	
	}

}

close FILE;

my $name = $ARGV[1];

$name =~ s/_taxid2name\.map//;

open(OUT, ">", $name."_species.ant");

print OUT @out;

close OUT;
