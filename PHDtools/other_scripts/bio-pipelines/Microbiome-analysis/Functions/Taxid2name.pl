#!/usr/bin/perl -w

use warnings;

use strict;

my %taxid2name;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	if($_ =~ m/scientific name/){
	
		$_ =~ s/\|//g;
		
		my @line = split/\t/,$_;
		
		$taxid2name{$line[0]} = $line[2];
	
	}

}

close FILE;

open(FILE, "<", $ARGV[1]);

my @taxid2name;

push @taxid2name, join("\t", ("Accession ID", "Species name", "Query length", "Aligned length"));

push @taxid2name, "\n";

while(<FILE>){
	
	if($_ !~ m/Accession ID/){
	
		my @line = split/\t/,$_;
		
		$line[1] = $taxid2name{$line[1]};
		
		push @taxid2name, join("\t", @line);
	
	}

}

close FILE;

my $name = $ARGV[1];

$name =~ s/_acc2taxid\.map//;

open(OUT, ">", $name."_taxid2name.map");

print OUT @taxid2name;

close OUT;
