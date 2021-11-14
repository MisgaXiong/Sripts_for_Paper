#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

my @out;

open(FILE, "<", $ARGV[0]);

for my $i (<FILE>){

	my @line = split/\t/,$i;
	
	if($line[0] eq "Virus Strain Name"){
	
		push @out, $i;
	
	}
	elsif($line[2] eq "GISAID" and $line[5] eq "Complete" and $line[7] eq "High"){
	
		push @out, $i;
	
	}
	elsif((grep{$_ =~ m/EPI_ISL_/} @line) and $line[5] eq "Complete" and $line[7] eq "High"){
	
		push @out, $i;
	
	}

}

close FILE;

print @out;
