#!/usr/bin/perl -w

use warnings;

use strict;

my @without;

my @depbase = ("R", "Y", "M", "K", "S", "W", "B", "V", "D", "H", "N");

open(FILE, "<", $ARGV[0]);

for my $i (<FILE>){
	
	chomp($i);
	
	my @line = split/\t/,$i;
	
	if($line[0] eq "Site"){
	
		push @without, $i."\n";
	
	}
	else{
	
		my @mut = split/->/,$line[1];
		
		if(grep{ $_ eq $mut[1] } @depbase){
		
			next;
		
		}
		else{
		
			push @without, $i."\n";
		
		}
	
	}

}

close FILE;

print @without;
