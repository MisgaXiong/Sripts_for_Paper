#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

my $total = `cat $ARGV[0] | wc -l`;

$total = $total -1;

my @out;

open(FILE, "<", $ARGV[0]);

my $i = 1;

my $filebatch = 1;

my $j = 1;

for my $row (<FILE>){

	my @line = split/\t/,$row;
	
	if($line[0] eq "Virus Strain Name"){
	
		next;
	
	}
	else{
	
		my @GISAID_info = grep{$_ =~ m/EPI_ISL_/} @line;
		
		my @GISAID_ACC;
		
		if($GISAID_info[0] =~ m/,/){
		
			my @sub_GISAID = split/,/,$GISAID_info[0];
			
			@GISAID_ACC = grep{$_ =~ m/EPI_ISL_/} @sub_GISAID;
		
		}
		else{
		
			push @GISAID_ACC, $GISAID_info[0];
		
		}
		
		if($i < 10000 and $j < $total){
		
			push @out, $GISAID_ACC[0]."\n";
			
			$j += 1;
		
			$i += 1;
		
		}
		elsif($i == 10000 and $j <= $total){
		
			push @out, $GISAID_ACC[0]."\n";
		
			open(OUT, ">", "GISAID_ACC_".$filebatch.".txt");
			
			print OUT @out;
			
			close OUT;
			
			undef @out;
			
			$filebatch += 1;
			
			$i = 1;
			
			$j += 1;
		
		}
		elsif($j == $total){
		
			push @out, $GISAID_ACC[0]."\n";
		
			open(OUT, ">", "GISAID_ACC_".$filebatch.".txt");
			
			print OUT @out;
			
			close OUT;
		
		}
	
	}

}
