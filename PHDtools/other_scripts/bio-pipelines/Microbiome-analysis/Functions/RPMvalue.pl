#!/usr/bin/perl -w

use warnings;

use strict;

my @RPM;

my %sample_TPM;

my %total_reads;

my @statsfile = glob($ARGV[0]."/*stat");

my @group;

for my $i (0..$#statsfile){

	open(FILE, "<", $statsfile[$i]);
	
	my $name = $statsfile[$i];
	
	$name =~ s/\.stat//;
	
	push @group, $name;
	
	while(<FILE>){
	
		if($_ !~ /\*/){
		
			my @line = split/\t/,$_;
			
			$total_reads{$name} += $line[2];
			
			$sample_TPM{$line[0]}{$name} = $line[2];
		
		}
	
	}
	
	close FILE;

}

for my $contig (keys %sample_TPM){

	for my $sample (keys %{$sample_TPM{$contig}}){
	
		$sample_TPM{$contig}{$sample} = $sample_TPM{$contig}{$sample}*(10**6)/$total_reads{$sample};
	
	}

}

push @RPM, join("\t", ("ContigsID", @group));

push @RPM, "\n";

for my $contig (keys %sample_TPM){

	my @line_RPM;
	
	push @line_RPM, $contig;
	
	foreach (@group){
	
		push @line_RPM, $sample_TPM{$contig}{$_};
	
	}
	
	push @RPM, join("\t", @line_RPM);
	
	push @RPM, "\n";

}

print @RPM;
