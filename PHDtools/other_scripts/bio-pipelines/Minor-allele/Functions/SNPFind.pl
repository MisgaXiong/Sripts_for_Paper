#!/usr/bin/perl -w

use warnings;

use strict;

use Bio::SeqIO;

my $Ref = Bio::SeqIO -> new(

	-file => $ARGV[0],
	
	-format => 'fasta',

);

my $mindep = 20;

my $SNP_Qual = 30;

my $MuAF = 0.5;

my $seq = $Ref -> next_seq;

my $seq_seq = $seq -> seq;

my @Refseq = split//,$seq_seq;

open(VCF, "<", $ARGV[1]);

my %VCF_SNP;

while(<VCF>){

	chomp;
	
	if($_ =~ m/#/){
	
		next;
	
	}
	else{
	
		my @line = split/\t/,$_;
		
		if(length($line[3]) == 1 and $line[5] >= $SNP_Qual){
		
			$VCF_SNP{$line[1]} = $line[3]."->".$line[4];
		
		}
	
	}

}

close VCF;

open(FILE, "<", $ARGV[2]);

my @out;

push @out, join("\t", ("Pos","Ref","Alt", "A", "T", "C", "G"));

push @out, "\n";

my @snps;

push @snps, join("\t", ("Pos","Ref->Alt", "A", "T", "C", "G"));

push @snps, "\n";

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] ne "Chr"){
	
		@line = @line[1..5];
	
		if(SUM(@line) >= 1){
	
			my $RefBase = $Refseq[($line[0]-1)];
		
			my ($AltBase, $AltFreq) = Major_Base(@line);
		
			if($RefBase ne $AltBase){
		
				push @out, join("\t", ($line[0], $RefBase, $AltBase, @line[1..$#line]));
			
				push @out, "\n";
		
			}
		
			if(SUM(@line) >= $mindep and $VCF_SNP{$line[0]} and $AltBase ne $RefBase and $AltFreq/SUM(@line) >= $MuAF){
		
				push @snps, join("\t", ($line[0], $RefBase."->".$AltBase, @line[1..$#line]));
			
				push @snps, "\n";
		
			}
	
		}
	
	}

}

close FILE;

my $name = $ARGV[2];

$name =~ s/-6base\.tsv//;

open(OUT, ">", $name."_poteintial_SNP.txt");

print OUT @out;

close OUT;

open(OUT, ">", $name."_SNP.txt");

print OUT @snps;

close OUT;

sub SUM {

	my @arr = @_;
	
	my $count = 0;
	
	for my $i (1..$#arr){
	
		$count += $arr[$i];
	
	}
	
	return($count);

}

sub Major_Base {

	my @arr = @_;
	
	my %num2base = (1=>"A", 2=>"T", 3=>"C", 4=>"G");
	
	my $major = "";
	
	my $max = 0;
	
	for my $i (1..$#arr){
	
		if($max <= $arr[$i]){
		
			$max = $arr[$i];
			
			$major = $num2base{$i};
		
		}
	
	}
	
	return($major, $max);

}
