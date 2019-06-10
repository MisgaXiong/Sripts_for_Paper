#!/usr/bin/perl

=head1 Name

	Step4_link_sort_allele.pl

=head1 Description

	After the Step3_sort_each_sample_allele.pl, the output requires to be linked.
    for i in `ls sort*.fasta`; do perl Step4_link_sort_allele.pl ${i} > link${i}; done

=head2 Usage

	perl Step4_link_sort_allele.pl sort.fasta > linked.fasta

=head2 Author

	Xiong Dongyan

=cut

use Bio::SeqIO;
my $input_sort_fa = $ARGV[0];
my $sort_fa_seq = Bio::SeqIO->new(
	-file => $input_sort_fa,
	-format => 'fasta',
);
my $link_seq;
while (my $seq = $sort_fa_seq->next_seq()) {
	my $id = $seq -> id;
	my @get_id = split/_/,$id;
	$seq_name = @get_id[0];
	my $seq_seq = $seq -> seq;
	$link_seq = $link_seq.$seq_seq;
};
my @total_seq = ">".$seq_name."\n".$link_seq."\n";
print @total_seq;
