#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);

# Feed __DATA__ content to stdin
*STDIN = *DATA;

ok $factory = Grinder->new(
   -reference_file   => '-'  ,
   -total_reads      => 100  ,
   -read_dist        => 48   ,
), 'Input from stdin';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, undef, $nof_reads);
};
is $nof_reads, 100;


sub ok_read {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   my $source = $read->reference->id;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $letters;
   if ( $source eq 'seq1' ) {
      $letters = 'a';
   } elsif ( $source eq 'seq2' ) {
      $letters = 'c';
   } elsif ( $source eq 'seq3' ) {
      $letters = 'g';
   } elsif ( $source eq 'seq4' ) {
      $letters = 't';
   } elsif ( $source eq 'seq5' ) {
      $letters = 'atg';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   like $read->seq, qr/[$letters]+/;
   is $read->id, $nof_reads;
   is $read->length, 48;
}

done_testing();




__DATA__
>seq1 this is the first sequence
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>seq2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
>seq3
gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
>seq4
tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
>seq5 last sequence, last comment
aaaaaaaaaattttttttttttttttttttttttttttttttttttttttttttttttttttttttttttgggggggggg
