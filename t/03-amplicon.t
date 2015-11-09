#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);

# Forward primer only, forward sequencing

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -read_dist       => 48                          ,
   -total_reads     => 100                         ,
), 'Forward primer only, forward sequencing';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, 1, $nof_reads);
};
is $nof_reads, 100;



# Forward and reverse primers

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -read_dist       => 48                                ,
   -total_reads     => 100                               ,
), 'Forward then reverse primers, forward sequencing';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, 1, $nof_reads);
};
is $nof_reads, 100;


# Reverse primer only, reverse sequencing

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('reverse_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => -1                          ,
   -read_dist       => 48                          ,
   -total_reads     => 100                         ,
), 'Reverse primer only, reverse sequencing';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, -1, $nof_reads);
};
is $nof_reads, 100;


# Reverse and forward primers, reverse sequencing

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('reverse_forward_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => -1                                ,
   -read_dist       => 48                                ,
   -total_reads     => 100                               ,
), 'Reverse then forward primers, reverse sequencing';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, -1, $nof_reads);
};
is $nof_reads, 100;

done_testing();



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
   if ( $source =~ m/^seq1/ ) {
      $letters = 'a';
   } elsif ( $source =~ m/^seq2/ ) {
      $letters = 'c';
   } elsif ( $source =~ m/^seq3/ ) {
      $letters = 'g';
   } elsif ( $source =~ m/^seq4/ ) {
      $letters = 't';
   } elsif ( $source =~ m/^seq5/ ) {
      $letters = 'atg';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   like $read->seq, qr/[$letters]+/;
   is $read->id, $nof_reads;
   is $read->length, 48;
}

