#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;


use Grinder;
my ($factory, $nof_reads, $read);


# Initialization with short argument

ok $factory = Grinder->new(
   -rf => data('shotgun_database_extended.fa'),
   -tr => 10                         ,
), 'Shotgun & short arguments';

ok $factory->next_read;

# Long argument

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -read_dist      => 48                                  ,
   -total_reads    => 100                                 ,
), 'Long arguments';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, undef, $nof_reads);
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
   } elsif ( $source eq 'seq7' ) {
      $letters = 'a';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   like $read->seq, qr/[$letters]+/;
   is $read->id, $nof_reads;
   if ($source eq 'seq7') {
      is $read->length, 1;
   } else {
      is $read->length, 48;
   }
}

