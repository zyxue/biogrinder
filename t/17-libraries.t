#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $lib, $nof_libs, $nof_reads, $read);


# Multiple shotgun libraries

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -read_dist      => 48                         ,
   -num_libraries  => 4                          ,
   -total_reads    => 99                         ,
), 'Multiple shotgun libraries';

$nof_libs = 0;
while ( $lib = $factory->next_lib ) {
   $nof_libs++;
   $nof_reads = 0;
   while ( $read = $factory->next_read ) {
      $nof_reads++;
      ok_read($read, undef, $nof_reads, $nof_libs);
   };
   is $nof_reads, 99;
}
is $nof_libs, 4;


# Multiple mate pair libraries

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 99                         ,
   -read_dist      => 48                         ,
   -num_libraries  => 4                          ,
   -insert_dist    => 250                        ,
), 'Multiple mate pair libraries';

$nof_libs = 0;
while ( $lib = $factory->next_lib ) {
   $nof_libs++;
   $nof_reads = 0;
   while ( $read = $factory->next_read ) {
      $nof_reads++;
      ok_mate($read, undef, $nof_reads, $nof_libs);
   };
   is $nof_reads, 99;
}
is $nof_libs, 4;

done_testing();




sub ok_read {
   my ($read, $req_strand, $nof_reads, $nof_libs) = @_;
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
   }
   like $read->seq, qr/[$letters]+/;
   is $read->id, $nof_libs.'_'.$nof_reads;
   is $read->length, 48;
}


sub ok_mate {
   my ($read, $req_strand, $nof_reads, $nof_libs) = @_;
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
   my $id = $nof_libs.'_'.round($nof_reads/2).'/'.($nof_reads%2?1:2);
   is $read->id, $id;
   is $read->length, 48;
}

