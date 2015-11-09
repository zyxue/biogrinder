#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;


my ($factory, $read, $nof_reads, %got_amplicons, %expected_amplicons);


# Template with several matching amplicons and forward primer only

ok $factory = Grinder->new(
   -reference_file  => data('multiple_amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')            ,
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 100                                  ,
   -total_reads     => 100                                  ,
), 'Forward primer only';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   $got_amplicons{$read->seq} = undef;
   ok_read_forward_only($read, 1, $nof_reads);
};
is $nof_reads, 100;

%expected_amplicons = (
   'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGT'      => undef,
   'AAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGTccccc' => undef,
   'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaggggggggaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGTggggg' => undef,
);

is_deeply( \%got_amplicons, \%expected_amplicons );
undef %got_amplicons;


# Template with several matching amplicons and forward and reverse primers

ok $factory = Grinder->new(
   -reference_file  => data('multiple_amplicon_database.fa'),
   -forward_reverse => data('forward_reverse_primers.fa')   ,
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 100                                  ,
   -total_reads     => 100                                  ,
), 'Forward and reverse primers';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   $got_amplicons{$read->seq} = undef;
   ok_read_forward_reverse($read, 1, $nof_reads);
};
is $nof_reads, 100;

%expected_amplicons = (
   'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT' => undef,
   'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaggggggggaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT' => undef,
   'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGT' => undef,
   'AAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT' => undef,
);
is_deeply( \%got_amplicons, \%expected_amplicons );
undef %got_amplicons;


# Template with several nested amplicons and forward and reverse primers

ok $factory = Grinder->new(
   -reference_file  => data('nested_amplicon_database.fa'),
   -forward_reverse => data('forward_reverse_primers.fa')   ,
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 100                                  ,
   -total_reads     => 100                                  ,
), 'Forward and reverse primers, nested amplicons';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   $got_amplicons{$read->seq} = undef;
   ok_read_forward_reverse($read, 1, $nof_reads);
};
is $nof_reads, 100;

%expected_amplicons = (
   'AAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT' => undef,
   'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGT' => undef,
   'AAACTTAAAGGAATTGACGGggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggGTACACACCGCCCGT' => undef,
);

is_deeply( \%got_amplicons, \%expected_amplicons );
undef %got_amplicons;


done_testing();




sub ok_read_forward_reverse {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   like $read->reference->id, qr/^seq\d+$/;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $readseq = $read->seq;
   is $read->id, $nof_reads;
   is $read->length, 95;
}

sub ok_read_forward_only {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   like $read->reference->id, qr/^seq\d+$/;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $readseq = $read->seq;
   is $read->id, $nof_reads;
   my $readlength = $read->length;
   ok ( ($readlength == 95) or ($readlength == 100) );
}
