#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read);


# Prepend a single multiplex identifier (MID), ACGT, to shotgun reads

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -multiplex_ids  => data('mids.fa')            ,
   -num_libraries  => 1                          ,
   -read_dist      => 52                         ,
   -total_reads    => 9                          ,
), 'Single MID - shotgun';

while ( $read = $factory->next_read ) {
   is $read->length, 52;
   is substr($read->seq, 0, 4), 'ACGT';
};


# Prepend two multiplex identifiers, ACGT and AAAATTTT, to shotgun reads

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -multiplex_ids  => data('mids.fa')            ,
   -num_libraries  => 2                          ,
   -read_dist      => 52                         ,
   -total_reads    => 10                         ,
), 'Two MIDs - shotgun';

while ( $read = $factory->next_read ) {
   is $read->length, 52;
   like $read->id, qr/^1_/;
   is substr($read->seq, 0, 4), 'ACGT';
};

$factory->next_lib;

while ( $read = $factory->next_read ) {
   like $read->id, qr/^2_/;
   is $read->length, 52;
   is substr($read->seq, 0, 8), 'AAAATTTT';
};


# Prepend a single multiplex identifier to amplicon reads

ok $factory = Grinder->new(
   -reference_file  => data('single_amplicon_database.fa'),
   -multiplex_ids   => data('mids.fa')                    ,
   -num_libraries   => 1                                  ,
   -read_dist       => 70                                 ,
   -total_reads     => 10                                 ,
   -forward_reverse => data('forward_reverse_primers.fa') ,
   -unidirectional  => 1                                  ,
), 'Single MID - amplicon';

while ( $read = $factory->next_read ) {
   is $read->seq, 'ACGTAAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGC';
};


# Request too long of a read

ok $factory = Grinder->new(
   -reference_file  => data('single_amplicon_database.fa'),
   -multiplex_ids   => data('mids.fa')                    ,
   -num_libraries   => 1                                  ,
   -read_dist       => 80                                 ,
   -total_reads     => 10                                 ,
   -forward_reverse => data('forward_reverse_primers.fa') ,
   -unidirectional  => 1                                  ,
), 'Single MID - amplicon too long';

while ( $read = $factory->next_read ) {
   is $read->seq, 'ACGTAAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
};


# Prepend two multiplex identifiers to amplicon reads

ok $factory = Grinder->new(
   -reference_file  => data('single_amplicon_database.fa'),
   -multiplex_ids   => data('mids.fa')                    ,
   -num_libraries   => 2                                  ,
   -shared_perc     => 100                                ,
   -read_dist       => 74                                 ,
   -total_reads     => 10                                 ,
   -forward_reverse => data('forward_reverse_primers.fa') ,
   -unidirectional  => 1                                  ,
), 'Two MIDs - amplicon';

while ( $read = $factory->next_read ) {
   like $read->id, qr/^1_/;
   is $read->seq, 'ACGTAAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
};


done_testing();
