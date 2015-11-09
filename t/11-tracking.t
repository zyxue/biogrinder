#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -total_reads    => 10                                  ,
   -unidirectional => 0                                   ,
   -desc_track     => 1                                   ,
), 'Bidirectional shotgun tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=.*position=(complement\()?\d+\.\.\d+(\))?/;
}


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -total_reads    => 10                                  ,
   -unidirectional => 1                                   ,
   -desc_track     => 1                                   ,
), 'Forward shotgun tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=.*position=\d+\.\.\d+/;
}


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -total_reads    => 10                                  ,
   -unidirectional => -1                                  ,
   -desc_track     => 1                                   ,
), 'Reverse shotgun tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=.*position=complement\(\d+\.\.\d+\)/;
}


ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -total_reads     => 10                          ,
   -desc_track      => 1                           ,
), 'Amplicon tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=\S+.*amplicon=\d+\.\.\d+.*position=.*/;
}


ok $factory = Grinder->new(
   -reference_file  => data('revcom_amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')          ,
   -length_bias     => 0                                  ,
   -unidirectional  => 1                                  ,
   -total_reads     => 10                                 ,
   -desc_track      => 1                                  ,
), 'Reverse-complemented amplicon tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=\S+.*amplicon=complement\(\d+\.\.\d+\).*position=.*/;
}


ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -total_reads     => 10                          ,
   -desc_track      => 1                           ,
   -chimera_perc    => 100                         ,
   -chimera_dist    => (1)                         ,
   -chimera_kmer    => 0                           ,
), 'Bimeric amplicon tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=\S+,\S+.*amplicon=\d+\.\.\d+,\d+\.\.\d+.*position=.*/;
}


ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -total_reads     => 10                          ,
   -desc_track      => 1                           ,
   -chimera_perc    => 100                         ,
   -chimera_dist    => (0, 1)                      ,
   -chimera_kmer    => 10                          ,
), 'Trimeric amplicon tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=\S+(,\S+){2}.*amplicon=\d+\.\.\d+(,\d+\.\.\d+){2}.*position=.*/;
}


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
   -desc_track     => 0                          ,
), 'No tracking';

ok $read = $factory->next_read;
while ($factory->next_read) {
   is $read->desc, undef;
}


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
), 'Tracking default';

ok $read = $factory->next_read;
while ($factory->next_read) {
   like $read->desc, qr/reference=.*position=.*(complement\()?\d+\.\.\d+(\))?/;
}

done_testing();
