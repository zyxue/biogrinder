#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read);


# No profile

ok $factory = Grinder->new(
    -reference_file => data('single_seq_database.fa'),
    -read_dist      =>  50                           ,
    -total_reads    =>  100                          ,
    -unidirectional =>  1                            ,
), 'No profile';

while ( $read = $factory->next_read ) {
   is $read->seq, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';   
};


# Grinder profile file that contains the same parameters as the previous test

ok $factory = Grinder->new(
   -profile_file => data('profile.txt'),
), 'Grinder profile';

while ( $read = $factory->next_read ) {
   is $read->seq, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';   
};


# A mix of profile and other command-line arguments

ok $factory = Grinder->new(
   -desc_track    => 0                  ,
   -num_libraries => 2                  ,
   -profile_file  => data('profile.txt'),
   -multiplex_ids => data('mids.fa')    ,
   -shared_perc   => 100                ,
), 'Mix of profile and manually-specified options';

while ( $read = $factory->next_read ) {
   is $read->seq, 'ACGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
   is $read->desc, undef;
};

done_testing();
