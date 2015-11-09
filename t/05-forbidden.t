#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read);


# Exclude forbidden characters

ok $factory = Grinder->new(
   -reference_file => data('dirty_database.fa'),
   -read_dist      => 80                       ,
   -random_seed    => 1233567890               ,
   -total_reads    => 10                       ,
), 'With dubious chars';

ok $read = $factory->next_read;
like $read->seq, qr/[N-]/i;


ok $factory = Grinder->new(
   -reference_file => data('dirty_database.fa'),
   -exclude_chars  => 'n-'                     , # case independent
   -read_dist      => 30                       ,
   -random_seed    => 1233567890               ,
   -total_reads    => 10                       ,
), 'Exclude chars';

while ( $read = $factory->next_read ) {
  unlike $read->seq, qr/[N-]/i;
}


ok $factory = Grinder->new(
   -reference_file => data('dirty_database.fa'),
   -exclude_chars  => 'N-'                     ,
   -read_dist      => 71                       ,
   -random_seed    => 1233567890               ,
   -total_reads    => 10                       ,
), 'Cannot generate read';

eval { $read = $factory->next_read };
like $@, qr/error/i;


# Delete forbidden characters

ok $factory = Grinder->new(
   -reference_file => data('dirty_database.fa'),
   -delete_chars   => 'N-'                     ,
   -read_dist      => 70                       ,
   -random_seed    => 1233567890               ,
   -total_reads    => 10                       ,
), 'Delete chars';

while ( $read = $factory->next_read ) {
  unlike $read->seq, qr/[N-]/i;
}

done_testing();
