#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $seed1, $seed2, $seed3, @dataset1, @dataset2);


# Seed the pseudo-random number generator

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -random_seed    => 1233567890                          ,
   -total_reads    => 10                                  ,
), 'Set the seed';
ok $seed1 = $factory->get_random_seed();
is $seed1, 1233567890;


# Get a seed automatically

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
), 'Get a seed automatically';
ok $seed2 = $factory->get_random_seed();
cmp_ok $seed2, '>', 0;
while (my $read = $factory->next_read) {
   push @dataset1, $read;
}


# Specify the same seed

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
   -random_seed    => $seed2                     ,
), 'Specify the same seed';
ok $seed3 = $factory->get_random_seed();
is $seed3, $seed2;
while (my $read = $factory->next_read) {
   push @dataset2, $read;
}

is_deeply \@dataset1, \@dataset2;



done_testing();
