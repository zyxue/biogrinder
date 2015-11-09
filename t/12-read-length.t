#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, @rlengths, $min, $max, $mean, $stddev, $hist,
    $ehist, $coeff);


# All sequences the same length

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -read_dist      => 50                         ,
   -total_reads    => 1000                       ,
), 'Same length reads';

while ( $read = $factory->next_read ) {
   push @rlengths, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@rlengths);
is $min, 50;
is $max, 50;
is $mean, 50;
is $stddev, 0;
@rlengths = ();


# Uniform distribution
ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -read_dist      => (50, 'uniform', 10)        ,
   -total_reads    => 1000                       ,
), 'Uniform distribution';

while ( $read = $factory->next_read ) {
   push @rlengths, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@rlengths);
is $min, 40;
is $max, 60;
is round($mean), 50;
between_ok( $stddev, 5.3, 6.3 ); # should be 5.79
$hist = hist(\@rlengths, 1, 100);
$ehist = uniform(1, 100, 40, 60, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

SKIP: {
   skip rfit_msg() if not can_rfit();
   test_uniform_dist(\@rlengths, 40, 60);
}

@rlengths = ();


# Normal distribution
ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -read_dist      => (50, 'normal', 5)          ,
   -total_reads    => 1000                       ,
), 'Normal distribution';

while ( $read = $factory->next_read ) {
   push @rlengths, $read->length;
}
($min, $max, $mean, $stddev) = stats(\@rlengths);
between_ok( $mean, 49, 51 ); # should be 50.0
between_ok( $stddev, 4.5, 5.5 ); # should be 5.0
$hist = hist(\@rlengths, 1, 100);
$ehist = normal(1, 100, $mean, $stddev**2, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

SKIP: {
   skip rfit_msg() if not can_rfit();
   test_normal_dist(\@rlengths, 50, 5);
}

@rlengths = ();

done_testing();
