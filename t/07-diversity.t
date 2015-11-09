#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, %sources);


# Single library, single diversity

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -random_seed    => 1233567880                 ,
   -total_reads    => 100                        ,
   -diversity      => 2                          ,
), 'Single library, single diversity';

while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   $sources{$source} = undef;
};

is scalar keys %sources, 2;
%sources = ();


# Two libraries, single diversity

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -random_seed    => 1233567880                 ,
   -total_reads    => 100                        ,
   -num_libraries  => 2                          ,
   -diversity      => 2                          ,
), 'Two libraries, single diversity';

$factory->next_lib;
while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   $sources{$source} = undef;
};

is scalar keys %sources, 2;
%sources = ();

$factory->next_lib;
while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   $sources{$source} = undef;
};

is scalar keys %sources, 2;
%sources = ();


# Two libraries, two diversities

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -random_seed    => 1233567880                 ,
   -total_reads    => 100                        ,
   -num_libraries  => 2                          ,
   -diversity      => (2, 3)                     ,
), 'Two libraries, two diversities';

$factory->next_lib;
while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   $sources{$source} = undef;
};

is scalar keys %sources, 2;
%sources = ();

$factory->next_lib;
while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   $sources{$source} = undef;
};

is scalar keys %sources, 3;
%sources = ();

done_testing();
