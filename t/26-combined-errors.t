#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read);


# Combined errors: indels, substitutions, homopolymers, chimeras

ok $factory = Grinder->new(
   -reference_file   => data('shotgun_database_extended.fa'),
   -unidirectional   => 1                                   ,
   -read_dist        => 48                                  ,
   -total_reads      => 1000                                ,
   -homopolymer_dist => 'balzer'                            ,
   -mutation_ratio   => (100, 0)                            ,
   -mutation_dist    => ('uniform', 10)                     ,
   -chimera_perc     => 10                                  ,
   -chimera_dist     => (100)                               ,
   -chimera_kmer     => 0                                   ,
), 'Combined errors (uniform)';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   is $read->id, $nof_reads;
};
is $nof_reads, 1000;


# Combined errors with linear model

ok $factory = Grinder->new(
   -reference_file   => data('shotgun_database_extended.fa'),
   -unidirectional   => 1                                   ,
   -read_dist        => (20, 'normal', 10)                  ,
   -total_reads      => 1000                                ,
   -homopolymer_dist => 'balzer'                            ,
   -mutation_ratio   => (85, 15)                            ,
   -mutation_dist    => ('linear', 2, 2)                    ,
   -chimera_perc     => 10                                  ,
   -chimera_dist     => (100)                               ,
   -chimera_kmer     => 0                                   ,
), 'Combined errors (linear)';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   is $read->id, $nof_reads;
};
is $nof_reads, 1000;


done_testing();
