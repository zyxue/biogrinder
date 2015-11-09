#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);
my %chim_sizes;
my %refs;
my %expected;
my $delta = 0.2;


# Bimeras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (1)             ,
   -chimera_kmer    => 8               ,
   -total_reads     => 300             ,
), 'Bimeras';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 2;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}

ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


# Trimeras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (0, 1)          ,
   -chimera_kmer    => 8               ,
   -total_reads     => 300             ,
), 'Trimeras';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 3;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


# Quadrameras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (0, 0, 1)       ,
   -chimera_kmer    => 8               ,
   -total_reads     => 300             ,
), 'Quadrameras';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 4;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


# 100% chimeras (bimeras, trimeras, quadrameras)

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (1, 1, 1)       ,
   -chimera_kmer    => 8               ,
   -total_reads     => 1000            ,
), '100% chimeras (bimeras, trimeras, quadrameras)';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   my $nof_refs = scalar @refs;
   $chim_sizes{$nof_refs}++;
   between_ok( $nof_refs, 2, 4 );
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
between_ok( $chim_sizes{2}, 333.3 * (1-$delta), 333.3 * (1+$delta) );
between_ok( $chim_sizes{3}, 333.3 * (1-$delta), 333.3 * (1+$delta) );
between_ok( $chim_sizes{4}, 333.3 * (1-$delta), 333.3 * (1+$delta) );
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


# From equal abundance sequences

ok $factory = Grinder->new(
   -reference_file  => data('kmers2.fa'),
   -length_bias     => 0                ,
   -unidirectional  => 1                ,
   -chimera_perc    => 100              ,
   -chimera_dist    => (0, 0, 0, 0, 1)  ,
   -chimera_kmer    => 8                ,
   -total_reads     => 1000             ,
), 'From equal abundance sequences';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 6;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}

%expected = ( 'seq1' => 6000 * 4/18,
              'seq2' => 6000 * 6/18,
              'seq3' => 6000 * 8/18, );

between_ok $refs{'seq1'}, $expected{'seq1'}*(1-$delta), $expected{'seq1'}*(1+$delta);
between_ok $refs{'seq2'}, $expected{'seq2'}*(1-$delta), $expected{'seq2'}*(1+$delta);
between_ok $refs{'seq3'}, $expected{'seq3'}*(1-$delta), $expected{'seq3'}*(1+$delta);


# From differentially abundant sequences

ok $factory = Grinder->new(
   -reference_file  => data('kmers2.fa')          ,
   -abundance_file  => data('abundance_kmers.txt'),
   -length_bias     => 0                          ,
   -unidirectional  => 1                          ,
   -chimera_perc    => 100                        ,
   -chimera_dist    => (0, 0, 0, 0, 1)            ,
   -chimera_kmer    => 8                          ,
   -total_reads     => 1000                       ,
), 'From differentially abundant sequences';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 6;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}

$delta = 0.2;
cmp_ok $refs{'seq2'}, '<', 1100;
# seq1 and seq3 should occur as frequently
cmp_ok $refs{'seq1'}, '>', 2300;
cmp_ok $refs{'seq3'}, '>', 2300;
between_ok $refs{'seq1'}, $expected{'seq3'}*(1-$delta), $expected{'seq3'}*(1+$delta);


done_testing();
