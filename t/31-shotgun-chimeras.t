#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;


# These tests are to exercise some edge cases of kmer chimeras.
# But let's take the opportunity to do shotgun chimeras

my ($factory, $read);
my %refs;


# Shotgun chimeras

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database_shared_kmers.fa'),
   -chimera_perc    => 100                                     ,
   -chimera_dist    => (1, 1, 1)                               ,
   -chimera_kmer    => 10                                      ,
   -total_reads     => 300                                     ,
   -diversity       => 5                                       ,
), 'Chimera from shotgun library';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   between_ok scalar @refs, 2, 4;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}

is $factory->next_lib, undef;


# Use only some of the sequences

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database_shared_kmers.fa'),
   -chimera_perc    => 100                                     ,
   -chimera_dist    => (1)                                     ,
   -chimera_kmer    => 2                                       ,
   -total_reads     => 300                                     ,
   -diversity       => 3                                       , # 3 out of 5 reference sequences
), 'Use only some of the reference sequences';

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 2;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
is scalar keys %refs, 3;
is $factory->next_lib, undef;


done_testing();
