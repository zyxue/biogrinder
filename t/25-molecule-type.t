#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read);

my $dna_want_chars = {
   'A' => undef,
   'C' => undef,
   'G' => undef,
   'T' => undef,
};

my $rna_want_chars = {
   'A' => undef,
   'C' => undef,
   'G' => undef,
   'U' => undef,
};

my $protein_want_chars = {
   'A' => undef,
   'R' => undef,
   'N' => undef,
   'D' => undef,
   'C' => undef,
   'Q' => undef,
   'E' => undef,
   'G' => undef,
   'H' => undef,
   'I' => undef,
   'L' => undef,
   'K' => undef,
   'M' => undef,
   'F' => undef,
   'P' => undef,
   'S' => undef,
   'T' => undef,
   'W' => undef,
   'Y' => undef,
   'V' => undef,
};


# DNA database

ok $factory = Grinder->new(
   -reference_file  => data('database_dna.fa'),
   -read_dist       => 240                    ,
   -total_reads     => 100                    ,
   -mutation_ratio => (100, 0)                ,
   -mutation_dist  => ('uniform', 20)         ,
), 'DNA';

is $factory->{alphabet}, 'dna';

while ($read = $factory->next_read) {
   my $got_chars = get_chars($read->seq);
   is_deeply $got_chars, $dna_want_chars;
}


# RNA

ok $factory = Grinder->new(
   -reference_file  => data('database_rna.fa'),
   -read_dist       => 240                    ,
   -total_reads     => 100                    ,
   -mutation_ratio => (100, 0)                ,
   -mutation_dist  => ('uniform', 20)         ,
), 'RNA';

is $factory->{alphabet}, 'rna';

while ($read = $factory->next_read) {
   my $got_chars = get_chars($read->seq);
   is_deeply $got_chars, $rna_want_chars;
}


# Protein

ok $factory = Grinder->new(
   -reference_file  => data('database_protein.fa'),
   -total_reads     => 100                        ,
   -read_dist       => 240                        ,
   -mutation_ratio => (100, 0)                    ,
   -mutation_dist  => ('uniform', 20)             ,
   -unidirectional  => +1                         ,
), 'Protein';

is $factory->{alphabet}, 'protein';

while ($read = $factory->next_read) {
   my $got_chars = get_chars($read->seq);
   is_deeply $got_chars, $protein_want_chars;
}


# Mixed

ok $factory = Grinder->new(
   -reference_file  => data('database_mixed.fa'),
   -total_reads     => 100                      ,
   -unidirectional  => +1                       ,
), 'Mixed';

is $factory->{alphabet}, 'protein';

done_testing();
