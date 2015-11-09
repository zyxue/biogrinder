#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, $lib_num, $ranks1, $ranks2, $ranks3, $rank1_perm);


# No species permuted

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 0                          ,
), 'No species permuted';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 5;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 0;
compare_ranks( $ranks1, $ranks2, $rank1_perm );


# Cannot have 20% permuted (1 species) because in this permutation method,
# top species are permuted amongst the top species


# 40% species permuted (2 species permuted)
# This test is very sensitive to the seed because when permuting the top 2
# species, 50% of the time, the answer will be (1,2) and the rest of the time,
# it will be (2,1)

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 2183567890                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 40                         ,
), '40% species permuted';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 5;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 2;
compare_ranks( $ranks1, $ranks2, $rank1_perm );


# 60% species permuted

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1095230708                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 60                         ,
), '60% species permuted';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 5;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 3;
compare_ranks( $ranks1, $ranks2, $rank1_perm );


# 80% species permuted

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1095230708                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 80                         ,
), '80% species permuted';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 5;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 4;
compare_ranks( $ranks1, $ranks2, $rank1_perm );


# All species permuted

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1933067890                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 100                        ,
), 'All species permuted';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 5;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 5;
compare_ranks( $ranks1, $ranks2, $rank1_perm );


# Inequal richness

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1243567820                 ,
   -abundance_model => ('powerlaw', 1.8)          ,
   -total_reads     => 1000                       ,
   -num_libraries   => 2                          ,
   -length_bias     => 0                          ,
   -diversity       => (3,5)                      ,
   -shared_perc     => 100                        ,
   -permuted_perc   => 100                        ,
), 'Inequal richness';

ok $factory->next_lib;
$ranks1 = get_ranks($factory);
is scalar @$ranks1, 3;
ok $factory->next_lib;
$ranks2 = get_ranks($factory);
is scalar @$ranks2, 5;

$rank1_perm = 3;
compare_ranks( $ranks1, $ranks2, $rank1_perm );



sub compare_ranks {
   # Compare genome ranks 2 to genome ranks 1
   my ($ranks1, $ranks2, $rank1_perm) = @_;

   # Top genomes that should be permuted
   my @perm_ids;
   if ($rank1_perm == 0) {
      # nothing to do
   } elsif ($rank1_perm > 0) {
      @perm_ids = @$ranks1[0..$rank1_perm-1];
   }   

   # Copy arrays
   my %refs1;
   for my $rank (1 .. scalar @$ranks1) {
      $refs1{$$ranks1[$rank-1]} = $rank;
   }
   my %refs2;
   for my $rank (1 .. scalar @$ranks2) {
      $refs2{$$ranks2[$rank-1]} = $rank;
   }

   # Test that permuted genomes have a different rank
   # Note: This is not foolproof because at high percentage permuted, on samples
   #       with a small richness, the high number of permutations can cause a 
   #       permuted genome to end up with the same rank as initially. Need to 
   #       play with the seed number to get it right.
   for my $perm_id ( @perm_ids ) {
      my $rank1 = $refs1{$perm_id};
      my $rank2 = $refs2{$perm_id};
      isnt $rank1, $rank2;
      delete $refs1{$perm_id};
      delete $refs2{$perm_id};
   }

   # Now, remaining genomes should have identical ranks (because it is 100% shared)
   my @refs1 = sort { $refs1{$a} <=> $refs1{$b} } (keys %refs1);  
   my @refs2 = sort { $refs2{$a} <=> $refs2{$b} } (keys %refs2);

   # Length of the smallest array (number of species shared is relative to )
   my $min_arr_len = scalar @refs1 < scalar @refs2 ? scalar @refs1 : scalar @refs2;
   for my $i (0 .. $min_arr_len - 1) {
      my $id1 = $refs1[$i];
      my $id2 = $refs2[$i];
      is $id1, $id2;
   } 

   return 1;
}

done_testing();



sub get_ranks {
   my ($factory) = @_;
   my %sources;
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      if (not exists $sources{$source}) {
         $sources{$source} = 1;
      } else {
         $sources{$source}++;
      }
   }
   my @ranks = sort { $sources{$b} <=> $sources{$a} } (keys %sources);
   return \@ranks
}
