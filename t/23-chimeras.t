#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads, $nof_chimeras, $nof_regulars);
my %chim_sizes;


# No Chimeras

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 0                                 ,
   -chimera_dist    => (1)                               ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 100                               ,
), 'No chimeras';

while ( $read = $factory->next_read ) {
   is scalar get_references($read), 1;
   # Remove forward and reverse primer
   my $seq = $read->seq;
   $seq = remove_primers($seq, 'AAACT.AAA.GAATTG.CGG', 'G.ACACACCGCCCGT');
   # Now the amplicon is simply long homopolymeric sequences
   like $seq, qr/^(a+|c+|g+|t+)+$/;
}


# 50% chimeras (bimeras)

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 50                                ,
   -chimera_dist    => (1)                               ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 100                               ,
), '50% chimeras (bimeras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   $nof_chimeras += scalar get_references($read);
   $nof_regulars += scalar get_references($read);
}
between_ok( $nof_chimeras / $nof_regulars, 0.9, 1.1 );


# 100% chimeras (bimeras)

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -chimera_dist    => (1)                               ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 100                               ,
), '100% chimeras (bimeras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   is scalar get_references($read), 2;
}


# 100% chimeras (trimeras)

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -chimera_dist    => (0, 1)                            ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 100                               ,
), '100% chimeras (trimeras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   is scalar get_references($read), 3;
}


# 100% chimeras (quadrameras)

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -chimera_dist    => (0, 0, 1)                         ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 100                               ,
), '100% chimeras (quadrameras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   is scalar get_references($read), 4;
}


# 100% chimeras (bimeras, trimeras, quadrameras)

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -chimera_dist    => (1, 1, 1)                         ,
   -chimera_kmer    => 0                                 ,
   -total_reads     => 1000                              ,
), '100% chimeras (bimeras, trimeras, quadrameras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   my $nof_refs = scalar get_references($read);
   $chim_sizes{$nof_refs}++;
   between_ok( $nof_refs, 2, 4 );
}
between_ok( $chim_sizes{2}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{3}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{4}, 333.3 * 0.9, 333.3 * 1.1 );


done_testing();





sub remove_primers {
   my ($seq, $forward_re, $reverse_re) = @_;
   $seq =~ s/$forward_re//i;
   $seq =~ s/$reverse_re//i;
   return $seq;
}


sub matches_ref {
   my ($read) = @_;
   my $read_seq = $read->seq;
   my $ref_seq  = $read->reference->seq;
   my $matches  = 0;
   if ($ref_seq =~ m/$read_seq/) {
      $matches = 1;
      print "$read_seq\nmatches\n$ref_seq\n\n";
   }
   return $matches;
}
