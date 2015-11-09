#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read1, $read2, $nof_reads);
my $total_reads = 100;


# FR-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  +1                         ,
   -mate_orientation => 'FR'                        ,
), 'FR-oriented mates';

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'FR';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  -1                         ,
   -mate_orientation => 'FR'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'FR';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>   0                         ,
   -mate_orientation => 'FR'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read2, $read1), 'FR';
};
is $nof_reads, $total_reads;


# FF-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  +1                         ,
   -mate_orientation => 'FF'                        ,
), 'FF-oriented mates';

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'FF';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  -1                         ,
   -mate_orientation => 'FF'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'RR';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>   0                         ,
   -mate_orientation => 'FF'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   like type($read1, $read2), qr/(FF|RR)/;
};
is $nof_reads, $total_reads;


# RF-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  +1                         ,
   -mate_orientation => 'RF'                        ,
), 'RF-oriented mates';

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'RF';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  -1                         ,
   -mate_orientation => 'RF'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'RF';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>   0                         ,
   -mate_orientation => 'RF'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'RF';
};
is $nof_reads, $total_reads;


# RR-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  +1                         ,
   -mate_orientation => 'RR'                        ,
), 'RR-oriented mates';

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'RR';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>  -1                         ,
   -mate_orientation => 'RR'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   is type($read1, $read2), 'FF';
};
is $nof_reads, $total_reads;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => $total_reads                ,
   -read_dist        =>  80                         ,
   -insert_dist      => 240                         ,
   -unidirectional   =>   0                         ,
   -mate_orientation => 'RR'                        ,
);

$nof_reads = 0;
while ( $read1 = $factory->next_read ) {
   $read2 = $factory->next_read;
   $nof_reads += 2;
   like type($read1, $read2), qr/(FF|RR)/;
};
is $nof_reads, $total_reads;


done_testing();


sub type {
   my ($read1, $read2) = @_;
   my $read1_t = {
      'CCCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' => 'F',
      'tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGGG' => 'R',
   };
   my $read2_t = {
      'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaTTT' => 'F',
      'AAAttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt' => 'R',
   };
   my $type = '';
   if ( (exists $read1_t->{$read1->seq}) && (exists $read2_t->{$read2->seq}) ) {
      $type = $read1_t->{$read1->seq}.$read2_t->{$read2->seq};
   } else {
      $type = $read1_t->{$read2->seq}.$read2_t->{$read1->seq};
   }
   return $type;
}

