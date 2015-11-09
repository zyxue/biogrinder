#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Bio::PrimarySeq;


use_ok 'Grinder::KmerCollection';
my ($col, $seq1, $seq2, $by_kmer, $by_seq, $file, $sources, $counts, $freqs,
    $kmers, $pos, $weights);


# Test the Grinder::KmerCollection module

$seq1 = Bio::PrimarySeq->new(
  -id => 'seq1',
  -seq => 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
);

$seq2 = Bio::PrimarySeq->new(
  -id => 'seq4',
  -seq => 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGCCCCCCCC'
);


ok $col = Grinder::KmerCollection->new( -k => 8 );

isa_ok $col, 'Grinder::KmerCollection';
is $col->k, 8;

ok $col->add_seqs([$seq1]);
ok $col->add_seqs([$seq2]);

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok exists $by_kmer->{'CCCCCCCC'}->{'seq4'};
ok exists $by_kmer->{'CCCCGGGG'}->{'seq4'};
ok exists $by_kmer->{'ACCCCCCC'}->{'seq4'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok exists $by_kmer->{'seq4'}->{'ACCCCCCC'};

ok $col = $col->filter_rare(2);
isa_ok $col, 'Grinder::KmerCollection';

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok exists $by_kmer->{'CCCCCCCC'}->{'seq4'};
ok not exists $by_kmer->{'CCCCGGGG'};
ok not exists $by_kmer->{'ACCCCCCC'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok not exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok not exists $by_kmer->{'seq4'}->{'ACCCCCCC'};


ok $col = Grinder::KmerCollection->new( -k => 8, -seqs => [$seq1, $seq2] );

# Count of all kmers
($kmers, $counts) = $col->counts();
$kmers  = [sort @$kmers];
$counts = [sort {$a <=> $b} @$counts];
is_deeply $kmers , [
          'AAAAAAAA',
          'AAAAAAAC',
          'AAAAAACC',
          'AAAAACCC',
          'AAAACCCC',
          'AAACCCCC',
          'AACCCCCC',
          'ACCCCCCC',
          'CAAAAAAA',
          'CCAAAAAA',
          'CCCAAAAA',
          'CCCCAAAA',
          'CCCCCAAA',
          'CCCCCCAA',
          'CCCCCCCA',
          'CCCCCCCC',
          'CCCCCCCG',
          'CCCCCCGG',
          'CCCCCGGG',
          'CCCCGGGG',
          'CCCGGGGG',
          'CCGGGGGG',
          'CGGGGGGG',
          'GCCCCCCC',
          'GGCCCCCC',
          'GGGCCCCC',
          'GGGGCCCC',
          'GGGGGCCC',
          'GGGGGGCC',
          'GGGGGGGC',
          'GGGGGGGG'
];
is_deeply $counts, [
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          43,
          74,
];

# Frequency of kmers from position >= 40
($kmers, $freqs) = $col->counts(undef, 40, 1); 
$kmers = [sort @$kmers];
$freqs = [sort {$a <=> $b} @$freqs];
is_deeply $kmers , [
          'AAAAAAAA',
          'AACCCCCC',
          'ACCCCCCC',
          'CCCCCCCC',
          'CCCCCCCG',
          'CCCCCCGG',
          'CCCCCGGG',
          'CCCCGGGG',
          'CCCGGGGG',
          'CCGGGGGG',
          'CGGGGGGG',
          'GCCCCCCC',
          'GGCCCCCC',
          'GGGCCCCC',
          'GGGGCCCC',
          'GGGGGCCC',
          'GGGGGGCC',
          'GGGGGGGC',
          'GGGGGGGG'
];
is_deeply $freqs, [
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.0147058823529412',
          '0.25',
          '0.5',
];

($kmers, $freqs) = $col->counts('seq1', 40, 1); 
is_deeply $kmers, [ 'AAAAAAAA' ];
is_deeply $freqs, [ 1 ];

($kmers, $freqs) = $col->counts('seq4', 40, 1); 
$kmers = [sort @$kmers];
$freqs = [sort {$a <=> $b} @$freqs];
is_deeply $kmers , [
          'AACCCCCC',
          'ACCCCCCC',
          'CCCCCCCC',
          'CCCCCCCG',
          'CCCCCCGG',
          'CCCCCGGG',
          'CCCCGGGG',
          'CCCGGGGG',
          'CCGGGGGG',
          'CGGGGGGG',
          'GCCCCCCC',
          'GGCCCCCC',
          'GGGCCCCC',
          'GGGGCCCC',
          'GGGGGCCC',
          'GGGGGGCC',
          'GGGGGGGC',
          'GGGGGGGG'
];
is_deeply $freqs, [
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.0294117647058824',
          '0.5',
];

ok $col = $col->filter_shared(2);
isa_ok $col, 'Grinder::KmerCollection';

($kmers, $counts) = $col->counts();
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [       74 ];

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok not exists $by_kmer->{'CCCCCCCC'};
ok not exists $by_kmer->{'CCCCGGGG'};
ok not exists $by_kmer->{'ACCCCCCC'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok not exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok not exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok not exists $by_kmer->{'seq4'}->{'ACCCCCCC'};

($sources, $counts) = $col->sources('AAAAAAAA');
is_deeply $sources, ['seq4', 'seq1'];
is_deeply $counts , [    1 ,    73 ];

($sources, $counts) = $col->sources('AAAAAAAA', 'seq1');
is_deeply $sources, ['seq4'];
is_deeply $counts , [    1 ];

($sources, $counts) = $col->sources('ZZZZZZZZ');
is_deeply $sources, [];
is_deeply $counts , [];

($kmers, $counts) = $col->kmers('seq1');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [       73 ];

($kmers, $counts) = $col->kmers('seq4');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [        1 ];

($kmers, $counts) = $col->kmers('asdf');
is_deeply $kmers , [];
is_deeply $counts, [];

$pos = $col->positions('AAAAAAAA', 'seq1');
is_deeply $pos, [1..73];

$pos = $col->positions('AAAAAAAA', 'seq4');
is_deeply $pos, [34];

$pos = $col->positions('CCCCGGGG', 'seq4');
is_deeply $pos, [];

$pos = $col->positions('AAAAAAAA', 'seq3');
is_deeply $pos, [];


ok $col = Grinder::KmerCollection->new(
   -k    => 8,
   -seqs => [$seq1, $seq2], 
   -ids  => ['abc', '123'],
)->filter_rare(2);

isa_ok $col, 'Grinder::KmerCollection';

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'abc'};
ok exists $by_kmer->{'AAAAAAAA'}->{'123'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'abc'}->{'AAAAAAAA'};
ok exists $by_kmer->{'123'}->{'AAAAAAAA'};

($sources, $counts) = $col->sources('AAAAAAAA');
is_deeply $sources, ['123', 'abc'];
is_deeply $counts , [    1 ,    73 ];

($sources, $counts) = $col->sources('AAAAAAAA', 'abc');
is_deeply $sources, ['123'];
is_deeply $counts , [    1 ];


# Using weights

ok $col = Grinder::KmerCollection->new(
   -k => 8, -seqs => [$seq1, $seq2],
)->filter_shared(2);

$weights = { 'seq1' => 10, 'seq4' => 0.1 };
ok $col->weights($weights);

($sources, $counts) = $col->sources('AAAAAAAA');
is_deeply $sources, ['seq4', 'seq1'];
is_deeply $counts , [  0.1 ,   730 ];

($kmers, $counts) = $col->counts();
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [    730.1 ];

($kmers, $counts) = $col->kmers('seq1');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [      730 ];

($kmers, $counts) = $col->kmers('seq4');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [      0.1 ];

ok $col->weights({});
is_deeply $col->weights, {};


# Read from file

$file = data('kmers.fa');
ok $col = Grinder::KmerCollection->new( -k => 8, -file => $file, );


done_testing;
