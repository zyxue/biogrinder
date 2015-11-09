#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, $hpols, $min, $max, $mean, $stddev,
    $expected_mean, $expected_stddev, $hist, $ehist, $coeff);

my $delta = 0.20; # 20%
my $min_coeff = 0.75;

# Balzer homopolymer distribution

ok $factory = Grinder->new(
   -reference_file   => data('homopolymer_database.fa'),
   -unidirectional   => 1                              ,
   -read_dist        => 220                            ,
   -total_reads      => 1000                           ,
   -homopolymer_dist => 'balzer'                       ,
), 'Balzer';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
   if ($error_str) {
      unlike $error_str, qr/%/;
      like   $error_str, qr/[+-]/;
   } else {
      ok 1;
      ok 1;
   }
}

for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
   last if $homo_len <= 3; ### TODO: go up to 2
   my $values = $$hpols{$homo_len};
   ($min, $max, $mean, $stddev) = stats($values);
   ($expected_mean, $expected_stddev) = balzer($homo_len);
   #print "Balzer homopolymer length: $homo_len\n";
   #print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";
   #print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
   between_ok( $mean, (1-$delta)*$expected_mean, (1+$delta)*$expected_mean );
   between_ok( $stddev, (1-$delta)*$expected_stddev, (1+$delta)*$expected_stddev );
   $hist = hist($$hpols{$homo_len}, 1, 20);
   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
   $coeff = corr_coeff($hist, $ehist, $mean);
   cmp_ok $coeff, '>', $min_coeff;
   #### TODO: Better test of normality
   #SKIP: {
   #   skip rfit_msg() if not can_rfit();
   #   test_normal_dist($values, $mean, $stddev);
   #}
}
$hpols = {};



# Richter homopolymer distribution

ok $factory = Grinder->new(
   -reference_file   => data('homopolymer_database.fa'),
   -unidirectional   => 1                              ,
   -read_dist        => 220                            ,
   -total_reads      => 1000                           ,
   -homopolymer_dist => 'richter'                      ,
), 'Richter';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
   if ($error_str) {
      unlike $error_str, qr/%/;
      like   $error_str, qr/[+-]/;
   } else {
      ok 1;
      ok 1;
   }
}

for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
   last if $homo_len <= 3; ### TODO: go up to 2
   my $values = $$hpols{$homo_len};
   ($min, $max, $mean, $stddev) = stats($values);
   ($expected_mean, $expected_stddev) = richter($homo_len);
   #print "Richter homopolymer length: $homo_len\n";
   #print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";
   #print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
   between_ok( $mean, (1-$delta)*$expected_mean, (1+$delta)*$expected_mean );
   between_ok( $stddev, (1-$delta)*$expected_stddev, (1+$delta)*$expected_stddev );
   $hist = hist($$hpols{$homo_len}, 1, 20);
   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
   $coeff = corr_coeff($hist, $ehist, $mean);
   cmp_ok $coeff, '>', $min_coeff;
   #### TODO: Better test of normality
   #SKIP: {
   #   skip rfit_msg() if not can_rfit();
   #   test_normal_dist($values, $mean, $stddev);
   #}
}
$hpols = {};



# Margulies homopolymer distribution

ok $factory = Grinder->new(
   -reference_file   => data('homopolymer_database.fa'),
   -unidirectional   => 1                              ,
   -read_dist        => 220                            ,
   -total_reads      => 1000                           ,
   -homopolymer_dist => 'margulies'                    ,
), 'Margulies';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
   if ($error_str) {
      unlike $error_str, qr/%/;
      like   $error_str, qr/[+-]/;
   } else {
      ok 1;
      ok 1;
   }
}

for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
   last if $homo_len <= 3; ### TODO: go up to 2
   ($expected_mean, $expected_stddev) = margulies($homo_len);
   my $values = $$hpols{$homo_len};
   ($min, $max, $mean, $stddev) = stats($values);
   #print "Margulies homopolymer length: $homo_len\n";
   #print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";
   #print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
   between_ok( $mean, (1-$delta)*$expected_mean, (1+$delta)*$expected_mean );
   between_ok( $stddev, (1-$delta)*$expected_stddev, (1+$delta)*$expected_stddev );
   $hist = hist($$hpols{$homo_len}, 1, 20);
   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
   $coeff = corr_coeff($hist, $ehist, $mean);
   cmp_ok $coeff, '>', $min_coeff;
   #### TODO: Better test of normality
   #SKIP: {
   #   skip rfit_msg() if not can_rfit();
   #   test_normal_dist($values, $mean, $stddev);
   #}
}
$hpols = {};

done_testing();




sub add_homopolymers {
   my ($err_str, $ref_seq, $err_h) = @_;

   # Record position and length of homopolymer errors
   my %errors;
   if (defined $err_str) {
      $err_str = combine_dels($err_str, $ref_seq);
      my @errors = split ',', $err_str;
      for my $error (@errors) {
         # Record homopolymer error
         my ($pos, $type, $repl) = ($error =~ m/(\d+)([%+-])(.*)/i);
         my $elen; # error length
         if ($type eq '-') {
            $repl .= '-';
            $elen = - length($repl);
         } elsif ($type eq '+') {
            $elen = + length($repl);
         }
         $errors{$pos} = $elen;
         # Test that proper residue was added to homopolymer
         my $hres    = substr $ref_seq, $pos-1, 1; # first residue of homopolymer
         my $new_res = substr $repl, 0, 1;         # residue to add to homopolymer
         if ( $type eq '+' ) { # in case of insertion (not deletion)
            die if not $hres eq $new_res;
         }
      }
   }

   # Record all homopolymer lengths
   while ( $ref_seq =~ m/(.)(\1+)/g ) {
      # Found a homopolymer
      my $hlen = length($2) + 1;               # length of the error-free homopolymer
      my $pos  = pos($ref_seq) - $hlen + 1;    # start of the homopolymer (residue no.)
      my $elen = $hlen + ($errors{$pos} || 0); # length of the error-containing homopolymer
      push @{$$err_h{$hlen}}, $elen;
   }

   return $err_h;
}



sub combine_dels {
   # Put homopolymer deletions at adjacent position into a single error entry
   # (but only if they affect the same homopolymer)
   # Ex: 45-,46-,47-    becomes 45---
   my ($err_str, $ref_seq) = @_;
   my %errors;

   for my $error (split ',', $err_str) {
      my ($pos, $type, $repl) = ($error =~ m/(\d+)([%+-])([a-z]*)/i);
      if ($type eq '-') {
         # Keep track of what was deleted
         $repl = substr $ref_seq, $pos-1, 1;
         $errors{$pos}{'-'} = [ $repl ];
      } else {
         push @{$errors{$pos}{$type}}, $repl;
      }
   }

   for my $pos (sort {$b <=> $a} (keys %errors)) {
      if (    exists $errors{$pos}{'-'}
          &&  exists $errors{$pos-1}
          &&  exists $errors{$pos-1}{'-'}
          && ($errors{$pos}{'-'}[0] eq $errors{$pos-1}{'-'}[0]) ) {
         push @{$errors{$pos-1}{'-'}}, @{$errors{$pos}{'-'}};
         delete $errors{$pos}{'-'};
         delete $errors{$pos} if scalar keys %{$errors{$pos}} == 0;
      }
   }

   $err_str = '';
   for my $pos (sort {$a <=> $b} (keys %errors)) {
      while ( my ($type, $repls) = each %{$errors{$pos}} ) {
         my $repl;
         if ($type eq '-') {
            $repl = '-' x (scalar @$repls - 1);
         } else {
            $repl = join '', @$repls;
         }
         $err_str .= $pos.$type.$repl.',';
      }
   }
   $err_str =~ s/,$//;
   return $err_str;
}


sub margulies {
   my ($homo_len) = @_;
   my $mean   = $homo_len;
   my $stddev = $homo_len * 0.15;
   return $mean, $stddev;
}


sub richter {
   my ($homo_len) = @_;
   my $mean   = $homo_len;
   my $stddev = sqrt($homo_len) * 0.15;
   return $mean, $stddev;
}


sub balzer {
   my ($homo_len) = @_;
   my $mean     = $homo_len;
   my $variance = 0.03494 + $homo_len * 0.06856;
   return $mean, $stddev;
}
