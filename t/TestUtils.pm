package t::TestUtils;


use strict;
use warnings;
use POSIX qw( floor ceil );
use Test::More;
use File::Spec::Functions;
use List::Util qw( min max );

use vars qw{@ISA @EXPORT};
BEGIN {
   @ISA     = 'Exporter';
   @EXPORT  = qw{
      PI
      round
      between_ok
      data
      get_references
      get_chars
      stats
      hist
      uniform
      normal
      corr_coeff
      write_data
      can_rfit
      rfit_msg
      error_positions
      test_normal_dist
      test_linear_dist
      test_uniform_dist      
   };
}

our $can_rfit;


#------------------------------------------------------------------------------#


# The Pi mathematical constant
use constant PI => 4 * atan2(1, 1);

sub round {
  # Round the number given as argument
  return int(shift() + 0.5);
}

sub between_ok {
   # Test that a value is in the given range (inclusive)
   my ($value, $min, $max) = @_;
   cmp_ok( $value, '>=', $min ) and
      cmp_ok( $value, '<=', $max ) or
      diag("Got $value but the allowed range was [$min, $max]");
}


sub data {
   # Get the complete filename of a test data file
   return catfile('t', 'data', @_);
}


sub get_references {
   # Get the number of references that a read comes from
   my ($read) = @_;
   my $desc = $read->desc;
   $desc =~ m/reference=(\S+)/;
   my $refs = $1;
   my @refs = split(',', $refs);
   return @refs;
}


sub get_chars {
   # Return a hashref where the keys are the characters seen in the specified
   # string
   my ($string) = @_;
   my %chars;
   for my $pos (0 .. length($string)-1) {
      my $char = substr $string, $pos, 1;
      $chars{$char} = undef;
   }
   my @chars = keys %chars;
   return \%chars;
}


sub stats {
   # Calculates min, max, mean, stddev
   my ($vals) = @_;
   my ($min, $max, $mean, $sum, $sqsum, $stddev) = (1E99, 0, 0, 0, 0, 0);
   my $num = scalar @$vals;
   for my $val (@$vals) {
      $min = $val if $val < $min;
      $max = $val if $val > $max;
      $sum += $val;
      $sqsum += $val**2
   }
   $mean = $sum / $num;
   $stddev = sqrt( $sqsum / $num - $mean**2 );
   return $min, $max, $mean, $stddev;
}


sub hist {
   # Count the number of occurence of each integer:
   # E.g. given the arrayref:
   #    [ 1, 1, 1, 3, 3, 4 ]
   # Return the arrayref:
   #    [ 3, 0, 2, 4 ]
   # The min and the max of the range to consider can be given as an option
   my ($data, $min, $max) = @_;
   if (not defined $data) {
      die "Error: no data provided to hist()\n";
   }
   my %hash;
   for my $val (@$data) {
      $hash{$val}++;
   }
   $min = min(@$data) if not defined $min;
   $max = max(@$data) if not defined $max;
   my @y_data;
   for my $x ($min .. $max) {
      push @y_data, $hash{$x} || 0;
   }
   return \@y_data;
}


sub normal {
   # Evaluate the normal function in the given integer range
   my ($x_min, $x_max, $mean, $variance, $num) = @_;
   my @ys;
   for my $x ($x_min .. $x_max) {
      my $proba = 1 / sqrt(2 * PI * $variance) * exp( - ($x - $mean)**2 / (2 * $variance));
      my $y = $proba * $num;
      push @ys, $y;
   }
   return \@ys;
}


sub uniform {
   # Evaluate the uniform function in the given integer range
   my ($x_min, $x_max, $min, $max, $num) = @_;
   my @ys;
   my $width = $max - $min + 1;
   for my $x ($x_min .. $x_max) {
      my $y;
      if ( ($x >= $min) and ($x <= $max) ) {
         $y = $num / $width;
      } else {
         $y = 0;
      }
      push @ys, $y;
   }
   return \@ys;
}


sub corr_coeff {
   # The correlation coefficient R2 is
   #    R2 = 1 - ( SSerr / SStot )
   # where
   #    SSerr = sum( (y - f)**2 )
   # and
   #    SStot = sum( (y - mean)**2 )
   my ($y, $f, $mean) = @_;
   my $SSerr = 0;
   my $SStot = 0;
   for my $i ( 0 .. scalar @$y - 1 ) {
      #print "  ".($i+1)."  ".$$y[$i]."  ".$$f[$i]."\n";
      $SSerr += ($$y[$i] - $$f[$i])**2;
      $SStot += ($$y[$i] - $mean)**2;
   }
   my $R2 = 1 - ($SSerr / $SStot);
   return $R2;
}


sub write_data {
   # Write a data series (array reference) to a file with the specified name, or
   # 'data.txt' by default
   my ($data, $filename) = @_;
   $filename = 'data.txt' if not defined $filename;
   open my $out, '>', $filename or die "Error: Could not write file $filename\n$!\n";
   for my $datum (@$data) {
      print $out "$datum\n";
   }
   close $out;
   return $filename;
}


sub can_rfit {
   # Determine if a system can run the fitdistrplus R module through the
   # Statistics::R Perl interface. Load Statistics::R if it can and return 1.
   # Return 0 otherwise.
   if (not defined $can_rfit) {
      eval {
         require Statistics::R;
         my $R = Statistics::R->new();
         my $ret = $R->run(q`library(fitdistrplus)`);
         $R->stop();
      };
      if ($@) {
         $can_rfit = 0;
         my $msg = "Note: The Statistics::R module for Perl, R (R-Project) ".
            "or the fitdistrplus module for R could not be found on this system.".
            " Some tests will be skipped...\n";
         warn $msg;
      } else {
         $can_rfit = 1;
      }
   }
   return $can_rfit;
}


sub rfit_msg {
   return "fitdistrplus not available...";
}


sub error_positions {
   my ($read) = @_;
   my ($err_str) = ($read->desc =~ /errors=(\S+)/);
   my @error_positions;
   if (defined $err_str) {
      for my $error (split ',', $err_str) {
         my ($pos, $type, $res) = ($error =~ m/(\d+)([%+-])([a-z]*)/i);
         push @error_positions, $pos;
      }
   }
   return @error_positions;
}


sub test_linear_dist {
   # Test that the datapoints provided follow a linear distribution
   my ($values, $want_min, $want_max, $want_slope) = @_;

   my ($min, $max, $ratio_lo, $ratio_hi, $slope, $chisqpvalue, $chisqtest) =
      fit_linear($values);

   is $want_min, $min, 'fitdist() linear';
   is $want_max, $max;

   between_ok( 2, $ratio_lo, $ratio_hi );
 
   between_ok( $slope, (1 - 0.05) * $want_slope, (1 + 0.05) * $want_slope ); # Allow a 5% standard deviation

   is( $chisqtest, 'not rejected', 'Chi square test') or
      diag("p-value was: $chisqpvalue");

   return 1;
}


sub test_uniform_dist {
   # Test that the integer series provided follow a uniform distribution with the 
   # specified minimum and maximum. Note that you probably need over 30-100
   # values for the statistical test to work!
   my ($values, $want_min, $want_max) = @_;

   my ($min_lo, $min_hi, $max_lo, $max_hi, $chisqpvalue, $chisqtest) =
      fit_uniform($values, $want_min, $want_max);

   # Need to be more lenient since fitdistrplus is not too good with integers
   #between_ok( $want_min, $min_lo, $min_hi );
   #between_ok( $want_max, $max_lo, $max_hi );
   between_ok( round($want_min), floor($min_lo), ceil($min_hi) );
   between_ok( round($want_max), floor($max_lo), ceil($max_hi) );

   is( $chisqtest, 'not rejected', 'Chi square test') or
      diag("p-value was: $chisqpvalue");

   return 1;
}


sub test_normal_dist {
   # Test that the integer series provided follow a normal distribution with the 
   # specified mean and standard deviation. Note that you probably need over
   # 30-100 values for the statistical test to work!
   my ($values, $want_mean, $want_sd, $filename) = @_;

   my ($mean_lo, $mean_hi, $sd_lo, $sd_hi, $chisqpvalue, $chisqtest) =
      fit_normal($values, $want_mean, $want_sd);

   # Need to be more lenient since fitdistrplus is not too good with integers
   #between_ok( $want_mean, $mean_lo, $mean_hi );
   #between_ok( $want_sd  , $sd_lo  , $sd_hi   );
   between_ok( round($want_mean), floor($mean_lo), ceil($mean_hi) );
   between_ok( round($want_sd  ), floor($sd_lo  ), ceil($sd_hi  ) );

   is( $chisqtest, 'not rejected', 'Chi square test') or 
      diag("p-value was: $chisqpvalue");

   return 1;
}


#------------------------------------------------------------------------------#


my $niter = 30; # number of iterations to fit the distributions


sub fit_linear {
   my ($values) = @_;

   # Fit a linear distribution. Since R does not have a linear distribution, use
   # the beta distribution:
   #    when beta shape1=1 & shape2=2, distribution is linearly decreasing (slope=-2)
   #    when beta shape1=2 & shape2=1, distribution is linearly increasing (slope=2)

   # Find min and max of series
   my $min = min(@$values);
   my $max = max(@$values);

   # Rescale values between in 0 and 1 instead of min and max
   my $rescaled_values;
   if ( ($min == 0) and ($max == 1) ) {
      $rescaled_values = $values;
   } else {
      for my $value (@$values) {
         push @$rescaled_values, ($value - $min) / ($max - $min);
      }
   }

   # Now we can run fit_beta()
   my ($shape1_lo, $shape1_hi, $shape2_lo, $shape2_hi, $chisqpvalue, $chisqtest)
      = fit_beta($values, 2, 1, $min, $max);

   my $ratio_hi = $shape1_hi / $shape2_lo;
   my $ratio_lo = $shape2_hi / $shape1_lo;

   # Calculate the slope
   my $slope = 2 / ($max - $min);

   return $min, $max, $ratio_lo, $ratio_hi, $slope, $chisqpvalue, $chisqtest;
}


sub fit_beta {
   # Try to fit a beta distribution to a series of data points using a maximum
   # goodness of fit method. Return the 95% confidence interval for the shape1
   # parameter, the shape2 parameter and the results of Chi square statistics.
   my ($values, $want_shape1, $want_shape2, $want_min, $want_max) = @_;
   my $break_num   = $want_max - $want_min;
   my $break_size  = 1 / $break_num;
   my $break_start = 0 - $break_size / 2;
   my $break_end   = 1 + $break_size / 2;
   my $start_p   = "start=list(shape1=$want_shape1, shape2=$want_shape2)";
   my $breaks_p  = "chisqbreaks=seq($break_start, $break_end, $break_size)";
   #my $fit_cmd   = "f  <- fitdist(x, distr='beta', method='mle', $start_p)";
   my $fit_cmd   = "f  <- fitdist(x, distr='beta', method='mge', gof='CvM', $start_p)";
   my $boot_cmd  = "fb <- bootdist(f, niter=$niter)";
   my $gof_cmd   = "g  <- gofstat(f, $breaks_p)";
   my $R = Statistics::R->new();
   $R->set('x', $values);
   $R->run('library(fitdistrplus)');
   $R->run($fit_cmd);
   $R->run($boot_cmd);
   $R->run($gof_cmd);
   my $shape1_lo   = $R->get('fb$CI[1,2]');
   my $shape1_hi   = $R->get('fb$CI[1,3]');
   my $shape2_lo   = $R->get('fb$CI[2,2]');
   my $shape2_hi   = $R->get('fb$CI[2,3]');
   my $chisqpvalue = $R->get('g$chisqpvalue');
   my $chisqtest   = test_result($chisqpvalue);
   $R->stop();
   return $shape1_lo, $shape1_hi, $shape2_lo, $shape2_hi, $chisqpvalue, $chisqtest;
}


sub fit_uniform {
   # Try to fit a uniform distribution to a series of integers using a maximum
   # goodness of fit method. Return the 95% confidence interval for the mean,
   # the standard deviation and the results of the Chi square statistics.
   my ($values, $want_min, $want_max) = @_;
   my $range_min = min(@$values) - 0.5;
   my $range_max = max(@$values) + 0.5;
   my $breaks_p  = "chisqbreaks=seq($range_min, $range_max)";
   my $start_p   = "start=list(min=$want_min, max=$want_max)";
   my $fit_cmd   = "f  <- fitdist(x, distr='unif', method='mge', gof='CvM', $start_p)";
   my $boot_cmd  = "fb <- bootdist(f, niter=$niter)";
   my $gof_cmd   = "g  <- gofstat(f, $breaks_p)";
   my $R = Statistics::R->new();
   $R->set('x', $values);
   $R->run('library(fitdistrplus)');
   $R->run($fit_cmd);
   $R->run($boot_cmd);
   $R->run($gof_cmd);
   my $min_lo      = $R->get('fb$CI[1,2]');
   my $min_hi      = $R->get('fb$CI[1,3]');
   my $max_lo      = $R->get('fb$CI[2,2]');
   my $max_hi      = $R->get('fb$CI[2,3]');
   my $chisqpvalue = $R->get('g$chisqpvalue');
   my $chisqtest   = test_result($chisqpvalue);
   $R->stop();
   return $min_lo, $min_hi, $max_lo, $max_hi, $chisqpvalue, $chisqtest;
}


sub fit_normal {
   # Try to fit a normal distribution to a series of integers using a maximum
   # likelihood method. Return the 95% confidence interval for the mean, the
   # standard deviation and the results of the Chi square statistics.
   my ($values, $want_mean, $want_sd) = @_;
   my $range_min = min(@$values) - 0.5;
   my $range_max = max(@$values) + 0.5;
   my $breaks_p  = "chisqbreaks=seq($range_min, $range_max)";
   my $start_p   = "start=list(mean=$want_mean, sd=$want_sd)";
   my $fit_cmd   = "f  <- fitdist(x, distr='norm', method='mle', $start_p)";
   my $boot_cmd  = "fb <- bootdist(f, niter=$niter)";
   my $gof_cmd   = "g  <- gofstat(f, $breaks_p)";
   my $R = Statistics::R->new();
   $R->set('x', $values);
   $R->run('library(fitdistrplus)');
   $R->run($fit_cmd);
   $R->run($boot_cmd);
   $R->run($gof_cmd);
   my $mean_lo     = $R->get('fb$CI[1,2]');
   my $mean_hi     = $R->get('fb$CI[1,3]');
   my $sd_lo       = $R->get('fb$CI[2,2]');
   my $sd_hi       = $R->get('fb$CI[2,3]');
   my $chisqpvalue = $R->get('g$chisqpvalue');
   my $chisqtest   = test_result($chisqpvalue);
   $R->stop();
   return $mean_lo, $mean_hi, $sd_lo, $sd_hi, $chisqpvalue, $chisqtest;
}


sub test_result {
   # Reject a statistical test if the p value is less than 0.05
   my ($p_value) = @_;
   my $test_result;
   if ( lc $p_value eq 'nan' ) {
      $p_value = 1; # probably a very large p value
   }
   my $thresh = 0.05;
   if ($p_value <= $thresh) {
      $test_result = 'rejected';
   } elsif ($p_value > $thresh) {
      $test_result = 'not rejected';
   } else {
      die "Error: '$p_value' is not a supported p-value\n";
   }
   return $test_result;
}


1;
