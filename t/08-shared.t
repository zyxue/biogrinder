#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, $lib_num, %sources, %shared);


# No species shared

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -length_bias     => 0                          ,
   -total_reads     => 100                        ,
   -num_libraries   => 3                          ,
   -shared_perc     => 0                          ,
), 'No species shared';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 3               &&
                                    exists $sources{1}{$source} &&
                                    exists $sources{2}{$source}   );
   }
};

is scalar keys %sources, 3;
is scalar keys %{$sources{1}}, 1;
is scalar keys %{$sources{2}}, 1;
is scalar keys %{$sources{3}}, 1;
is scalar keys %shared, 0;
%sources = ();
%shared  = ();


# 50% species shared

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -length_bias     => 0                          ,
   -total_reads     => 100                        ,
   -num_libraries   => 3                          ,
   -shared_perc     => 50                         ,
), '50% species shared';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 3               &&
                                    exists $sources{1}{$source} &&
                                    exists $sources{2}{$source}   );
   }
};

is scalar keys %sources, 3;
is scalar keys %{$sources{1}}, 2;
is scalar keys %{$sources{2}}, 2;
is scalar keys %{$sources{3}}, 2;
is scalar keys %shared, 1;
%sources = ();
%shared  = ();


# 66% species shared

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -length_bias     => 0                          ,
   -total_reads     => 100                        ,
   -num_libraries   => 3                          ,
   -shared_perc     => 66                         ,
), '66% species shared';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 3               &&
                                    exists $sources{1}{$source} &&
                                    exists $sources{2}{$source}   );
   }
};

is scalar keys %sources, 3;
is scalar keys %{$sources{1}}, 2;
is scalar keys %{$sources{2}}, 2;
is scalar keys %{$sources{3}}, 2;
is scalar keys %shared, 1;
%sources = ();
%shared  = ();


# 67% species shared

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -length_bias     => 0                          ,
   -total_reads     => 100                        ,
   -num_libraries   => 3                          ,
   -shared_perc     => 67                         ,
), '67% species shared';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 3               &&
                                    exists $sources{1}{$source} &&
                                    exists $sources{2}{$source}   );
   }
};

is scalar keys %sources, 3;
is scalar keys %{$sources{1}}, 3;
is scalar keys %{$sources{2}}, 3;
is scalar keys %{$sources{3}}, 3;
is scalar keys %shared, 2;
%sources = ();
%shared  = ();


# All species shared

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -length_bias     => 0                          ,
   -total_reads     => 100                        ,
   -num_libraries   => 3                          ,
   -shared_perc     => 100                        ,
), 'All species shared';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 3               &&
                                    exists $sources{1}{$source} &&
                                    exists $sources{2}{$source}   );
   }
};

is scalar keys %sources, 3;
is scalar keys %{$sources{1}}, 5;
is scalar keys %{$sources{2}}, 5;
is scalar keys %{$sources{3}}, 5;
is scalar keys %shared, 5;
%sources = ();
%shared  = ();


# Inequal richness

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -random_seed     => 1233567890                 ,
   -abundance_model => 'uniform'                  ,
   -total_reads     => 100                        ,
   -length_bias     => 0                          ,
   -num_libraries   => 2                          ,
   -diversity       => (3,5)                      ,
   -shared_perc     => 100                        ,
), 'Inequal richness';

while ($factory->next_lib) {
   $lib_num = $factory->{cur_lib};
   while ( $read = $factory->next_read ) {
      my $source = $read->reference->id;
      $sources{$lib_num}{$source} = undef;
      # Is this source genome shared?
      $shared{$source} = undef if ( $lib_num == 2               &&
                                    exists $sources{1}{$source}    );
   }
};

is scalar keys %sources, 2;
is scalar keys %{$sources{1}}, 3;
is scalar keys %{$sources{2}}, 5;
is scalar keys %shared, 3;
%sources = ();
%shared  = ();

done_testing();

