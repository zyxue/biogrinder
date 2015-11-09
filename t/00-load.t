#! perl

use strict;
use warnings;
use Test::More;


BEGIN {
   use_ok('Grinder');
   use_ok('Grinder::KmerCollection');
   use_ok('Grinder::Database');
   use_ok('t::TestUtils');
}

diag( "Testing Grinder $Grinder::VERSION, Perl $], $^X" );



done_testing();
