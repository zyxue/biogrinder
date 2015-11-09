#line 1
package Module::Install::AutoLicense;

use strict;
use warnings;
use base qw(Module::Install::Base);
use vars qw($VERSION);

$VERSION = '0.08';

my %licenses = (
    perl         => 'Software::License::Perl_5',
    apache       => 'Software::License::Apache_2_0',
    artistic     => 'Software::License::Artistic_1_0',
    artistic_2   => 'Software::License::Artistic_2_0',
    lgpl2        => 'Software::License::LGPL_2_1',
    lgpl3        => 'Software::License::LGPL_3_0',
    bsd          => 'Software::License::BSD',
    gpl          => 'Software::License::GPL_1',
    gpl2         => 'Software::License::GPL_2',
    gpl3         => 'Software::License::GPL_3',
    mit          => 'Software::License::MIT',
    mozilla      => 'Software::License::Mozilla_1_1',
);

sub auto_license {
  my $self = shift;
  return unless $Module::Install::AUTHOR;
  my %opts = @_;
  $opts{lc $_} = delete $opts{$_} for keys %opts;
  my $holder = $opts{holder} || _get_authors( $self );
  #my $holder = $opts{holder} || $self->author;
  my $license = $self->license();
  unless ( defined $licenses{ $license } ) {
     warn "No license definition for '$license', aborting\n";
     return 1;
  }
  my $class = $licenses{ $license };
  eval "require $class";
  my $sl = $class->new( { holder => $holder } );
  open LICENSE, '>LICENSE' or die "$!\n";
  print LICENSE $sl->fulltext;
  close LICENSE;
  $self->postamble(<<"END");
distclean :: license_clean

license_clean:
\t\$(RM_F) LICENSE
END

  return 1;
}

sub _get_authors {
  my $self = shift;
  my $joined = join ', ', @{ $self->author() || [] };
  return $joined;
}

'Licensed to auto';
__END__

#line 125
