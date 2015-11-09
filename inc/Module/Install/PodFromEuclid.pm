#line 1
package Module::Install::PodFromEuclid;


#line 72


use 5.006;
use strict;
use warnings;
use File::Spec;
use Env qw(@INC);
use base qw(Module::Install::Base);

our $VERSION = '0.01';


sub pod_from {
   my ($self, $in_file) = @_;
   return unless $self->is_admin;
   if (not defined $in_file) {
      $in_file = $self->_all_from or die "Error: Could not determine file to make pod_from";
   }
   my @inc = map { ( '-I', File::Spec->rel2abs($_) ) } @INC;
   # use same -I included modules as caller
   my @args = ($^X, @inc, $in_file, '--podfile');
   system(@args) == 0 or die "Error: Could not run command ".join(' ',@args).": $?\n";
   return 1;
}


sub _all_from {
   my $self = shift;
   return unless $self->admin->{extensions};
   my ($metadata) = grep {
      ref($_) eq 'Module::Install::Metadata';
   } @{$self->admin->{extensions}};
   return unless $metadata;
   return $metadata->{values}{all_from} || '';
}


1;

