package Grinder::Database;

use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root); # using throw() and _rearrange() methods


sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($fasta_file, $unidirectional, $primers, $abundance_file, $delete_chars,
      $minimum_length) = $self->_rearrange([qw(FASTA_FILE UNIDIRECTIONAL PRIMERS
      ABUNDANCE_FILE DELETE_CHARS MINIMUM_LENGTH)], @args);

   $minimum_length = 1  if not defined $minimum_length;
   $self->_set_minimum_length($minimum_length);

   $delete_chars = '' if not defined $delete_chars;
   $self->_set_delete_chars($delete_chars);

   # Index file, filter sequences and get IDs
   $self->_init_db($fasta_file, $abundance_file, $delete_chars, $minimum_length);

   $unidirectional = 0 if not defined $unidirectional; # bidirectional
   $self->_set_unidirectional($unidirectional);

   # Read amplicon primers
   $self->_set_primers($primers) if defined $primers;

   # Error if trying to reverse complement a protein database
   if ( ($self->get_alphabet eq 'protein') && ($self->get_unidirectional != 1) ) {
      $self->throw("Got <unidirectional> = $unidirectional but can only use ".
         "<unidirectional> = 1 with proteic reference sequences\n");
   }

   return $self;
}


sub _init_db {
   # Read and import sequences
   # Parameters:
   #   * FASTA file containing the sequences or '-' for stdin. REQUIRED
   #   * Abundance file (optional): To avoid registering unwanted sequences
   #   * Delete chars (optional): Characters to delete from the sequences.
   #   * Minimum sequence size: Skip sequences smaller than that
   my ($self, $fasta_file, $abundance_file, $delete_chars, $min_len) = @_;

   # Get list of all IDs with a manually-specified abundance
   my %ids_to_keep;
   my $nof_ids_to_keep = 0;
   if ($abundance_file) {
      my ($ids) = community_read_abundances($abundance_file);
      for my $comm_num (0 .. $#$ids) {
         for my $gen_num ( 0 .. scalar @{$$ids[$comm_num]} - 1 ) {
            my $id = $$ids[$comm_num][$gen_num];
            $ids_to_keep{$id} = undef;
            $nof_ids_to_keep++;
         }
      }
   }

   # Index input file
   my $db = Bio::DB::Fasta->new($fasta_file, -reindex => 1);
   $self->_set_database($db);

   # List sequences that are ok to use
   my %seq_ids;
   my $nof_seqs;
   my %mol_types;
   
   my $stream = $db->get_PrimarySeq_stream;
   while (my $seq = $stream->next_seq) {

      # Skip empty sequences
      next if not $seq->seq;

      # Record molecule type
      $mol_types{$seq->alphabet}++;

      # Skip unwanted sequences
      my $seq_id = $seq->id;
      next if ($nof_ids_to_keep > 0) && (not exists $ids_to_keep{$seq_id});

      # Remove specified characters
      $seq = $self->_remove_chars($seq, $delete_chars);

      # Skip sequence if is not empty
      next if not defined $seq;

      # Skip the sequence if it is too small
      next if $seq->length < $min_len;

      # Record this sequence
      $seq_ids{$seq->id} = undef;
      $nof_seqs++;
   }

   # Error if no usable sequences in the database
   if ($nof_seqs == 0) {
      $self->throw("No genome sequences could be used. If you specified a file ".
         "of abundances for the genome sequences, make sure that their ID match".
         " the ID in the FASTA file. If you specified amplicon primers, verify ".
         "that they match some genome sequences.\n");
   }

   # Determine database type: dna, rna, protein
   my $db_alphabet = $self->_set_alphabet( $self->_get_mol_type(\%mol_types) );

   # Record the sequence IDs
   $self->_set_ids( \%seq_ids );

   return $db;
}


#sub get_primers {
#   my ($self) = @_;
#   return $self->{'primers'};
#}


#sub _set_primers {
#  my ($self, $forward_reverse_primers) = @_;
#  # Read primer file and convert primers into regular expressions to catch
#  # amplicons present in the database
#  if (defined $forward_reverse_primers) {

#    # Read primers from FASTA file
#    my $primer_in = Bio::SeqIO->newFh(
#      -file   => $forward_reverse_primers,
#      -format => 'fasta',
#    );

#    # Mandatory first primer
#    my $primer = <$primer_in>;
#    if (not defined $primer) {
#      $self->throw("The file '$forward_reverse_primers' contains no primers\n");
#    }
#    $primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences
#    $self->_set_forward_regexp( iupac_to_regexp($primer->seq) );
#    $primer = undef;

#    # Take reverse-complement of optional reverse primers
#    $primer = <$primer_in>;
#    if (defined $primer) {
#      $primer->alphabet('dna');
#      $primer = $primer->revcom;
#      $self->_set_reverse_regexp( iupac_to_regexp($primer->seq) );
#    }

#  }
#  $self->{'primers'} = $forward_reverse_primers;
#  return $self->get_primers;
#}


#sub get_forward_regexp {
#   my ($self) = @_;
#   return $self->{'forward_regexp'};
#}


#sub _set_forward_regexp {
#   my ($self, $val) = @_;
#   $self->{'forward_regexp'} = $val;
#   return $self->get_forward_regexp;
#}


#sub get_reverse_regexp {
#   my ($self) = @_;
#   return $self->{'reverse_regexp'};
#}


#sub _set_reverse_regexp {
#   my ($self, $val) = @_;
#   $self->{'reverse_regexp'} = $val;
#   return $self->get_reverse_regexp;
#}


sub get_alphabet {
   my ($self) = @_;
   return $self->{'alphabet'};
}


sub _set_alphabet {
   my ($self, $val) = @_;
   $self->{'alphabet'} = $val;
   return $self->get_alphabet;
}


sub get_ids {
   my ($self) = @_;
   my @ids = keys %{$self->{'ids'}};
   return \@ids;
}


sub _set_ids {
   my ($self, $val) = @_;
   $self->{'ids'} = $val;
   return $self->get_ids;
}


sub get_unidirectional {
   my ($self) = @_;
   return $self->{'unidirectional'};
}


sub _set_unidirectional {
   my ($self, $val) = @_;
   # Error if using wrong direction on protein database
   if ( ($self->get_alphabet eq 'protein') && ($val != 1) ) {
      $self->throw("Got <unidirectional> = $val but can only use ".
         "<unidirectional> = 1 with proteic reference sequences\n");
   }
   $self->{'unidirectional'} = $val;
   return $self->get_unidirectional;
}


sub get_minimum_length {
   my ($self) = @_;
   return $self->{'minimum_length'};
}


sub _set_minimum_length {
   my ($self, $val) = @_;
   $self->{'minimum_length'} = $val;
   return $self->get_minimum_length;
}


sub get_delete_chars {
   my ($self) = @_;
   return $self->{'delete_chars'};
}


sub _set_delete_chars {
   my ($self, $val) = @_;
   $self->{'delete_chars'} = $val;
   return $self->get_delete_chars;
}


sub get_database {
   my ($self) = @_;
   return $self->{'database'};
}


sub _set_database {
   my ($self, $val) = @_;
   $self->{'database'} = $val;
   return $self->get_database;
}


sub get_seq {
   my ($self, $id) = @_;
   # Get a sequence from the database. The query format is: id:start..end/strand
   # Only the id is mandatory. Start and end default to the full-length sequence
   # and strand defaults to 1.

   # Extract id, start, stop, and strand
   $id =~ s/\/(.+)$//i;
   my $strand = $1 || 1;
   ($id =~ s/:(\d+)..(\d+)$//i);
   my ($start, $stop) = ($1, $2);

   # Check that sequence is allowed
   if (not exists $self->{'ids'}->{$id}) {
      return undef;
   }

   # Invert start and stop for sequences on reverse strand
   if ($start && $stop && ($strand < 0) ) {
      ($start, $stop) = ($stop, $start);
   }

   #### if forbidden chars, start and stop provided, probably need to remove
   #### forbidden chars first

   # Get sequence from database
   my $seq = Bio::PrimarySeq->new(
      -id     => $id,
      -seq    => $self->{'database'}->seq($id, $start, $stop),
   );

   if ( ((not $start) || (not $stop)) && ($strand < 0) ) {
      $seq = $seq->revcom;
   }

   return $seq;
}


#sub next_seq {
#   my ($self) = @_;

#   # Get the database sequence stream, or set it the first time
#   my $stream = $self->get_stream ||
#                $self->_set_stream($self->get_database->get_PrimarySeq_stream);

#   my $seq = $stream->next_seq;

#   if (not defined $seq) {
#      # End of stream
#      return undef;
#   }

#   # If we are sequencing from the reverse strand, reverse complement now
#   if ($self->get_unidirectional == -1) {
#      $seq = $seq->revcom;
#   }

#   # then delete chars
#   #my $delete_chars   = $self->get_delete_chars;

#   # then fetch amplicons

#   # finally remove seqs < min_len

#   

#   # Extract amplicons if needed
##   my $amp_seqs;
##   if (defined $self->get_forward_regexp) {
##      $amp_seqs = $self->database_extract_amplicons($seq, $self->get_forward_regexp,
##      $self->get_reverse_regexp, \%ids_to_keep);
##      next if scalar @$amp_seqs == 0;
##   } else {
##      $amp_seqs = [$seq];
##   }

##   for my $amp_seq (@$amp_seqs) {
##      # Remove forbidden chars
##      if ( (defined $delete_chars) && (not $delete_chars eq '') ) {
##         my $clean_seq = $amp_seq->seq;
##         $clean_seq =~ s/[$delete_chars]//gi;
##         $amp_seq->seq($clean_seq);
##      }
##      # Skip the sequence if it is too small
##      next if $amp_seq->length < $min_len;
##      # Save amplicon sequence and identify them by their unique object reference
##      $seq_db{$amp_seq} = $amp_seq;
##      $seq_ids{$ref_seq_id}{$amp_seq} = undef;
##   }

#}


sub _remove_chars {
   # Remove forbidden chars
   my ($self, $seq, $chars) = @_;
   if ( defined($chars) && not($chars eq '') ) {

      my $seq_string = $seq->seq;
      my $count = ($seq_string =~ s/[$chars]//gi);

      if ( length $seq_string == 0 ) {
         # All characters were removed
         $seq = undef;
      } else {
         if ($count > 0) {
            # Some characters were removed.
            # Cannot modify a sequence from Bio::DB::Fasta. Create a new one if needed.
            $seq = Bio::PrimarySeq->new(
               -id  => $seq->id,
               -seq => $seq_string,
            );
         }
      }

   }

   return $seq;
}


sub _get_mol_type {
  # Given a count of the different molecule types in the database, determine
  # what molecule type it is.
  my ($self, $mol_types) = @_;
  my $max_count = 0;
  my $max_type  = '';

  while (my ($type, $count) = each %$mol_types) {
    if ($count > $max_count) {
      $max_count = $count;
      $max_type  = $type;
    }
  }
  my $other_count = 0;
  while (my ($type, $count) = each %$mol_types) {
    if (not $type eq $max_type) {
      $other_count += $count;
    }
  }
  if ($max_count < $other_count) {
    $self->throw("Cannot determine to what type of molecules the reference ".
        "sequences belong. Got $max_count sequences of type '$max_type' and ".
        "$other_count others.\n");
  }
  if ( (not $max_type eq 'dna') &&
       (not $max_type eq 'rna') &&
       (not $max_type eq 'protein') ) {
    $self->throw("Reference sequences have an unknown alphabet '$max_type'.\n");
  }
  return $max_type;
}


###sub _extract_amplicons {
###  my ($self, $seq, $forward_regexp, $reverse_regexp, $ids_to_keep) = @_;
###  # A database sequence can have several amplicons, e.g. a genome can have 
###  # several 16S rRNA genes. Extract all amplicons from a sequence (both strands)
###  # but take only the shortest when amplicons are nested.
###  # Fetch amplicons from both strands

###  # Get amplicons from forward and reverse strand
###  my $fwd_amplicons = _extract_amplicons_from_strand($seq, $forward_regexp, $reverse_regexp, 1);
###  my $rev_amplicons = _extract_amplicons_from_strand($seq, $forward_regexp, $reverse_regexp, -1);

###  # Deal with nested amplicons by removing the longest of the two
###  my $re = qr/(\d+)\.\.(\d+)/;
###  for (my $rev = 0; $rev < scalar @$rev_amplicons; $rev++) {
###    my ($rev_start, $rev_end) = ( $rev_amplicons->[$rev]->{_amplicon} =~ m/$re/ );
###    for (my $fwd = 0; $fwd < scalar @$fwd_amplicons; $fwd++) {
###      my ($fwd_start, $fwd_end) = ( $fwd_amplicons->[$fwd]->{_amplicon} =~ m/$re/ );
###      if ( ($fwd_start < $rev_start) && ($rev_end < $fwd_end) ) {
###        splice @$fwd_amplicons, $fwd, 1; # Remove forward amplicon
###        $fwd--;
###        next;
###      }
###      if ( ($rev_start < $fwd_start) && ($fwd_end < $rev_end) ) {
###        splice @$rev_amplicons, $rev, 1; # Remove reverse amplicon
###        $rev--;
###      }
###    }
###  }
###  
###  my $amplicons = [ @$fwd_amplicons, @$rev_amplicons ];

###  # Complain if primers did not match explicitly specified reference sequence
###  my $seqid = $seq->id;
###  if ( (scalar keys %{$ids_to_keep} > 0) &&
###       (exists $$ids_to_keep{$seqid}   ) &&
###       (scalar @$amplicons == 0         ) ) {
###    die "Error: Requested sequence $seqid did not match the specified forward primer.\n";
###  }

###  return $amplicons;
###}


###sub _extract_amplicons_from_strand {
###  # Get amplicons from the given strand (orientation) of the given sequence.
###  # For nested amplicons, only the shortest is returned to mimic PCR.
###  my ($seq, $forward_regexp, $reverse_regexp, $orientation) = @_;

###  # Reverse-complement sequence if looking at a -1 orientation
###  my $seqstr;
###  if ($orientation == 1) {
###    $seqstr = $seq->seq;
###  } elsif ($orientation == -1) {
###    $seqstr = $seq->revcom->seq;
###  } else {
###    die "Error: Invalid orientation '$orientation'\n";
###  }

###  # Get amplicons from sequence string
###  my $amplicons = [];
###  if ( (defined $forward_regexp) && (not defined $reverse_regexp) ) {
###    while ( $seqstr =~ m/($forward_regexp)/g ) {
###      my $start = pos($seqstr) - length($1) + 1;
###      my $end   = $seq->length;
###      push @$amplicons, _create_amplicon($seq, $start, $end, $orientation);
###    }
###  } elsif ( (defined $forward_regexp) && (defined $reverse_regexp) ) {
###    while ( $seqstr =~ m/($forward_regexp.*?$reverse_regexp)/g ) {
###      my $end   = pos($seqstr);
###      my $start = $end - length($1) + 1;
###      # Now trim the left end to obtain the shortest amplicon
###      my $ampliconstr = substr $seqstr, $start - 1, $end - $start + 1;
###      if ($ampliconstr =~ m/$forward_regexp.*($forward_regexp)/g) {
###         $start += pos($ampliconstr) - length($1);
###      }
###      push @$amplicons, _create_amplicon($seq, $start, $end, $orientation);
###    }
###  } else {
###    die "Error: Need to provide at least a forward primer\n";
###  }

###  return $amplicons;
###}


###sub _create_amplicon {
###  # Create an amplicon sequence and register its coordinates
###  my ($seq, $start, $end, $orientation) = @_;
###  my $amplicon;
###  my $coord;

###  if ($orientation == -1) {
###    # Calculate coordinates relative to forward strand. For example, given a
###    # read starting at 10 and ending at 23 on the reverse complement of a 100 bp
###    # sequence, return complement(77..90).
###    $amplicon = $seq->revcom->trunc($start, $end);
###    my $seq_len = $seq->length;
###    $start = $seq_len - $start + 1;
###    $end   = $seq_len - $end + 1;
###    ($start, $end) = ($end, $start);
###    $coord = "complement($start..$end)";
###  } else {
###    $amplicon = $seq->trunc($start, $end);
###    $coord = "$start..$end";
###  }
###  $amplicon->{_amplicon} = $coord;

###  return $amplicon
###}


####
#sub DESTROY {
# remove indexed files
#}
####

1;
