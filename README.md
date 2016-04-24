*This is a copy from
http://sourceforge.net/projects/biogrinder/files/biogrinder/Grinder-0.5.3/Grinder-0.5.3.tar.gz
with better formatted README, no source code has been changed.*

<!-- it's odd that in the table of contents that some words need to be
connected with dash while others with underscore in order to make linking work
when rendered on github page -->

1. [NAME](#name)
1. [DESCRIPTION](#description)
1. [CITATION](#citation)
1. [VERSION](#version)
1. [AUTHOR](#author)
1. [INSTALLATION](#installation)
   1. [Dependencies](#dependencies)
   1. [Procedure](#procedure)
   1. [No administrator privileges?](#no-administrator-privileges)
1. [RUNNING GRINDER](#running-grinder)
1. [REFERENCE SEQUENCE DATABASE](#reference-sequence-database)
1. [CLI EXAMPLES](#cli-examples)
1. [CLI REQUIRED ARGUMENTS](#cli-required-arguments)
1. [CLI OPTIONAL ARGUMENTS](#cli-optional-arguments)
1. [CLI OUTPUT](#cli-output)
1. [API EXAMPLES](#api-examples)
1. [API METHODS](#api-methods)
   1. [new](#new)
   1. [next\_lib](#next_lib)
   1. [next\_read](#next_read)
   1. [get\_random\_seed](#get_random_seed)
1. [COPYRIGHT](#copyright)
1. [BUGS](#bugs)

## Name

Grinder - A versatile omics shotgun and amplicon sequencing read
simulator

## Description

Grinder is a versatile program to create random shotgun and amplicon
sequence libraries based on DNA, RNA or proteic reference sequences
provided in a FASTA file.

Grinder can produce genomic, metagenomic, transcriptomic,
metatranscriptomic, proteomic, metaproteomic shotgun and amplicon
datasets from current sequencing technologies such as Sanger, 454,
Illumina. These simulated datasets can be used to test the accuracy of
bioinformatic tools under specific hypothesis, e.g. with or without
sequencing errors, or with low or high community diversity. Grinder may
also be used to help decide between alternative sequencing methods for a
sequence-based project, e.g. should the library be paired-end or not,
how many reads should be sequenced.

Grinder features include:

* shotgun or amplicon read libraries
* omics support to generate genomic, transcriptomic, proteomic, metagenomic,
  metatranscriptomic or metaproteomic datasets
* arbitrary read length distribution and number of reads
* simulation of PCR and sequencing errors (chimeras, point mutations,
  homopolymers)
* support for paired-end (mate pair) datasets
* specific rank-abundance settings or manually given abundance for each genome,
  gene or protein
* creation of datasets with a given richness (alpha diversity)
* independent datasets can share a variable number of genomes (beta diversity)
* modeling of the bias created by varying genome lengths or gene copy number
* profile mechanism to store preferred options
* available to biologists or power users through multiple interfaces: GUI, CLI
  and API

Briefly, given a FASTA file containing reference sequence (genomes,
genes, transcripts or proteins), Grinder performs the following steps:

1. Read the reference sequences, and for amplicon datasets, extracts full-length
   reference PCR amplicons using the provided degenerate PCR primers.

1. Determine the community structure based on the provided alpha diversity
   (number of reference sequences in the library), beta diversity (number of
   reference sequences in common between several independent libraries) and
   specified rank- abundance model.

1. Take shotgun reads from the reference sequences or amplicon reads from the
   full- length reference PCR amplicons. The reads may be paired-end reads when
   an insert size distribution is specified. The length of the reads depends on
   the provided read length distribution and their abundance depends on the
   relative abundance in the community structure. Genome length may also biases
   the number of reads to take for shotgun datasets at this step. Similarly, for
   amplicon datasets, the number of copies of the target gene in the reference
   genomes may bias the number of reads to take.

1. Alter reads by inserting sequencing errors (indels, substitutions and
   homopolymer errors) following a position-specific model to simulate reads
   created by current sequencing technologies (Sanger, 454, Illumina). Write the
   reads and their quality scores in FASTA, QUAL and FASTQ files.

## Citation

If you use Grinder in your research, please cite:

```
Angly FE, Willner D, Rohwer F, Hugenholtz P, Tyson GW (2012), Grinder: a versatile amplicon and shotgun sequence simulator, Nucleic Acids Reseach
```
Available from <http://dx.doi.org/10.1093/nar/gks251>.

## Author

Florent Angly \<<florent.angly@gmail.com>\>

## Installation

You need to install these dependencies first:

* [Perl](http://www.perl.com/download.csp) (>= 5.6)
* [make](http://www.gnu.org/s/make/)
* The following CPAN Perl modules are dependencies that will be installed
  automatically for you:
  * `Bioperl` modules (>=1.6.901). Note that some unreleased Bioperl modules
    have been included in Grinder.
  * `Getopt::Euclid` (>= 0.3.4)
  * `List::Util`. First released with Perl v5.7.3
  * `Math::Random::MT` (>= 1.13)
  * `version` (>= 0.77). First released with Perl v5.9.0

To install Grinder globally on your system, run the following commands
in a terminal or command prompt:

##### On Linux, Unix, MacOS:

```
perl Makefile.PL
make
# The following step needs administrator privileges
make install
```

##### On Windows:

```
perl Makefile.PL
nmake
# The following step needs administrator privileges
nmake install
```

If you don't have administrator privileges, then Grinder needs to be installed
in your home directory.

Follow the instructions to install `local::lib`
[here](http://search.cpan.org/~apeiron/local-lib-1.008004/lib/local/lib.pm#The_bootstrapping_technique).
Then, every Perl module that you install manually or through the CPAN
command-line application will be installed in your home directory. At last,
install Grinder by following the instructions detailed as above.

## Running Grinder

You can run Grinder using the

* command-line interface (CLI). see ```grinder --help```
* programming interface (API). see ```perldoc Grinder```
* or graphical user interface (GUI) in Galaxy. See the Galaxy documentation
  [here](http://wiki.g2.bx.psu.edu/FrontPage)


##### The `utils` folder

It's included in the Grinder package and contains the following utilities:

* `average genome size`: This calculates the average genome size (in bp) of a
    simulated random library produced by Grinder.

* `change_paired_read_orientation`: This reverses the orientation of each
    second mate-pair read (ID ending in /2) in a FASTA file.


## Reference Sequence Database

A variety of FASTA databases can be used as input for Grinder. For example,
[GreenGenes database](http://greengenes.lbl.gov/Download),
[RDP](http://rdp.cme.msu.edu) or [Silva](http://www.arb-silva.de). They all
provide many 16S rRNA sequences, but Silva includes eukaryotic sequences, as
well.

While 16S rRNA is a popular gene, datasets containing any type of gene could can
used in the same fashion to generate simulated amplicon datasets, provided that
appropriate primers are used.

The over 2,400 curated microbial genome sequences in the
[NCBI RefSeq microbial collection](ftp://ftp.ncbi.nih.gov/refseq/release/microbial/)
would also be suitable for producing 16S rRNA simulated datasets (using the
adequate primers). However, the lower diversity of this database compared to the
previous two makes it more appropriate for producing artificial microbial
metagenomes. Individual genomes from this database are also very suitable for
the simulation of single or double-barreled shotgun libraries. Similarly, the
[RefSeq viral collection](ftp://ftp.ncbi.nih.gov/refseq/release/viral/) contains
over 3,100 curated viral sequences which can be used to produce artificial viral
metagenomes.

Quite a few eukaryotic organisms have been sequenced and their genome or
genes can be the basis for simulating genomic, transcriptomic (RNA-seq)
or proteomic datasets. For example, you can use 
* the human genome available at ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/
* the human transcripts at ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz
* or the human proteome at ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz

## CLI Examples

A shotgun DNA library with a coverage of 0.1X

```
grinder -reference_file genomes.fna -coverage_fold 0.1
```

Same thing but save the result files in a specific folder and with a specific
name

```
grinder -reference_file genomes.fna -coverage_fold 0.1 -base_name my_name -output_dir my_dir
```

A DNA shotgun library with 1000 reads

```
grinder -reference_file genomes.fna -total_reads 1000
```

A DNA shotgun library where species are distributed according to a power law

```
grinder -reference_file genomes.fna -abundance_model powerlaw 0.1
```

A DNA shotgun library with 123 genomes taken random from the given genomes

```
grinder -reference_file genomes.fna -diversity 123
```

Two DNA shotgun libraries that have 50% of the species in common

```
grinder -reference_file genomes.fna -num_libraries 2 -shared_perc 50
```

Two DNA shotgun library with no species in common and distributed according to a
exponential rank-abundance model. Note that because the parameter value for the
exponential model is omitted, each library uses a different randomly chosen
value:

```
grinder -reference_file genomes.fna -num_libraries 2 -abundance_model exponential
```

A DNA shotgun library where species relative abundances are manually specified

```
grinder -reference_file genomes.fna -abundance_file my_abundances.txt
```

A DNA shotgun library with Sanger reads

```
grinder -reference_file genomes.fna -read_dist 800 -mutation_dist linear 1 2 -mutation_ratio 80 20
```

A DNA shotgun library with first-generation 454 reads

```
grinder -reference_file genomes.fna -read_dist 100 normal 10 -homopolymer_dist balzer
```

A paired-end DNA shotgun library, where the insert size is normally distributed
around 2.5 kbp and has 0.2 kbp standard deviation

```
grinder -reference_file genomes.fna -insert_dist 2500 normal 200
```

A transcriptomic dataset

```
grinder -reference_file transcripts.fna
```

A unidirectional transcriptomic dataset

```
grinder -reference_file transcripts.fna -unidirectional 1
```

**Note:** the use of `-unidirectional 1` to prevent reads to be taken from the
reverse-complement of the reference sequences.

A proteomic dataset

```
grinder -reference_file proteins.faa -unidirectional 1
```

A 16S rRNA amplicon library
```
grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1
```

**Note:** the use of `-length_bias 0` because reference sequence length should
not affect the relative abundance of amplicons.

The same amplicon library with 20% of chimeric reads (90% bimera, 10% trimera)

```
grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -chimera_perc 20 -chimera_dist 90 10
```

Three 16S rRNA amplicon libraries with specified MIDs and no reference sequences
in common

```
grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -num_libraries 3 -multiplex_ids MIDs.fna
```

Reading reference sequences from the standard input, which allows you to
decompress FASTA files on the fly:

```
zcat microbial_db.fna.gz | grinder -reference_file - -total_reads 100
```


## CLI Output

For each shotgun or amplicon read library requested, the following files are
generated:

1. A rank-abundance file, tab-delimited, that shows the relative abundance of
   the different reference sequences
1. A file containing the read sequences in FASTA format. The read headers
   contain information necessary to track from which reference sequence each
   read was taken and what errors it contains. This file is not generated if
   `fastq_output` option was provided.
1. If the `qual_levels` option was specified, a file containing the quality
   scores of the reads (in QUAL format).
1. If the `fastq_output`` option was provided, a file containing the read
   sequences in FASTQ format.


## API Examples

The Grinder API allows to conveniently use Grinder within Perl scripts. Here is
a synopsis:

```
use Grinder;

# Set up a new factory (see the OPTIONS section for a complete list of parameters)
my $factory = Grinder->new( -reference_file => 'genomes.fna' );

# Process all shotgun libraries requested
while ( my $struct = $factory->next_lib ) {

  # The ID and abundance of the 3rd most abundant genome in this community
  my $id = $struct->{ids}->[2];
  my $ab = $struct->{abs}->[2];

  # Create shotgun reads
  while ( my $read = $factory->next_read) {

    # The read is a Bioperl sequence object with these properties:
    my $read_id     = $read->id;     # read ID given by Grinder
    my $read_seq    = $read->seq;    # nucleotide sequence
    my $read_mid    = $read->mid;    # MID or tag attached to the read
    my $read_errors = $read->errors; # errors that the read contains

    # Where was the read taken from? The reference sequence refers to the
    # database sequence for shotgun libraries, amplicon obtained from the
    # database sequence, or could even be a chimeric sequence
    my $ref_id     = $read->reference->id; # ID of the reference sequence
    my $ref_start  = $read->start;         # start of the read on the reference
    my $ref_end    = $read->end;           # end of the read on the reference
    my $ref_strand = $read->strand;        # strand of the reference

  }
}

# Similarly, for shotgun mate pairs
my $factory = Grinder->new( -reference_file => 'genomes.fna',
                            -insert_dist    => 250            );
while ( $factory->next_lib ) {
  while ( my $read = $factory->next_read ) {
    # The first read is the first mate of the mate pair
    # The second read is the second mate of the mate pair
    # The third read is the first mate of the next mate pair
    # ...
  }
}

# To generate an amplicon library
my $factory = Grinder->new( -reference_file  => 'genomes.fna',
                            -forward_reverse => '16Sgenes.fna',
                            -length_bias     => 0,
                            -unidirectional  => 1              );
while ( $factory->next_lib ) {
  while ( my $read = $factory->next_read) {
    # ...
  }
}
```

API METHODS
===========

The rest of the documentation details the available Grinder API methods.

new
---

Title : new

Function: Create a new Grinder factory initialized with the passed
arguments. Available parameters described in the OPTIONS section.

Usage : my \$factory = Grinder-\>new( -reference\_file =\> 'genomes.fna'
);

Returns : a new Grinder object

next\_lib
---------

Title : next\_lib

Function: Go to the next shotgun library to process.

Usage : my \$struct = \$factory-\>next\_lib;

Returns : Community structure to be used for this library, where
\$struct-\>{ids} is an array reference containing the IDs of the genome
making up the community (sorted by decreasing relative abundance) and
\$struct-\>{abs} is an array reference of the genome abundances (in the
same order as the IDs).

next\_read
----------

Title : next\_read

Function: Create an amplicon or shotgun read for the current library.

Usage : my \$read = \$factory-\>next\_read; \# for single read my
\$mate1 = \$factory-\>next\_read; \# for mate pairs my \$mate2 =
\$factory-\>next\_read;

Returns : A sequence represented as a Bio::Seq::SimulatedRead object

get\_random\_seed
-----------------

Title : get\_random\_seed

Function: Return the number used to seed the pseudo-random number
generator

Usage : my \$seed = \$factory-\>get\_random\_seed;

Returns : seed number

------------------------------------------------------------------------

COPYRIGHT
=========

Copyright 2009-2012 Florent ANGLY \<<florent.angly@gmail.com>\>

Grinder is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version. Grinder is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. You should have received a
copy of the GNU General Public License along with Grinder. If not, see
\<http://www.gnu.org/licenses/\>.

------------------------------------------------------------------------

BUGS
====

All complex software has bugs lurking in it, and this program is no
exception. If you find a bug, please report it on the SourceForge
Tracker for Grinder:
[http://sourceforge.net/tracker/](http://sourceforge.net/tracker/?group_id=244196&atid=1124737)

Bug reports, suggestions and patches are welcome. Grinder's code is
developed on Sourceforge
([http://sourceforge.net/scm/](http://sourceforge.net/scm/?type=git&group_id=244196))
and is under Git revision control. To get started with a patch, do:

       git clone git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder
