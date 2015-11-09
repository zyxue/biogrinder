# This file is part of the Grinder package, copyright 2009-2013
# Florent Angly <florent.angly@gmail.com>, under the GPLv3 license

package Grinder;

use 5.006;
use strict;
use warnings;
use File::Spec;
use List::Util qw(max);
use Bio::SeqIO;
use Grinder::KmerCollection;
use Bio::Location::Split;
use Bio::Seq::SimulatedRead;
use Bio::SeqFeature::SubSeq;
use Bio::Tools::AmpliconSearch;
use Math::Random::MT qw(srand rand);
use Getopt::Euclid qw(:minimal_keys :defer);

use version; our $VERSION = version->declare('0.5.3');


#---------- GRINDER POD DOC ---------------------------------------------------#

=head1 NAME

Grinder - A versatile omics shotgun and amplicon sequencing read simulator

=head1 DESCRIPTION

Grinder is a versatile program to create random shotgun and amplicon sequence
libraries based on DNA, RNA or proteic reference sequences provided in a FASTA
file.

Grinder can produce genomic, metagenomic, transcriptomic, metatranscriptomic,
proteomic, metaproteomic shotgun and amplicon datasets from current sequencing
technologies such as Sanger, 454, Illumina. These simulated datasets can be used
to test the accuracy of bioinformatic tools under specific hypothesis, e.g. with
or without sequencing errors, or with low or high community diversity. Grinder
may also be used to help decide between alternative sequencing methods for a
sequence-based project, e.g. should the library be paired-end or not, how many
reads should be sequenced.

Grinder features include:

=over

=item *

shotgun or amplicon read libraries

=item *

omics support to generate genomic, transcriptomic, proteomic,
metagenomic, metatranscriptomic or metaproteomic datasets

=item *

arbitrary read length distribution and number of reads

=item *

simulation of PCR and sequencing errors (chimeras, point mutations, homopolymers)

=item *

support for paired-end (mate pair) datasets

=item *

specific rank-abundance settings or manually given abundance for each genome, gene or protein

=item *

creation of datasets with a given richness (alpha diversity)

=item *

independent datasets can share a variable number of genomes (beta diversity)

=item *

modeling of the bias created by varying genome lengths or gene copy number

=item *

profile mechanism to store preferred options

=item *

available to biologists or power users through multiple interfaces: GUI, CLI and API

=back

Briefly, given a FASTA file containing reference sequence (genomes, genes,
transcripts or proteins), Grinder performs the following steps:

=over

=item 1.

Read the reference sequences, and for amplicon datasets, extracts full-length
reference PCR amplicons using the provided degenerate PCR primers.

=item 2.

Determine the community structure based on the provided alpha diversity (number
of reference sequences in the library), beta diversity (number of reference
sequences in common between several independent libraries) and specified rank-
abundance model.

=item 3.

Take shotgun reads from the reference sequences or amplicon reads from the full-
length reference PCR amplicons. The reads may be paired-end reads when an insert
size distribution is specified. The length of the reads depends on the provided
read length distribution and their abundance depends on the relative abundance
in the community structure. Genome length may also biases the number of reads to
take for shotgun datasets at this step. Similarly, for amplicon datasets, the
number of copies of the target gene in the reference genomes may bias the number
of reads to take.

=item 4.

Alter reads by inserting sequencing errors (indels, substitutions and homopolymer
errors) following a position-specific model to simulate reads created by current
sequencing technologies (Sanger, 454, Illumina). Write the reads and their
quality scores in FASTA, QUAL and FASTQ files.

=back

=head1 CITATION

If you use Grinder in your research, please cite:

   Angly FE, Willner D, Rohwer F, Hugenholtz P, Tyson GW (2012), Grinder: a
   versatile amplicon and shotgun sequence simulator, Nucleic Acids Reseach

Available from L<http://dx.doi.org/10.1093/nar/gks251>.

=head1 VERSION

0.5.3

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 INSTALLATION

=head2 Dependencies

You need to install these dependencies first:

=over

=item *

Perl (>= 5.6)

L<http://www.perl.com/download.csp>

=item *

make

Many systems have make installed by default. If your system does not, you should
install the implementation of make of your choice, e.g. GNU make: L<http://www.gnu.org/s/make/>

=back

The following CPAN Perl modules are dependencies that will be installed automatically
for you:

=over

=item *

Bioperl modules (>=1.6.901).

Note that some unreleased Bioperl modules have been included in Grinder.

=item *

Getopt::Euclid (>= 0.3.4)

=item *

List::Util

First released with Perl v5.7.3

=item *

Math::Random::MT (>= 1.13)

=item *

version (>= 0.77)

First released with Perl v5.9.0

=back

=head2 Procedure

To install Grinder globally on your system, run the following commands in a
terminal or command prompt:

On Linux, Unix, MacOS:

   perl Makefile.PL
   make

And finally, with administrator privileges:

   make install

On Windows, run the same commands but with nmake instead of make.

=head2 No administrator privileges?

If you do not have administrator privileges, Grinder needs to be installed in
your home directory.

First, follow the instructions to install local::lib
at L<http://search.cpan.org/~apeiron/local-lib-1.008004/lib/local/lib.pm#The_bootstrapping_technique>. After local::lib is installed, every Perl
module that you install manually or through the CPAN command-line application
will be installed in your home directory.

Then, install Grinder by following the instructions detailed in the "Procedure"
section.

=head1 RUNNING GRINDER

After installation, you can run Grinder using a command-line interface (CLI), 
an application programming interface (API) or a graphical user interface (GUI)
in Galaxy.

To get the usage of the CLI, type:

  grinder --help

More information, including the documentation of the Grinder API, which allows
you to run Grinder from within other Perl programs, is available by typing:

  perldoc Grinder

To run the GUI, refer to the Galaxy documentation at L<http://wiki.g2.bx.psu.edu/FrontPage>.

The 'utils' folder included in the Grinder package contains some utilities:

=over

=item average genome size:

This calculates the average genome size (in bp) of a simulated random library
produced by Grinder.

=item change_paired_read_orientation:

This reverses the orientation of each second mate-pair read (ID ending in /2)
in a FASTA file.

=back

=head1 REFERENCE SEQUENCE DATABASE

A variety of FASTA databases can be used as input for Grinder. For example, the
GreenGenes database (L<http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/Isolated_named_strains_16S_aligned.fasta>)
contains over 180,000 16S rRNA clone sequences from various species which would
be appropriate to produce a 16S rRNA amplicon dataset. A set of over 41,000 OTU
representative sequences and their affiliation in seven different taxonomic
sytems can also be used for the same purpose (L<http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta>
and L<http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010/taxonomies/>). The
RDP (L<http://rdp.cme.msu.edu/download/release10_27_unaligned.fa.gz>) and Silva
(L<http://www.arb-silva.de/no_cache/download/archive/release_108/Exports/>)
databases also provide many 16S rRNA sequences and Silva includes eukaryotic
sequences. While 16S rRNA is a popular gene, datasets containing any type of gene
could be used in the same fashion to generate simulated amplicon datasets, provided
appropriate primers are used.

The >2,400 curated microbial genome sequences in the NCBI RefSeq collection
(L<ftp://ftp.ncbi.nih.gov/refseq/release/microbial/>) would also be suitable for
producing 16S rRNA simulated datasets (using the adequate primers). However, the
lower diversity of this database compared to the previous two makes it more
appropriate for producing artificial microbial metagenomes. Individual genomes
from this database are also very suitable for the simulation of single or
double-barreled shotgun libraries. Similarly, the RefSeq database contains
over 3,100 curated viral sequences (L<ftp://ftp.ncbi.nih.gov/refseq/release/viral/>)
which can be used to produce artificial viral metagenomes.

Quite a few eukaryotic organisms have been sequenced and their genome or genes
can be the basis for simulating genomic, transcriptomic (RNA-seq) or proteomic 
datasets. For example, you can use the human genome available at
L<ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/>, the human transcripts
downloadable from L<ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz>
or the human proteome at L<ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz>.

=head1 CLI EXAMPLES

Here are a few examples that illustrate the use of Grinder in a terminal:

=over

=item 1.

A shotgun DNA library with a coverage of 0.1X

   grinder -reference_file genomes.fna -coverage_fold 0.1

=item 2.

Same thing but save the result files in a specific folder and with a specific name

   grinder -reference_file genomes.fna -coverage_fold 0.1 -base_name my_name -output_dir my_dir

=item 3.

A DNA shotgun library with 1000 reads

   grinder -reference_file genomes.fna -total_reads 1000

=item 4.

A DNA shotgun library where species are distributed according to a power law

   grinder -reference_file genomes.fna -abundance_model powerlaw 0.1

=item 5.

A DNA shotgun library with 123 genomes taken random from the given genomes

   grinder -reference_file genomes.fna -diversity 123

=item 6.

Two DNA shotgun libraries that have 50% of the species in common

   grinder -reference_file genomes.fna -num_libraries 2 -shared_perc 50

=item 7.

Two DNA shotgun library with no species in common and distributed according to a
exponential rank-abundance model. Note that because the parameter value for the
exponential model is omitted, each library uses a different randomly chosen value:

   grinder -reference_file genomes.fna -num_libraries 2 -abundance_model exponential

=item 8.

A DNA shotgun library where species relative abundances are manually specified

   grinder -reference_file genomes.fna -abundance_file my_abundances.txt

=item 9.

A DNA shotgun library with Sanger reads

   grinder -reference_file genomes.fna -read_dist 800 -mutation_dist linear 1 2 -mutation_ratio 80 20

=item 10.

A DNA shotgun library with first-generation 454 reads

   grinder -reference_file genomes.fna -read_dist 100 normal 10 -homopolymer_dist balzer

=item 11.

A paired-end DNA shotgun library, where the insert size is normally distributed
around 2.5 kbp and has 0.2 kbp standard deviation

   grinder -reference_file genomes.fna -insert_dist 2500 normal 200

=item 12.

A transcriptomic dataset

   grinder -reference_file transcripts.fna

=item 13.

A unidirectional transcriptomic dataset

   grinder -reference_file transcripts.fna -unidirectional 1

Note the use of -unidirectional 1 to prevent reads to be taken from the reverse-
complement of the reference sequences.

=item 14.

A proteomic dataset

   grinder -reference_file proteins.faa -unidirectional 1

=item 15.

A 16S rRNA amplicon library

   grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1

Note the use of -length_bias 0 because reference sequence length should not affect
the relative abundance of amplicons.

=item 16.

The same amplicon library with 20% of chimeric reads (90% bimera, 10% trimera)

   grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -chimera_perc 20 -chimera_dist 90 10

=item 17.

Three 16S rRNA amplicon libraries with specified MIDs and no reference sequences
in common

   grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -num_libraries 3 -multiplex_ids MIDs.fna

=item 18.

Reading reference sequences from the standard input, which allows you to
decompress FASTA files on the fly:

   zcat microbial_db.fna.gz | grinder -reference_file - -total_reads 100

=back

=head1 CLI REQUIRED ARGUMENTS

=over

=item -rf <reference_file> | -reference_file <reference_file> | -gf <reference_file> | -genome_file <reference_file>

FASTA file that contains the input reference sequences (full genomes, 16S rRNA
genes, transcripts, proteins...) or '-' to read them from the standard input. See the
README file for examples of databases you can use and where to get them from. 
Default: reference_file.default

=for Euclid:
   reference_file.type: readable
   reference_file.default: '-'

=back

=head1 CLI OPTIONAL ARGUMENTS

Basic parameters

=over

=item -tr <total_reads> | -total_reads <total_reads>

Number of shotgun or amplicon reads to generate for each library. Do not specify
this if you specify the fold coverage. Default: total_reads.default

=for Euclid:
   total_reads.type: +integer
   total_reads.default: 100

=item -cf <coverage_fold> | -coverage_fold <coverage_fold>

Desired fold coverage of the input reference sequences (the output FASTA length
divided by the input FASTA length). Do not specify this if you specify the number
of reads directly.

=for Euclid:
   coverage_fold.type: +number
   coverage_fold.excludes: total_reads

=back

Advanced shotgun and amplicon parameters

=over

=item -rd <read_dist>... | -read_dist <read_dist>...

Desired shotgun or amplicon read length distribution specified as:
   average length, distribution ('uniform' or 'normal') and standard deviation.

Only the first element is required. Examples:

  All reads exactly 101 bp long (Illumina GA 2x): 101
  Uniform read distribution around 100+-10 bp: 100 uniform 10
  Reads normally distributed with an average of 800 and a standard deviation of 100
    bp (Sanger reads): 800 normal 100
  Reads normally distributed with an average of 450 and a standard deviation of 50
    bp (454 GS-FLX Ti): 450 normal 50

Reference sequences smaller than the specified read length are not used. Default:
read_dist.default

=for Euclid:
   read_dist.type: string
   read_dist.default: [100]

=item -id <insert_dist>... | -insert_dist <insert_dist>...

Create paired-end or mate-pair reads spanning the given insert length.
Important: the insert is defined in the biological sense, i.e. its length includes
the length of both reads and of the stretch of DNA between them:
   0 : off,
   or: insert size distribution in bp, in the same format as the read length
       distribution (a typical value is 2,500 bp for mate pairs)
Two distinct reads are generated whether or not the mate pair overlaps. Default:
insert_dist.default

=for Euclid:
   insert_dist.type: string
   insert_dist.default: [0]

=item -mo <mate_orientation> | -mate_orientation <mate_orientation>

When generating paired-end or mate-pair reads (see <insert_dist>), specify the
orientation of the reads (F: forward, R: reverse):

   FR:  ---> <---  e.g. Sanger, Illumina paired-end, IonTorrent mate-pair
   FF:  ---> --->  e.g. 454
   RF:  <--- --->  e.g. Illumina mate-pair
   RR:  <--- <---

Default: mate_orientation.default

=for Euclid:
   mate_orientation.type: string, mate_orientation eq 'FF' || mate_orientation eq 'FR' || mate_orientation eq 'RF' || mate_orientation eq 'RR'
   mate_orientation.type.error: <mate_orientation> must be FR, FF, RF or RR (not mate_orientation)
   mate_orientation.default: 'FR'

=item -ec <exclude_chars> | -exclude_chars <exclude_chars>

Do not create reads containing any of the specified characters (case insensitive).
For example, use 'NX' to prevent reads with ambiguities (N or X). Grinder will
error if it fails to find a suitable read (or pair of reads) after 10 attempts.
Consider using <delete_chars>, which may be more appropriate for your case.
Default: 'exclude_chars.default'

=for Euclid:
   exclude_chars.type: string
   exclude_chars.default: ''

=item -dc <delete_chars> | -delete_chars <delete_chars>

Remove the specified characters from the reference sequences (case-insensitive),
e.g. '-~*' to remove gaps (- or ~) or terminator (*). Removing these characters
is done once, when reading the reference sequences, prior to taking reads. Hence
it is more efficient than <exclude_chars>. Default: delete_chars.default

=for Euclid:
   delete_chars.type: string
   delete_chars.default: ''

=item -fr <forward_reverse> | -forward_reverse <forward_reverse>

Use DNA amplicon sequencing using a forward and reverse PCR primer sequence
provided in a FASTA file. The reference sequences and their reverse complement
will be searched for PCR primer matches. The primer sequences should use the
IUPAC convention for degenerate residues and the reference sequences that that
do not match the specified primers are excluded. If your reference sequences are
full genomes, it is recommended to use <copy_bias> = 1 and <length_bias> = 0 to
generate amplicon reads. To sequence from the forward strand, set <unidirectional>
to 1 and put the forward primer first and reverse primer second in the FASTA
file. To sequence from the reverse strand, invert the primers in the FASTA file
and use <unidirectional> = -1. The second primer sequence in the FASTA file is
always optional. Example: AAACTYAAAKGAATTGRCGG and ACGGGCGGTGTGTRC for the 926F
and 1392R primers that target the V6 to V9 region of the 16S rRNA gene.

=for Euclid:
   forward_reverse.type: readable

=item -un <unidirectional> | -unidirectional <unidirectional>

Instead of producing reads bidirectionally, from the reference strand and its
reverse complement, proceed unidirectionally, from one strand only (forward or
reverse). Values: 0 (off, i.e. bidirectional), 1 (forward), -1 (reverse). Use
<unidirectional> = 1 for amplicon and strand-specific transcriptomic or
proteomic datasets. Default: unidirectional.default

=for Euclid:
   unidirectional.type: integer, unidirectional >= -1 && unidirectional <= 1
   unidirectional.type.error: <unidirectional> must be 0, 1 or -1 (not unidirectional)
   unidirectional.default: 0

=item -lb <length_bias> | -length_bias <length_bias>

In shotgun libraries, sample reference sequences proportionally to their length.
For example, in simulated microbial datasets, this means that at the same
relative abundance, larger genomes contribute more reads than smaller genomes
(and all genomes have the same fold coverage).
0 = no, 1 = yes. Default: length_bias.default

=for Euclid:
   length_bias.type: integer, length_bias == 0 || length_bias == 1
   length_bias.type.error: <length_bias> must be 0 or 1 (not length_bias)
   length_bias.default: 1

=item -cb <copy_bias> | -copy_bias <copy_bias>

In amplicon libraries where full genomes are used as input, sample species
proportionally to the number of copies of the target gene: at equal relative
abundance, genomes that have multiple copies of the target gene contribute more
amplicon reads than genomes that have a single copy. 0 = no, 1 = yes. Default:
copy_bias.default

=for Euclid:
   copy_bias.type: integer, copy_bias == 0 || copy_bias == 1
   copy_bias.type.error: <copy_bias> must be 0 or 1 (not copy_bias)
   copy_bias.default: 1

=back

Aberrations and sequencing errors

=over

=item -md <mutation_dist>... | -mutation_dist <mutation_dist>...

Introduce sequencing errors in the reads, under the form of mutations
(substitutions, insertions and deletions) at positions that follow a specified
distribution (with replacement): model (uniform, linear, poly4), model parameters.
For example, for a uniform 0.1% error rate, use: uniform 0.1. To simulate Sanger
errors, use a linear model where the errror rate is 1% at the 5' end of reads and
2% at the 3' end: linear 1 2. To model Illumina errors using the 4th degree
polynome 3e-3 + 3.3e-8 * i^4 (Korbel et al 2009), use: poly4 3e-3 3.3e-8.
Use the <mutation_ratio> option to alter how many of these mutations are
substitutions or indels. Default: mutation_dist.default

=for Euclid:
   mutation_dist.type: string
   mutation_dist.default: ['uniform', 0, 0]

=item -mr <mutation_ratio>... | -mutation_ratio <mutation_ratio>...

Indicate the percentage of substitutions and the number of indels (insertions
and deletions). For example, use '80 20' (4 substitutions for each indel) for
Sanger reads. Note that this parameter has no effect unless you specify the
<mutation_dist> option. Default: mutation_ratio.default

=for Euclid:
   mutation_ratio.type: num, mutation_ratio >= 0
   mutation_ratio.default: [80, 20]

=item -hd <homopolymer_dist> | -homopolymer_dist <homopolymer_dist>

Introduce sequencing errors in the reads under the form of homopolymeric
stretches (e.g. AAA, CCCCC) using a specified model where the homopolymer length
follows a normal distribution N(mean, standard deviation) that is function of
the homopolymer length n:

  Margulies: N(n, 0.15 * n)              ,  Margulies et al. 2005.
  Richter  : N(n, 0.15 * sqrt(n))        ,  Richter et al. 2008.
  Balzer   : N(n, 0.03494 + n * 0.06856) ,  Balzer et al. 2010.

Default: homopolymer_dist.default

=for Euclid:
   homopolymer_dist.type: string
   homopolymer_dist.default: 0

=item -cp <chimera_perc> | -chimera_perc <chimera_perc>

Specify the percent of reads in amplicon libraries that should be chimeric
sequences. The 'reference' field in the description of chimeric reads will
contain the ID of all the reference sequences forming the chimeric template.
A typical value is 10% for amplicons. This option can be used to generate
chimeric shotgun reads as well. Default: chimera_perc.default %

=for Euclid:
   chimera_perc.type: number, chimera_perc >= 0 && chimera_perc <= 100
   chimera_perc.type.error: <chimera_perc> must be a number between 0 and 100 (not chimera_perc)
   chimera_perc.default: 0

=item -cd <chimera_dist>... | -chimera_dist <chimera_dist>...

Specify the distribution of chimeras: bimeras, trimeras, quadrameras and
multimeras of higher order. The default is the average values from Quince et al.
2011: '314 38 1', which corresponds to 89% of bimeras, 11% of trimeras and 0.3%
of quadrameras. Note that this option only takes effect when you request the
generation of chimeras with the <chimera_perc> option. Default: chimera_dist.default

=for Euclid:
   chimera_dist.type: number, chimera_dist >= 0
   chimera_dist.type.error: <chimera_dist> must be a positive number (not chimera_dist)
   chimera_dist.default: [314, 38, 1]

=item -ck <chimera_kmer> | -chimera_kmer <chimera_kmer>

Activate a method to form chimeras by picking breakpoints at places where k-mers
are shared between sequences. <chimera_kmer> represents k, the length of the
k-mers (in bp). The longer the kmer, the more similar the sequences have to be
to be eligible to form chimeras. The more frequent a k-mer is in the pool of
reference sequences (taking into account their relative abundance), the more
often this k-mer will be chosen. For example, CHSIM (Edgar et al. 2011) uses this
method with a k-mer length of 10 bp. If you do not want to use k-mer information
to form chimeras, use 0, which will result in the reference sequences and
breakpoints to be taken randomly on the "aligned" reference sequences. Note that
this option only takes effect when you request the generation of chimeras with
the <chimera_perc> option. Also, this options is quite memory intensive, so you
should probably limit yourself to a relatively small number of reference sequences
if you want to use it. Default: chimera_kmer.default bp

=for Euclid:
   chimera_kmer.type: number, chimera_kmer == 0 || chimera_kmer >= 2
   chimera_kmer.type.error: <chimera_kmer> must be 0 or an integer larger than 1 (not chimera_kmer)
   chimera_kmer.default: 10

=back

Community structure and diversity

=over

=item -af <abundance_file> | -abundance_file <abundance_file>

Specify the relative abundance of the reference sequences manually in an input
file. Each line of the file should contain a sequence name and its relative
abundance (%), e.g. 'seqABC 82.1' or 'seqABC 82.1 10.2' if you are specifying two
different libraries.

=for Euclid:
   abundance_file.type: readable

=item -am <abundance_model>... | -abundance_model <abundance_model>...

Relative abundance model for the input reference sequences: uniform, linear, powerlaw,
logarithmic or exponential. The uniform and linear models do not require a
parameter, but the other models take a parameter in the range [0, infinity). If
this parameter is not specified, then it is randomly chosen. Examples:

  uniform distribution: uniform
  powerlaw distribution with parameter 0.1: powerlaw 0.1
  exponential distribution with automatically chosen parameter: exponential

Default: abundance_model.default

=for Euclid:
   abundance_model.type: string
   abundance_model.default: ['uniform', 1]

=item -nl <num_libraries> | -num_libraries <num_libraries>

Number of independent libraries to create. Specify how diverse and similar they
should be with <diversity>, <shared_perc> and <permuted_perc>. Assign them
different MID tags with <multiplex_mids>. Default: num_libraries.default

=for Euclid:
   num_libraries.type: +integer
   num_libraries.default: 1

=item -mi <multiplex_ids> | -multiplex_ids <multiplex_ids>

Specify an optional FASTA file that contains multiplex sequence identifiers
(a.k.a MIDs or barcodes) to add to the sequences (one sequence per library, in
the order given). The MIDs are included in the length specified with the
-read_dist option and can be altered by sequencing errors. See the MIDesigner or
BarCrawl programs to generate MID sequences.

=for Euclid:
   multiplex_ids.type: readable

=item -di <diversity>... | -diversity <diversity>...

This option specifies alpha diversity, specifically the richness, i.e. number of
reference sequences to take randomly and include in each library. Use 0 for the
maximum richness possible (based on the number of reference sequences available).
Provide one value to make all libraries have the same diversity, or one richness
value per library otherwise. Default: diversity.default

=for Euclid:
   diversity.type: 0+integer
   diversity.default: [ 0 ]

=item -sp <shared_perc> | -shared_perc <shared_perc>

This option controls an aspect of beta-diversity. When creating multiple
libraries, specify the percent of reference sequences they should have in common
(relative to the diversity of the least diverse library). Default: shared_perc.default %

=for Euclid:
   shared_perc.type: number, shared_perc >= 0 && shared_perc <= 100
   shared_perc.type.error: <shared_perc> must be a number between 0 and 100 (not shared_perc)
   shared_perc.default: 0

=item -pp <permuted_perc> | -permuted_perc <permuted_perc>

This option controls another aspect of beta-diversity. For multiple libraries,
choose the percent of the most-abundant reference sequences to permute (randomly
shuffle) the rank-abundance of. Default: permuted_perc.default %

=for Euclid:
   permuted_perc.type: number, permuted_perc >= 0 && permuted_perc <= 100
   permuted_perc.type.error: <permuted_perc> must be a number between 0 and 100 (not permuted_perc)
   permuted_perc.default: 100

=back

Miscellaneous

=over

=item -rs <random_seed> | -random_seed <random_seed>

Seed number to use for the pseudo-random number generator.

=for Euclid:
   random_seed.type: +integer

=item -dt <desc_track> | -desc_track <desc_track>

Track read information (reference sequence, position, errors, ...) by writing
it in the read description. Default: desc_track.default

=for Euclid:
   desc_track.type: integer, desc_track == 0 || desc_track == 1
   desc_track.type.error: <desc_track> must be 0 or 1 (not desc_track)
   desc_track.default: 1

=item -ql <qual_levels>... | -qual_levels <qual_levels>...

Generate basic quality scores for the simulated reads. Good residues are given a
specified good score (e.g. 30) and residues that are the result of an insertion
or substitution are given a specified bad score (e.g. 10). Specify first the
good score and then the bad score on the command-line, e.g.: 30 10. Default:
qual_levels.default

=for Euclid:
   qual_levels.type: 0+integer
   qual_levels.default: [ ]

=item -fq <fastq_output> | -fastq_output <fastq_output>

Whether to write the generated reads in FASTQ format (with Sanger-encoded
quality scores) instead of FASTA and QUAL or not (1: yes, 0: no).
<qual_levels> need to be specified for this option to be effective. Default: fastq_output.default

=for Euclid:
   fastq_output.type: integer, fastq_output == 0 || fastq_output == 1
   fastq_output.type.error: <fastq_output> must be 0 or 1 (not fastq_output)
   fastq_output.default: 0

=item -bn <base_name> | -base_name <base_name>

Prefix of the output files. Default: base_name.default

=for Euclid:
   base_name.type: string
   base_name.default: 'grinder'

=item -od <output_dir> | -output_dir <output_dir>

Directory where the results should be written. This folder will be created if
needed. Default: output_dir.default

=for Euclid:
   output_dir.type: writable
   output_dir.default: '.'

=item -pf <profile_file> | -profile_file <profile_file>

A file that contains Grinder arguments. This is useful if you use many options
or often use the same options. Lines with comments (#) are ignored. Consider the
profile file, 'simple_profile.txt':

  # A simple Grinder profile
  -read_dist 105 normal 12
  -total_reads 1000

Running: grinder -reference_file viral_genomes.fa -profile_file simple_profile.txt

Translates into: grinder -reference_file viral_genomes.fa -read_dist 105 normal 12 -total_reads 1000

Note that the arguments specified in the profile should not be specified again on the command line.

=back

=head1 CLI OUTPUT

For each shotgun or amplicon read library requested, the following files are
generated:

=over

=item *

A rank-abundance file, tab-delimited, that shows the relative abundance of the
different reference sequences

=item *

A file containing the read sequences in FASTA format. The read headers
contain information necessary to track from which reference sequence each read
was taken and what errors it contains. This file is not generated if <fastq_output>
option was provided.

=item *

If the <qual_levels> option was specified, a file containing the quality scores
of the reads (in QUAL format).

=item *

If the <fastq_output> option was provided, a file containing the read sequences
in FASTQ format.

=back

=head1 API EXAMPLES

The Grinder API allows to conveniently use Grinder within Perl scripts. Here is
a synopsis:

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

=head1 API METHODS

The rest of the documentation details the available Grinder API methods.

=head2 new

Title   : new

Function: Create a new Grinder factory initialized with the passed arguments.
          Available parameters described in the OPTIONS section.

Usage   : my $factory = Grinder->new( -reference_file => 'genomes.fna' );

Returns : a new Grinder object

=head2 next_lib

Title   : next_lib

Function: Go to the next shotgun library to process.

Usage   : my $struct = $factory->next_lib;

Returns : Community structure to be used for this library, where $struct->{ids}
          is an array reference containing the IDs of the genome making up the
          community (sorted by decreasing relative abundance) and $struct->{abs}
          is an array reference of the genome abundances (in the same order as
          the IDs).

=head2 next_read

Title   : next_read

Function: Create an amplicon or shotgun read for the current library.

Usage   : my $read  = $factory->next_read; # for single read
          my $mate1 = $factory->next_read; # for mate pairs
          my $mate2 = $factory->next_read;

Returns : A sequence represented as a Bio::Seq::SimulatedRead object

=head2 get_random_seed

Title   : get_random_seed

Function: Return the number used to seed the pseudo-random number generator

Usage   : my $seed = $factory->get_random_seed;

Returns : seed number


=head1 COPYRIGHT

Copyright 2009-2013 Florent ANGLY <florent.angly@gmail.com>

Grinder is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPL) as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Grinder is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Grinder.  If not, see <http://www.gnu.org/licenses/>.

=head1 BUGS

All complex software has bugs lurking in it, and this program is no exception.
If you find a bug, please report it on the SourceForge Tracker for Grinder:
L<http://sourceforge.net/tracker/?group_id=244196&atid=1124737>

Bug reports, suggestions and patches are welcome. Grinder's code is developed
on Sourceforge (L<http://sourceforge.net/scm/?type=git&group_id=244196>) and is
under Git revision control. To get started with a patch, do:

   git clone git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder

=cut

#---------- GRINDER FUNCTIONAL API --------------------------------------------#

sub Grinder {
  # This is the main function and is called by the script 'grinder'
  my (@args) = @_;

  # Create Grinder object
  my $factory = Grinder->new(@args);

  # Print diversity and percent shared and permuted
  diversity_report( $factory->{num_libraries}, $factory->{shared_perc},
    $factory->{permuted_perc}, $factory->{overall_diversity} );

  # Create the output directory if needed
  if ( not -d $factory->{output_dir} ) {
    mkdir $factory->{output_dir} or die "Error: Could not create output folder ".
      $factory->{output_dir}."\n$!\n";
  }

  # Generate sequences
  while ( my $c_struct = $factory->next_lib ) {
    my $cur_lib = $factory->{cur_lib};

    # Output filenames
    my $lib_str = '';
    if ($factory->{num_libraries} > 1) {
      $lib_str = '-'.sprintf('%0'.length($factory->{num_libraries}).'d', $cur_lib);
    }
    my $out_reads_basename = File::Spec->catfile($factory->{output_dir}, 
      $factory->{base_name}.$lib_str.'-reads.');
    my $out_fasta_file;
    my $out_qual_file;
    my $out_fastq_file;
    if ( $factory->{fastq_output} ) {
      $out_fastq_file = $out_reads_basename . 'fastq';
    } else {
      $out_fasta_file = $out_reads_basename . 'fa';
      if (scalar @{$factory->{qual_levels}} > 0) {
        $out_qual_file = $out_reads_basename . 'qual';
      }
    }
    my $out_ranks_file = File::Spec->catfile($factory->{output_dir},
      $factory->{base_name}.$lib_str."-ranks.txt");

    # Write community structure file
    $factory->write_community_structure($c_struct, $out_ranks_file);

    # Prepare output FASTA file
    my $out_fastq;
    if ( defined $out_fastq_file ) {
       $out_fastq = Bio::SeqIO->new( -format  => 'fastq',
                                     -variant => 'sanger',
                                     -flush   => 0,
                                     -file    => ">$out_fastq_file" );
    }

    my $out_fasta;
    if ( defined $out_fasta_file ) {
       $out_fasta = Bio::SeqIO->new( -format => 'fasta',
                                     -flush  => 0,
                                     -file   => ">$out_fasta_file" );
    }

    my $out_qual;
    if ( defined $out_qual_file ) {
       $out_qual = Bio::SeqIO->new( -format => 'qual',
                                    -flush  => 0,
                                    -file   => ">$out_qual_file" );
    }

    # Library report
    my $diversity = $factory->{diversity}[$cur_lib-1];
    library_report( $cur_lib, $factory->{alphabet}, $factory->{forward_reverse},
      $out_ranks_file, $out_fastq_file, $out_fasta_file, $out_qual_file,
      $factory->{cur_coverage_fold}, $factory->{cur_total_reads}, $diversity);

    # Generate shotgun or amplicon reads and write them to a file
    while ( my $read = $factory->next_read ) {
      $out_fastq->write_seq($read) if defined $out_fastq;
      $out_fasta->write_seq($read) if defined $out_fasta;
      $out_qual->write_seq($read)  if defined $out_qual
    }
    $out_fastq->close if defined $out_fastq;
    $out_fasta->close if defined $out_fasta;
    $out_qual->close  if defined $out_qual;
  }

  return 1;
}


sub diversity_report {
  my ($num_libraries, $perc_shared, $perc_permuted, $overall_diversity) = @_;
  my $format = '%.1f';
  print "Overall diversity = $overall_diversity genomes\n";
  if ($num_libraries > 1) {
    my $nof_shared  = $perc_shared / 100 * $overall_diversity;
    $perc_shared = sprintf($format, $perc_shared);
    print "Percent shared   = $perc_shared % ($nof_shared genomes)\n";
    my $nof_permuted  = $perc_permuted / 100 * $overall_diversity;
    $perc_permuted = sprintf($format, $perc_permuted);
    print "Percent permuted = $perc_permuted % ($nof_permuted top genomes)\n";
  }
  return 1;
}


sub write_community_structure {
  my ($self, $c_struct, $filename) = @_;
  open(OUT, ">$filename") || die("Error: Could not write in file $filename: $!\n");
  print OUT "# rank\tseq_id\trel_abund_perc\n";
  my $diversity = scalar @{$c_struct->{ids}};
  my %species_abs;
  for my $rank ( 1 .. $diversity ) {
    my $oid        = $c_struct->{'ids'}->[$rank-1];
    my $species_id = $self->database_get_parent_id($oid);
    my $seq_ab     = $c_struct->{'abs'}->[$rank-1];
    $species_abs{$species_id} += $seq_ab;
  }
  my $rank = 0;
  for my $species_id ( sort { $species_abs{$b} <=> $species_abs{$a} } keys %species_abs ) {
    $rank++;
    my $species_ab = $species_abs{$species_id};
    $species_ab *= 100; # in percentage
    print OUT "$rank\t$species_id\t$species_ab\n";
  }
  close OUT;
  return 1;
}


sub library_report {
  my ($cur_lib, $alphabet, $forward_reverse, $ranks_file, $fastq_file, $fasta_file,
    $qual_file, $coverage, $nof_seqs, $diversity) = @_;
  my $format = '%.3f';
  $coverage = sprintf($format, $coverage);
  my $lib_alphabet = uc $alphabet;
  $lib_alphabet =~ s/protein/Proteic/i;
  my $lib_type = defined $forward_reverse ? 'amplicon' : 'shotgun';
  print "$lib_alphabet $lib_type library $cur_lib:\n";
  print "   Community structure  = $ranks_file\n";
  print "   FASTQ file           = $fastq_file\n" if defined $fastq_file;
  print "   FASTA file           = $fasta_file\n" if defined $fasta_file;
  print "   QUAL file            = $qual_file\n"  if defined $qual_file;
  print "   Library coverage     = $coverage x\n";
  print "   Number of reads      = $nof_seqs\n";
  print "   Diversity (richness) = $diversity\n";
  return 1;
}


#---------- GRINDER OO API ----------------------------------------------------#


sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, ref($class) || $class;
  $self->argparse(\@args);
  $self->initialize();
  return $self;
}


sub next_lib {
  my ($self) = @_;
  $self->{cur_lib}++;
  $self->{cur_read} = 0;
  $self->{cur_total_reads} = 0;
  $self->{cur_coverage_fold} = 0;
  $self->{next_mate} = undef;
  $self->{positions} = undef;

  # Prepare sampling from this community
  my $c_struct = $self->{c_structs}[$self->{cur_lib}-1];
  if ( defined $c_struct ) {

    # Create probabilities of picking genomes from community structure
    $self->{positions} = $self->proba_create($c_struct, $self->{length_bias},
      $self->{copy_bias});

    # Calculate needed number of sequences based on desired coverage
    ($self->{cur_total_reads}, $self->{cur_coverage_fold}) = $self->lib_coverage($c_struct);

    # If chimeras are needed, update the kmer collection with sequence abundance
    my $kmer_col = $self->{chimera_kmer_col};
    if ($kmer_col) {
      my $weights;
      for (my $i = 0; $i < scalar @{$c_struct->{'ids'}}; $i++) {
        my $id     = $c_struct->{'ids'}->[$i];
        my $weight = $c_struct->{'abs'}->[$i];
        $weights->{$id} = $weight;
      }
      $kmer_col->weights($weights);
      my ($kmers, $freqs) = $kmer_col->counts(undef, undef, 1);
      $self->{chimera_kmer_arr} = $kmers;
      $self->{chimera_kmer_cdf} = $self->proba_cumul($freqs);
    }

  }

  return $c_struct;
}


sub next_read {
  my ($self) = @_;
  $self->next_lib if not $self->{cur_lib};
  $self->{cur_read}++;
  my $read;
  if ( $self->{cur_read} <= $self->{cur_total_reads} ) {
    # Generate the next read
    if ($self->{mate_length}) {
      # Generate a mate pair read
      if ( not $self->{next_mate} ) {
        # Generate a new pair of reads
        ($read, my $read2) = $self->next_mate_pair( );
        # Save second read of the pair for later
        $self->{next_mate} = $read2;
      } else {
        # Use saved read
        $read = $self->{next_mate};
        $self->{next_mate} = undef;
      }
    } else {
      # Generate a single shotgun or amplicon read
      $read = $self->next_single_read( );
    }
  }
  return $read;
}


sub get_random_seed {
  my ($self) = @_;
  return $self->{random_seed};
}


#---------- GRINDER INTERNALS -------------------------------------------------#


sub argparse {
  # Process arguments
  my ($self, $args) = @_;
  my @old_args = @$args;
  # Read profile file
  $args = process_profile_file($args);
  # Parse and validate arguments with Getopt::Euclid
  Getopt::Euclid->process_args($args);
  # Check that Euclid worked, i.e. that there is at least one parameter in %ARGV
  if ( scalar keys %ARGV == 0 ) {
    die "Error: the command line arguments could not be parsed because of an ".
      "internal problem\n";
  }
  # Get parsed arguments from %ARGV and put them in $self
  for my $arg (keys %ARGV) {
    # Skip short argument names (they are also represented with long names)
    next if length($arg) <= 2;
    # Process long argument names. Copy their value into $self
    my $ref = ref($ARGV{$arg});
    if (not $ref) {
      $self->{$arg} = $ARGV{$arg};
    } elsif ($ref eq 'ARRAY') {
      @{$self->{$arg}} = @{$ARGV{$arg}};
    } else {
      die "Error: unsupported operation on argument '$arg' which is a reference".
         "of type $ref\n";
    }
  }
  return 1;
}


sub process_profile_file {
  # Find profile file in arguments and read the profiles. The profile file
  # only contains Grinder arguments, and lines starting with a '#' are comments.
  my ($args) = @_;
  my $file;
  for (my $i = 0; $i < scalar @$args; $i++) {
    my $arg = $$args[$i];
    if ($arg =~ m/^-profile_file/ || $arg =~ m/-pf/) {
      $file = $$args[$i+1];
      if ( (not defined $file) || ($file =~ m/^-/) ) {
        die "Error: no value was given to --profile_file\n";
      }
    }
  }
  if (defined $file) {
    open my $in, '<', $file or die "Error: Could not read file '$file'\n$!\n";
    my $profile = '';
    while (my $line = <$in>) {
      chomp $line;
      next if $line =~ m/^\s*$/;
      next if $line =~ m/^\s*#/;
      $profile .= "$line ";
    }
    close $in;
    push @$args, split /\s+/, $profile;
  }
  return $args;
}


sub initialize {
  my ($self) = @_;
  # Returns:

  # Parameter processing: read length distribution
  if ( (not ref $self->{read_dist}) or (ref $self->{read_dist} eq 'SCALAR') ){
    $self->{read_dist}   = [$self->{read_dist}];
  }
  $self->{read_length} = $self->{read_dist}[0] || 100;
  $self->{read_model}  = $self->{read_dist}[1] || 'uniform';
  $self->{read_delta}  = $self->{read_dist}[2] || 0;
  delete $self->{read_dist};

  # Parameter processing: mate insert length distribution
  if ( (not ref $self->{insert_dist}) or (ref $self->{insert_dist} eq 'SCALAR') ){
    $self->{insert_dist} = [$self->{insert_dist}];
  }
  $self->{mate_length} = $self->{insert_dist}[0] || 0;
  $self->{mate_model}  = $self->{insert_dist}[1] || 'uniform';
  $self->{mate_delta}  = $self->{insert_dist}[2] || 0;
  delete $self->{insert_dist};

  # Parameter processing: genome abundance distribution
  if ( (not ref $self->{abundance_model}) or (ref $self->{abundance_model} eq 'SCALAR') ){
    $self->{abundance_model} = [$self->{abundance_model}];
  }
  $self->{distrib} = $self->{abundance_model}[0] || 'uniform';
  $self->{param}   = $self->{abundance_model}[1];
  delete $self->{abundance_model};

  # Parameter processing: point sequencing error distribution
  if ( (not ref $self->{mutation_dist}) or (ref $self->{mutation_dist}  eq 'SCALAR') ) {
    $self->{mutation_dist} = [$self->{mutation_dist}];
  }
  $self->{mutation_model} = $self->{mutation_dist}[0] || 'uniform';
  $self->{mutation_para1} = $self->{mutation_dist}[1] || 0;
  $self->{mutation_para2} = $self->{mutation_dist}[2] || 0;
  delete $self->{mutation_dist};

  # Parameter processing: mutation ratio
  $self->{mutation_ratio}[0] = $self->{mutation_ratio}[0] || 0;
  $self->{mutation_ratio}[1] = $self->{mutation_ratio}[1] || 0;
  my $sum = $self->{mutation_ratio}[0] + $self->{mutation_ratio}[1];
  if ($sum == 0) {
    $self->{mutation_ratio}[0] = 50;
    $self->{mutation_ratio}[1] = 50;
  } else {
    $self->{mutation_ratio}[0] = $self->{mutation_ratio}[0] *100 / $sum;
    $self->{mutation_ratio}[1] = $self->{mutation_ratio}[1] *100 / $sum;
  }

  # Parameter processing: homopolymer model
  $self->{homopolymer_dist} = lc $self->{homopolymer_dist} if defined $self->{homopolymer_dist};
  
  # Parameter processing: chimera_distribution
  if ( (not ref $self->{chimera_dist}) or (ref $self->{chimera_dist}  eq 'SCALAR') ) {
    $self->{chimera_dist} = [$self->{chimera_dist}];
  }
  if ($self->{chimera_dist}) {
    # Normalize to 1
    my $total = 0;
    for my $multimera_abundance (@{$self->{chimera_dist}}) {
      $total += $multimera_abundance;
    }
    $self->{chimera_dist} = undef if $total == 0;
    $self->{chimera_dist} = normalize($self->{chimera_dist}, $total);
    # Calculate cdf
    if ($self->{chimera_perc}) {
      $self->{chimera_dist_cdf} = $self->proba_cumul( $self->{chimera_dist} );
    }
  }

  # Parameter processing: fastq_output requires qual_levels
  if ( ($self->{fastq_output}) &&
       (not scalar @{$self->{qual_levels}} > 0) ) {
    die "Error: <qual_levels> needs to be specified to output FASTQ reads\n";
  }

  # Random number generator: seed or be auto-seeded
  if (defined $self->{random_seed}) {
    srand( $self->{random_seed} );
  } else {
    $self->{random_seed} = srand( );
  }

  # Sequence length check
  my $max_read_length = $self->{read_length} + $self->{read_delta}; # approximation
  if ($self->{mate_length}) {
    my $min_mate_length = $self->{mate_length} - $self->{mate_delta};
    if ($max_read_length > $min_mate_length) {
      die("Error: The mate insert length cannot be smaller than read length. ".
        "Try increasing the mate insert length or decreasing the read length\n");
    }
  }

  # Pre-compile regular expression to check if reads are valid
  if ( (defined $self->{exclude_chars}) && (not $self->{exclude_chars} eq '') ) {
    $self->{exclude_re} = qr/[${\$self->{exclude_chars}}]/i; # Match any of the chars
  }
  
  # Read MIDs
  $self->{multiplex_ids} = $self->read_multiplex_id_file($self->{multiplex_ids}, 
    $self->{num_libraries}) if defined $self->{multiplex_ids};

  # Import reference sequences
  my $min_seq_len;
  if ($self->{chimera_dist_cdf}) {
    # Each chimera needs >= 1 bp. Use # sequences required by largest chimera.
    $min_seq_len = scalar @{$self->{chimera_dist}} + 1;
  } else {
    $min_seq_len = 1;
  }
  $self->{database} = $self->database_create( $self->{reference_file},
    $self->{unidirectional}, $self->{forward_reverse}, $self->{abundance_file},
    $self->{delete_chars}, $min_seq_len );

  $self->initialize_alphabet;
  if ( ($self->{alphabet} eq 'protein')     &&
       ($self->{mate_length} != 0)          &&
       (not $self->{mate_orientation} eq 'FF') ) {
    die "Error: Can only use <mate_orientation> FF with proteic reference sequences\n";
  }

  # Genome relative abundance in the different independent libraries to create
  $self->{c_structs} = $self->community_structures( $self->{database}->{ids},
    $self->{abundance_file}, $self->{distrib}, $self->{param},
    $self->{num_libraries}, $self->{shared_perc}, $self->{permuted_perc},
    $self->{diversity}, $self->{forward_reverse} );

  # Count kmers in the database if we need to form kmer-based chimeras
  if ($self->{chimera_perc} && $self->{chimera_kmer}) {

    # Get all wanted sequences (not all the sequences in the database)
    my %ids_hash;
    my @ids;
    my @seqs;
    for my $c_struct ( @{ $self->{c_structs} } ) {
      for my $id (@{$c_struct->{ids}}) {
        if (not exists $ids_hash{$id}) {
          $ids_hash{$id} = undef;
          push @ids, $id;
          push @seqs, $self->database_get_seq($id);
        }
      }
    }
    %ids_hash = ();

    # Now create a collection of kmers
    $self->{chimera_kmer_col} = Grinder::KmerCollection->new(
      -k    => $self->{chimera_kmer},
      -seqs => \@seqs,
      -ids  => \@ids,
    )->filter_shared(2);
  }

  # Markers to keep track of computation progress
  $self->{cur_lib}  = 0;
  $self->{cur_read} = 0;

  return $self;
}


sub initialize_alphabet {
  # Store the characters of the alphabet to use and calculate their cdf so that
  # we can easily pick them at random later
  my ($self) = @_;
  my $alphabet = $self->{alphabet};
  # Characters available in alphabet
  my %alphabet_hash;
  if ($alphabet eq 'dna') {
    %alphabet_hash = (
      'A' => undef,
      'C' => undef,
      'G' => undef,
      'T' => undef,
    );
  } elsif ($alphabet eq 'rna') {
    %alphabet_hash = (
      'A' => undef,
      'C' => undef,
      'G' => undef,
      'U' => undef,
    );
  } elsif ($alphabet eq 'protein') {
    %alphabet_hash = (
      'A' => undef,
      'R' => undef,
      'N' => undef,
      'D' => undef,
      'C' => undef,
      'Q' => undef,
      'E' => undef,
      'G' => undef,
      'H' => undef,
      'I' => undef,
      'L' => undef,
      'K' => undef,
      'M' => undef,
      'F' => undef,
      'P' => undef,
      'S' => undef,
      'T' => undef,
      'W' => undef,
      'Y' => undef,
      'V' => undef,
      #'B' => undef, # D or N
      #'Z' => undef, # Q or E
      #'X' => undef, # any amino-acid
      # J, O and U are the only unused letters
    );
  } else {
    die "Error: unknown alphabet '$alphabet'\n";
  }
  my $num_chars = scalar keys %alphabet_hash;
  $self->{alphabet_hash} = \%alphabet_hash;
  # CDF for this alphabet
  $self->{alphabet_complete_cdf} = $self->proba_cumul([(1/$num_chars) x $num_chars]);
  $self->{alphabet_truncated_cdf} = $self->proba_cumul([(1/($num_chars-1)) x ($num_chars-1)]);
  return 1;
}


sub read_multiplex_id_file {
  my ($self, $file, $nof_indep) = @_;
  my @mids;
  # Read FASTA file containing the MIDs
  my $in = Bio::SeqIO->newFh(
    -file   => $file,
    -format => 'fasta',
  );
  while (my $mid = <$in>) {
    push @mids, $mid->seq;
  }
  undef $in;
  # Sanity check
  my $nof_mids = scalar @mids;
  if ($nof_mids < $nof_indep) {
    die "Error: $nof_indep communities were requested but the MID file ".
      "had only $nof_mids sequences.\n"; 
  } elsif ($nof_mids > $nof_indep) {
    warn "Warning: $nof_indep communities were requested but the MID file ".
      "contained $nof_mids sequences. Ignoring extraneous MIDs...\n";
  }
  return \@mids;
}


sub community_structures {
  # Create communities with a specified structure, alpha and beta-diversity
  my ($self, $seq_ids, $abundance_file, $distrib, $param, $nof_indep,
    $perc_shared, $perc_permuted, $diversities, $forward_reverse) = @_;

  # Calculate community structures
  my $c_structs;
  if ($abundance_file) {
    # Sanity check
    if ( (scalar @$diversities > 1) || $$diversities[0] ) {
      warn "Warning: Diversity cannot be specified when an abundance file is specified. Ignoring it...\n";
    }
    if ( ($perc_shared > 0) || ($perc_permuted < 100) ) {
      warn "Warning: Percent shared and percent permuted cannot be specified when an abundance file is specified. Ignoring them...\n";
    }
    # One or several communities with specified rank-abundances
    $c_structs = community_given_abundances($abundance_file, $seq_ids);

    # Calculate number of libraries
    my $got_indep = scalar @$c_structs;
    if ($nof_indep != 1) { # 1 is the default value
      if ($nof_indep > $got_indep) {
        die "Error: $nof_indep communities were requested but the abundance file".
          " specified the abundances for only $got_indep.\n";
      } elsif ($nof_indep < $got_indep) {
        warn "Warning: $nof_indep communities were requested by the abundance ".
          "file specified the abundances for $got_indep. Ignoring extraneous ".
          "communities specified in the file.\n";
      }
    }
    $nof_indep = $got_indep;
    $self->{num_libraries} = $nof_indep;
    # Calculate diversities based on given community abundances
    ($self->{diversity}, $self->{overall_diversity}, $self->{shared_perc},
      $self->{permuted_perc}) = community_calculate_diversities($c_structs);
  } else {
    # One or several communities with rank-abundance to be calculated
    # Sanity check
    if ($nof_indep == 1) { # 1 is the default value
      $nof_indep = scalar @$diversities;
    }
    if ($nof_indep != scalar @$diversities) {
      if (scalar @$diversities == 1) {
        # Use same diversity for all libraries
        my $diversity = $$diversities[0];
        for my $i (1 .. $nof_indep-1) {
          push @$diversities, $diversity;
        }
      } else {
        die "Error: The number of richness values provided (".(scalar @$diversities).
           ") did not match the requested number of libraries ($nof_indep).\n";
      }
    }
    $self->{num_libraries} = $nof_indep;
    # Select shared species
    my $c_ids;
    my $overall_diversity = 0;
    ($c_ids, $overall_diversity, $diversities, $perc_shared) = community_shared(
      $seq_ids, $nof_indep, $perc_shared, $diversities );

    # Shuffle the abundance-ranks of the most abundant genomes
    ($c_ids, $perc_permuted) = community_permuted($c_ids, $perc_permuted);
    # Update values in $self object
    $self->{overall_diversity} = $overall_diversity;
    $self->{diversity} = $diversities;
    $self->{shared_perc} = $perc_shared;
    $self->{permuted_perc} = $perc_permuted;
    # Put results in a community structure "object"
    for my $c (1 .. $nof_indep) {
      # Assign a random parameter if needed
      my $comm_param = defined $param ? $param : randig(1,0.05);
      # Calculate relative abundance of the community members
      my $diversity = $self->{diversity}[$c-1];
      my $c_abs = community_calculate_species_abundance($distrib, $comm_param,
         $diversity);
      my $c_ids = $$c_ids[$c-1];
      my $c_struct;
      $c_struct->{'ids'}   = $c_ids;
      $c_struct->{'abs'}   = $c_abs;
      $c_struct->{'param'} = $comm_param;
      $c_struct->{'model'} = $distrib;
      push @$c_structs, $c_struct;
    }
  }

  # Convert sequence IDs to object IDs
  for my $c_struct (@$c_structs) {
    ($c_struct->{'abs'}, $c_struct->{'ids'}) = community_calculate_amplicon_abundance(
      $c_struct->{'abs'}, $c_struct->{'ids'}, $seq_ids );
  }

  return $c_structs;
}


sub community_calculate_diversities {
  my ($c_structs) = @_;
  my ($diversities, $overall_diversity, $perc_shared, $perc_permuted) = (0, 0, 0, 0);

  # Calculate diversity (richness) based on given community abundances
  my $nof_libs = scalar @$c_structs;
  my %all_ids;
  my @richnesses;
  for my $c_struct (@$c_structs) {
    my $richness = 0;
    for my $i (0 .. scalar @{$$c_struct{ids}} - 1) {
      my $id = $$c_struct{ids}[$i];
      my $ab = $$c_struct{abs}[$i];
      next if not $ab;
      $richness++;
      
      if (defined $all_ids{$id}) {
        $all_ids{$id}++;
      } else {
        $all_ids{$id} = 1;
      }
    }
    push @richnesses, $richness;
  }
  $overall_diversity = scalar keys %all_ids;


  # Calculate percent shared
  my $nof_non_shared = 0;
  for my $id (keys %all_ids) {
    $nof_non_shared++ if $all_ids{$id} < $nof_libs;
  }
  $perc_shared = ($overall_diversity - $nof_non_shared) * 100 / $overall_diversity;

  # TODO: Could calculate percent permuted

  return \@richnesses, $overall_diversity, $perc_shared, $perc_permuted;
}


sub community_given_abundances {
  # Read a file of genome abundances. The file should be space or tab-delimited. 
  # The first column should be the IDs of genomes, and the subsequent columns is
  # for their relative abundance in different communities. An optional list of
  # valid IDs can be provided. Then the abundances are normalized so that their
  # sum is 1.
  my ($file, $seq_ids) = @_;

  # Read abundances
  my ($ids, $abs) = community_read_abundances($file);
  # Remove genomes with unknown IDs and calculate cumulative abundance
  my $totals;
  for my $comm_num (0 .. $#$ids) {
    my $i = 0;
    while ( $i < scalar @{$$ids[$comm_num]} ) {
      my $id = $$ids[$comm_num][$i];
      my $ab = $$abs[$comm_num][$i];
      if ( (scalar keys %$seq_ids == 0) || (exists $$seq_ids{$id}) ) {
        $$totals[$comm_num] += $ab;
        $i++;
      } else {
        die "Error: Requested reference sequence '$id' in file '$file' does not".
          " exist in the input database.\n";
        splice @{$$ids[$comm_num]}, $i, 1;
        splice @{$$abs[$comm_num]}, $i, 1;
      }
    }
  }
  # Process the communities
  my @c_structs;
  for my $comm_num (0 .. scalar @$ids - 1) {
    my $comm_ids   = $$ids[$comm_num];
    my $comm_abs   = $$abs[$comm_num];
    my $comm_total = $$totals[$comm_num];
    if ($comm_total == 0) {
      warn "Warning: The abundance of all the genomes for community ".($comm_num+1)." was zero. Skipping this community...\n";
      next;
    }
    # Normalize the abundances
    $comm_abs = normalize($comm_abs, $comm_total);
    # Sort relative abundances by decreasing 
    ($comm_abs, $comm_ids) = two_array_sort($comm_abs, $comm_ids);
    $comm_abs = [reverse(@$comm_abs)];
    $comm_ids = [reverse(@$comm_ids)];
    # Save community structure
    my $c_struct = { 'ids' => $comm_ids, 'abs' => $comm_abs };
    push @c_structs, $c_struct;
  }
  return \@c_structs;
}


sub community_read_abundances {
  my ($file) = @_;
  # Read abundances of genomes from a file
  my $ids; # genome IDs
  my $abs; # genome relative abundance
  open my $io, '<', $file or die "Error: Could not read file '$file'\n$!\n";
  while ( my $line = <$io> ) {
    # Ignore comment or empty lines
    if ( $line =~ m/^\s*$/ || $line =~ m/^#/ ) {
      next;
    }
    # Read abundance info from line
    my ($id, @rel_abs) = ($line =~ m/(\S+)/g);
    if (defined $id) {
      for my $comm_num (0 .. $#rel_abs) {
        my $rel_ab = $rel_abs[$comm_num];
        push @{$$ids[$comm_num]}, $id;
        push @{$$abs[$comm_num]}, $rel_ab;
      }
    } else {
      warn "Warning: Line $. of file '$file' has an unknown format. Skipping it...\n";
    }
  }
  close $io;
  return $ids, $abs;
}


sub community_permuted {
  # Change the abundance rank of species in all but the first community.
  # The number of species changed in abundance is determined by the percent
  # permuted, i.e. a given percentage of the most abundant species in this community.
  my ($c_ids, $perc_permuted) = @_;
  my $nof_indep = scalar @$c_ids;

  # Leave the first community alone, but permute the ones after
  for my $c ( 2 .. $nof_indep ) {
    my $ids = $$c_ids[$c-1];
    my $diversity = scalar @$ids;
    # Number of top genomes to permute
    # Percent permuted is relative to diversity in this community
    my $nof_permuted = $perc_permuted / 100 * $diversity;
    $nof_permuted = int($nof_permuted + 0.5); # round number

    # Method published in Angly et al 2006 PLOS Biology    
    # Take the $nof_permuted first ranks (most abundant genomes) and shuffle
    # (permute) their ranks amongst the $nof_permuted first ranks.
    # Caveat: cannot permute only 1 genome
    my $idxs;
    if ($nof_permuted > 0) {
      # Add shuffled top genomes
      my $permuted_idxs = randomize( [0 .. $nof_permuted-1] );
      push @$idxs, @$permuted_idxs;
    }
    if ($diversity - $nof_permuted > 0) {
      # Add other genomes in same order
      my $non_permuted_idxs = [$nof_permuted .. $diversity-1];
      push @$idxs, @$non_permuted_idxs;
    }
    @$ids = @$ids [ @$idxs ];

  }

  return $c_ids, $perc_permuted;
}


sub community_shared {
  # Randomly split a library of sequences into a given number of groups that
  # share a specified percent of their genomes.
  # The % shared is the number of species shared / the total diversity in all communities
  # Input: arrayref of sequence ids
  #        number of communities to produce
  #        percentage of genomes shared between the communities
  #        diversity (optional, will use all genomes if not specified) 
  # Return: arrayref of IDs that are shared
  #         arrayref of arrayref with the unique IDs for each community
  my ($seq_ids, $nof_indep, $perc_shared, $diversities) = @_;

  # If diversity is not specified (is '0'), use the maximum value possible
  my $nof_refs = scalar keys %$seq_ids;
  my $min_diversity = 1E99;
  for my $i (0 .. scalar @$diversities - 1) {
    if ($$diversities[$i] == 0) {
      $$diversities[$i] = $nof_refs / ( $perc_shared/100 + $nof_indep*(1-$perc_shared/100) );
      $$diversities[$i] = int( $$diversities[$i] );
      if ( ($i > 0) && ($$diversities[$i-1] != $$diversities[$i]) ) {
        die "Error: Define either all the diversities or none.\n";
      }
    }
    if ($$diversities[$i] < $min_diversity) {
      $min_diversity = $$diversities[$i];
    }
  }

  if ($min_diversity == 0) {
    die "Error: Cannot make $nof_indep libraries sharing $perc_shared % species".
      " from $nof_refs references\n";
  }

  # Calculate the number of sequences to share, noting that the percent shared
  # is relative to the diversity of the least abundant library
  my $nof_shared = int($min_diversity * $perc_shared / 100);
  $perc_shared = $nof_shared * 100 / $min_diversity;

  # Unique sequences
  my @nof_uniques;
  my $sum_not_uniques = 0;
  for my $diversity (@$diversities) {
    my $nof_unique = $diversity - $nof_shared;
    $sum_not_uniques += $nof_unique;
    push @nof_uniques, $nof_unique;
  }

  # Overall diversity
  my $overall_diversity = $nof_shared + $sum_not_uniques;
  if ($nof_refs <  $overall_diversity) {
    die "Error: The number of reference sequences available ($nof_refs) is not".
      " large enough to support the requested diversity ($overall_diversity ".
      "genomes overall with $perc_shared % genomes shared between $nof_indep ".
      "libraries)\n";
  }

  # Add shared sequences
  my @ids = keys %$seq_ids;
  my @shared_ids;
  for (0 .. $nof_shared - 1) {
    # Pick a random sequence
    my $rand_offset = int(rand($nof_refs));
    my $rand_id = splice @ids, $rand_offset, 1;
    $nof_refs = scalar(@ids);
    # Add this sequence in all independent libraries
    push @shared_ids, $rand_id;
  }

  # Add sequences not shared
  my @unique_ids;
  for my $lib_num (0 .. $nof_indep-1) {
    my $nof_unique = $nof_uniques[$lib_num];
    for (0 .. $nof_unique - 1) {
      # Pick a random sequence
      my $rand_offset = int(rand($nof_refs));
      my $rand_id = splice @ids, $rand_offset, 1;
      $nof_refs = scalar(@ids);
      # Add this sequence in this independent library only
      push @{$unique_ids[$lib_num]}, $rand_id;
    }
  }

  # Randomly pick the rank of the shared IDs
  my $shared_ranks = randomize( [1 .. $min_diversity] );
  @$shared_ranks = splice @$shared_ranks, 0, $nof_shared;

  # Construct community ranks
  my @c_ranks;
  for my $lib_num (0 .. $nof_indep-1) {
    my $diversity = $$diversities[$lib_num];
    my @ranks = (undef) x $diversity;
    # Add shared IDs
    for my $i (0 .. $nof_shared-1) {
      my $id   = $shared_ids[$i];
      my $rank = $$shared_ranks[$i];
      $ranks[$rank-1] = $id;
    }
    # Add unique IDs
    my $ids = $unique_ids[$lib_num];
    for my $rank (1 .. $diversity) {
      next if defined $ranks[$rank-1];
      $ranks[$rank-1] = pop @$ids;   
    }
    push @c_ranks, \@ranks;
  }

  return \@c_ranks, $overall_diversity, $diversities, $perc_shared;
}


sub community_calculate_species_abundance {
  # Calculate relative abundance based on a distribution and its parameters.
  # Input is a model, its 2 parameters, and the number of values to generate
  # Output is a reference to a list of relative abundance. The abundance adds up
  # to 1
  my ($distrib, $param, $diversity) = @_;
  # First calculate rank-abundance values
  my $rel_ab;
  my $total = 0;
  if ($distrib eq 'uniform') {
    # no parameter
    my $val = 1 / $diversity;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = $val;
    }
    $total = 1;
  } elsif ($distrib eq 'linear') {
    # no parameter
    my $slope = 1 / $diversity;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = 1 - $slope * $index;
      $total += $$rel_ab[$index];
    }
  } elsif ($distrib eq 'powerlaw') {
    # 1 parameter
    die "Error: The powerlaw model requires an input parameter (-p option)\n"
      if not defined $param;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = ($index+1)**-$param;
      $total += $$rel_ab[$index];
    }
  } elsif ($distrib eq 'logarithmic') {
    # 1 parameter
    die "Error: The logarithmic model requires an input parameter (-p option)\n"
      if not defined $param;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = log($index+2)**-$param;
      $total += $$rel_ab[$index];
    }
  } elsif ($distrib eq 'exponential') {
    # 1 parameter
    die "Error: The exponential model requires an input parameter (-p option)\n"
      if not defined $param;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = exp(-($index+1)*$param);
      $total += $$rel_ab[$index];
    }
  } else {
    die "Error: $distrib is not a valid rank-abundance distribution\n";
  }
  # Normalize to 1 if needed
  if ($total != 1) {
    $rel_ab = normalize($rel_ab, $total);
  }
  return $rel_ab;
}


sub community_calculate_amplicon_abundance {
  my ($r_spp_abs, $r_spp_ids, $seq_ids) = @_;
  # Convert abundance of species into abundance of their amplicons because there
  # can be multiple amplicon per species and the amplicons have a different ID
  # from the species. The r_spp_ids and r_spp_abs arrays are the ID and abundance
  # of the species, sorted by decreasing abundance.

  # Give amplicons from the same species the same sampling probability
  for (my $i  = 0; $i < scalar @$r_spp_ids; $i++) {
    my $species_ab    = $$r_spp_abs[$i];
    my $species_id    = $$r_spp_ids[$i];
    my @amplicon_ids  = keys %{$seq_ids->{$species_id}};
    my $nof_amplicons = scalar @amplicon_ids;
    my @amplicon_abs  = ($species_ab / $nof_amplicons) x $nof_amplicons;
    splice @$r_spp_abs, $i, 1, @amplicon_abs;
    splice @$r_spp_ids, $i, 1, @amplicon_ids;
    $i += $nof_amplicons - 1;
  }

  return $r_spp_abs, $r_spp_ids;
}


sub next_single_read {
  # Generate a single shotgun or amplicon read
  my ($self) = @_;
  my $oids           = $self->{c_structs}->[$self->{cur_lib}-1]->{ids};
  my $mid            = $self->{multiplex_ids}->[$self->{cur_lib}-1] || '';
  my $lib_num        = $self->{num_libraries} > 1 ? $self->{cur_lib} : undef;
  my $max_nof_tries  = $self->{forward_reverse} ? 1 : 10;

  # Choose a random genome or amplicon
  my $genome = $self->rand_seq($self->{positions}, $oids);
  my $nof_tries = 0;
  my $shotgun_seq;
  do {
    # Error if we have exceeded the maximum number attempts
    $nof_tries++;
    if ($nof_tries > $max_nof_tries) {
      my $message = "Error: Could not take a random shotgun read without ".
        "forbidden characters from reference sequence ".$genome->seq->id;
      $message .= " ($max_nof_tries attempts made)" if ($max_nof_tries > 1);
      $message .= ".\n";
      die $message;
    }

    # Chimerize the template sequence if needed
    $genome = $self->rand_seq_chimera($genome, $self->{chimera_perc},
      $self->{positions}, $oids) if $self->{chimera_perc};

    # Take a random orientation if needed
    my $orientation = ($self->{unidirectional} != 0) ? 1 : rand_seq_orientation();
    # Choose a read size according to the specified distribution
    my $length = rand_seq_length($self->{read_length}, $self->{read_model},
      $self->{read_delta});

    # Shorten read length if too long
    my $max_length = $genome->length + length($mid);
    if ( $length > $max_length) {
      $length = $max_length;
    }

    # Read position on genome or amplicon
    my ($start, $end) = rand_seq_pos($genome, $length, $self->{forward_reverse},
      $mid);

    # New sequence object
    $shotgun_seq = new_subseq($self->{cur_read}, $genome, $self->{unidirectional},
      $orientation, $start, $end, $mid, undef, $lib_num, $self->{desc_track},
      $self->{qual_levels});

    # Simulate sequence aberrations and sequencing error if needed
    $shotgun_seq = $self->rand_seq_errors($shotgun_seq)
      if ($self->{homopolymer_dist} || $self->{mutation_para1});
  } while (
    $self->{exclude_re} && not $self->is_valid($shotgun_seq)
  );
  return $shotgun_seq;
}


sub next_mate_pair {
  # Generate a shotgun mate pair
  my ($self) = @_;
  my $oids           = $self->{c_structs}->[$self->{cur_lib}-1]->{ids};
  my $mid            = $self->{multiplex_ids}->[$self->{cur_lib}-1] || '';
  my $lib_num        = $self->{num_libraries} > 1 ? $self->{cur_lib} : undef;
  my $pair_num       = int( $self->{cur_read} / 2 + 0.5 );
  my $max_nof_tries  = $self->{forward_reverse} ? 1 : 10;

  # Deal with mate orientation
  my @mate_orientations = split('', $self->{mate_orientation} );
  my $mate_1_orientation = $mate_orientations[0] eq 'F' ? 1 : -1;
  my $mate_2_orientation = $mate_orientations[1] eq 'F' ? 1 : -1;

  # Choose a random genome
  my $genome = $self->rand_seq($self->{positions}, $oids);

  my $nof_tries = 0;
  my ($shotgun_seq_1, $shotgun_seq_2);
  while (1) {
    # Error if we have exceeded the maximum number of attempts
    $nof_tries++;
    if ($nof_tries > $max_nof_tries) {
      my $message = "Error: Could not take a pair of random shotgun read ".
        "without forbidden characters from reference sequence ".$genome->seq->id;
      $message .= " ($max_nof_tries attempts made)" if ($max_nof_tries > 1);
      $message .= ".\n";
      die $message;
    }
    # Chimerize the template sequence if needed
    $genome = $self->rand_seq_chimera($genome, $self->{chimera_perc},
      $self->{positions}, $oids) if $self->{chimera_perc};
    # Take from a random strand if needed
    my $orientation = ($self->{unidirectional} != 0) ? 1 : rand_seq_orientation();
    # Choose a mate pair length according to the specified distribution
    my $mate_length = rand_seq_length($self->{mate_length}, $self->{mate_model},
      $self->{mate_delta});
    # Shorten mate length if too long
    my $max_length = $genome->length + length($mid);
    if ( $mate_length > $max_length) {
      $mate_length = $max_length;
    }
    # Mate position on genome or amplicon
    my ($mate_start, $mate_end) = rand_seq_pos($genome, $mate_length,
      $self->{forward_reverse}, $mid);
    # Determine mate-pair position
    my $read_length = rand_seq_length($self->{read_length}, $self->{read_model}, $self->{read_delta});
    my $seq_1_start = $mate_start;
    my $seq_1_end   = $mate_start + $read_length - 1;
    $read_length    = rand_seq_length($self->{read_length}, $self->{read_model}, $self->{read_delta});
    my $seq_2_start = $mate_end - $read_length + 1;
    my $seq_2_end   = $mate_end;
    if ($orientation == -1) {
       $mate_1_orientation = $orientation * $mate_1_orientation;
       $mate_2_orientation = $orientation * $mate_2_orientation;
       ($seq_1_start, $seq_2_start) = ($seq_2_start, $seq_1_start);
       ($seq_1_end  , $seq_2_end  ) = ($seq_2_end  , $seq_1_end  );
    }
    # Generate first mate read
    $shotgun_seq_1 = new_subseq($pair_num, $genome, $self->{unidirectional},
      $mate_1_orientation, $seq_1_start, $seq_1_end, $mid, '1', $lib_num, $self->{desc_track},
      $self->{qual_levels});
    $shotgun_seq_1 = $self->rand_seq_errors($shotgun_seq_1)
      if ($self->{homopolymer_dist} || $self->{mutation_para1});
    if ($self->{exclude_re} && not $self->is_valid($shotgun_seq_1)) {
      next;
    }
    # Generate second mate read
    $shotgun_seq_2 = new_subseq($pair_num, $genome, $self->{unidirectional},
      $mate_2_orientation, $seq_2_start, $seq_2_end, $mid, '2', $lib_num, $self->{desc_track},
      $self->{qual_levels});
    $shotgun_seq_2 = $self->rand_seq_errors($shotgun_seq_2)
      if ($self->{homopolymer_dist} || $self->{mutation_para1});
    if ($self->{exclude_re} && not $self->is_valid($shotgun_seq_2)) {
      next;
    }
    # Both shotgun reads were valid
    last;
  }
  return $shotgun_seq_1, $shotgun_seq_2;
}


sub is_valid {
  # Return 1 if the sequence object is valid (is not empty and does not have any
  # of the specified forbidden characters), 0 otherwise. Specify the forbidden
  # characters as a single string, e.g. 'N-' to prevent any reads to have 'N' or
  # '-'. The search is case-insensitive.
  my ($self, $seq) = @_;
  if ($seq->seq =~ $self->{exclude_re}) {
    return 0;
  }
  return 1;
}


sub proba_create {
  my ($self, $c_struct, $size_dep, $copy_bias) = @_;
  # 1/ Calculate size-dependent, copy number-dependent probabilities
  my $probas = $self->proba_bias_dependency($c_struct, $size_dep, $copy_bias);
  # 2/ Generate proba starting position
  my $positions = $self->proba_cumul($probas);
  return $positions;
}


sub proba_bias_dependency {
  # Affect probability of picking a species by considering genome length or gene
  # copy number bias
  my ($self, $c_struct, $size_dep, $copy_bias) = @_;

  # Calculate probability
  my $probas;
  my $totproba = 0;
  my $diversity = scalar @{$c_struct->{'ids'}};
  for my $i (0 .. scalar $diversity - 1) {
    my $proba = $c_struct->{'abs'}[$i];

    if ( defined $self->{forward_reverse} ) {
      # Gene copy number bias
      if ($copy_bias) {
        my $refseq_id = $self->database_get_parent_id($c_struct->{'ids'}[$i]);
        my $nof_amplicons = scalar @{ $self->database_get_children_seq($refseq_id) };
        $proba *= $nof_amplicons;
      }
    } else {
      # Genome length bias
      if ($size_dep) {
        my $id  = $c_struct->{'ids'}[$i];
        my $seq = $self->database_get_seq($id);
        my $len = $seq->length;
        $proba *= $len;
      }
    }

    push @$probas, $proba;
    $totproba += $proba;
  }

  # Normalize if necessary
  if ($totproba != 1) {
    $probas = normalize($probas, $totproba);
  }

  return $probas;
}


sub proba_cumul {
  # Put the probas end to end on a line and generate their start position on the
  # line (cumulative distribution). This will help with picking genomes or 
  # nucleotides at random using the rand_weighted() subroutine.
  my ($self, $probas) = @_;
  my $sum = 0;
  return [ 0, map { $sum += $_ } @$probas ];
}


sub rand_weighted {
  # Pick a random number based on the given cumulative probabilities.
  # Cumulative weights can be obtained from the proba_cumul() subroutine.
  my ($cum_probas, $pick, $index) = (shift, rand, -1);
  map { $pick >= $_ ? $index++ : return $index } @$cum_probas;
}


sub rand_seq {
  # Choose a sequence object randomly using a probability distribution
  my ($self, $positions, $oids) = @_;
  return $self->database_get_seq( $$oids[rand_weighted($positions)] ); 
}


sub rand_seq_chimera {
  my ($self, $sequence, $chimera_perc, $positions, $oids) = @_;
  # Produce an amplicon that is a chimera of multiple sequences

  my $chimera;
  # Sanity check
  if ( (scalar @$oids < 2) && ($chimera_perc > 0) ) {
    die "Error: Not enough sequences to produce chimeras\n";
  }
  # Fate now decides to produce a chimera or not
  if ( rand(100) <= $chimera_perc ) {

    # Pick multimera size
    my $m = $self->rand_chimera_size();

    # Pick chimera fragments
    my @pos;
    if ($self->{chimera_kmer}) {
      @pos = $self->kmer_chimera_fragments($m);
    } else {
      # TODO: try to not provide $positions and $oids
      @pos = $self->rand_chimera_fragments($m, $sequence, $positions, $oids);
    }

    # Join chimera fragments
    $chimera = assemble_chimera(@pos);

  } else {
    # No chimera needed
    $chimera = $sequence;
  }
  return $chimera;
}


sub rand_chimera_size {
  # Decide of the number of sequences that the chimera will have, based on the
  # user-defined chimera distribution
  my ($self) = @_;
  return rand_weighted( $self->{chimera_dist_cdf} ) + 2;
}


sub kmer_chimera_fragments {
  # Return a kmer-based chimera of the required size. It is impossible to
  # randomly make one that will meet the required size. So, make multiple
  # attempts and save failed attempts in a pool for later reuse.
  my ($self, $m) = @_;

  my $frags;
  my $pool = $self->{chimera_kmer_pool}->{$m};
  if ( (defined $pool) && (scalar @$pool > 0) ) {
    # Pick a chimera from the pool if possible
    $frags = shift @$pool;

  } else {
    # Attempt multiple times to generate a suitable chimera
    my $actual_m = 0;
    my $nof_tries = 0;
    my $max_nof_tries = 100;
    while ( ($actual_m < $m) && ($nof_tries <= $max_nof_tries) ) {
      $nof_tries++;
      $frags = [ $self->kmer_chimera_fragments_backend($m) ];
      my $actual_m = scalar @$frags / 3;

      if ($nof_tries >= $max_nof_tries) {
        # Could not make a suitable chimera, accept the current chimera
        warn "Warning: Could not make a chimera of $m sequences after ".
          "$max_nof_tries attempts. Accepting a chimera of $actual_m sequences".
          " instead...\n";
        $actual_m = $m;
      }

      if ($actual_m < $m) {
        # Add unsuitable chimera to the pool
        $pool = $self->{chimera_kmer_pool}->{$actual_m};
        push @$pool, $frags;
        # Prevent the pool from growing too big
        my $max_pool_size = 100;
        shift @$pool if scalar @$pool > $max_pool_size;

      } else {
        # We got a suitable chimera... done
        last;
      }
    }
  }

  return @$frags;
}


sub kmer_chimera_fragments_backend {
  # Pick sequence fragments for multimeras where breakpoints are located on
  # shared kmers. A smaller chimera than requested may be returned.
  my ($self, $m) = @_;

  # Initial pair of fragments
  my @pos = $self->rand_kmer_chimera_initial();

  # Append sequence to chimera
  for my $i (3 .. $m) {
    my ($seqid1, $start1, $end1, $seqid2, $start2, $end2) =
      $self->rand_kmer_chimera_extend($pos[-3], $pos[-2], $pos[-1]);

    if (not defined $seqid2) {
      # Could not find a sequence that shared a suitable kmer      
      last;
    }

    @pos[-3..-1] = ($seqid1, $start1, $end1);
    push @pos, ($seqid2, $start2, $end2);
  }

  # Put sequence objects instead of sequence IDs
  for (my $i = 0; $i < scalar @pos; $i = $i+3) {
     my $seqid = $pos[$i];
     my $seq   = $self->database_get_seq($seqid);
     $pos[$i]  = $seq;
  }

  return @pos;
}


sub rand_kmer_chimera_extend {
  # Pick another fragment to add to a kmer-based chimera. Return undef if none
  # can be found
  my ($self, $seqid1, $start1, $end1) = @_;  
  my ($seqid2, $start2, $end2);

  # Get kmer frequencies in the end part of sequence 1
  my ($kmer_arr, $freqs) = $self->{chimera_kmer_col}->counts($seqid1, $start1, 1);

  if (defined $kmer_arr) {

    # Pick a random kmer
    my $kmer_cdf = $self->proba_cumul($freqs);
    my $kmer = $self->rand_kmer_from_collection($kmer_arr, $kmer_cdf);

    # Get a sequence that has the same kmer as the first but is not the first
    $seqid2 = $self->rand_seq_with_kmer( $kmer, $seqid1 );

    # Pick a suitable kmer start on that sequence
    if (defined $seqid2) {

      # Pick a random breakpoint
      # TODO: can we prefer a position not too crazy?
      my $pos1 = $self->rand_kmer_start( $kmer, $seqid1, $start1 );
      my $pos2 = $self->rand_kmer_start( $kmer, $seqid2 );

      # Place breakpoint about the middle of the kmer (kmers are at least 2 bp long) 
      my $middle = int($self->{chimera_kmer} / 2);
      #$start1 = $start1;
      $end1    = $pos1 + $middle - 1;
      $start2  = $pos2 + $middle;
      $end2    = $self->database_get_seq($seqid2)->length;

    }
  }

  return $seqid1, $start1, $end1, $seqid2, $start2, $end2;
}


sub rand_kmer_chimera_initial {
  # Pick two sequences and start points to assemble a kmer-based bimera.
  # An optional starting sequence can be provided.
  my ($self, $seqid1) = @_;

  my $kmer;
  if (defined $seqid1) {
    # Try to pick a kmer from the requested sequence
    $kmer = $self->rand_kmer_of_seq( $seqid1 );
    if (not defined $kmer) {
      die "Error: Sequence $seqid1 did not contain a suitable kmer\n";
    }
  } else {
    # Pick a random kmer and sequence containing that kmer
    $kmer   = $self->rand_kmer_from_collection();
    $seqid1 = $self->rand_seq_with_kmer( $kmer );
  }

  # Get a sequence that has the same kmer as the first but is not the first
  my $seqid2 = $self->rand_seq_with_kmer( $kmer, $seqid1 );
  if (not defined $seqid2) {
    die "Error: Could not find another sequence that contains kmer $kmer\n";
  }

  # Pick random breakpoint positions
  my $pos1 = $self->rand_kmer_start( $kmer, $seqid1 );
  my $pos2 = $self->rand_kmer_start( $kmer, $seqid2 );

  # Swap sequences so that pos1 < pos2
  if ($pos1 > $pos2) {
    ($seqid1, $seqid2) = ($seqid2, $seqid1);
    ($pos1, $pos2) = ($pos2, $pos1);
  }

  # Place breakpoint about the middle of the kmer (kmers are at least 2 bp long) 
  my $middle = int($self->{chimera_kmer} / 2);
  my $start1 = 1;
  my $end1   = $pos1 + $middle - 1;
  my $start2 = $pos2 + $middle;
  my $end2   = $self->database_get_seq($seqid2)->length;

  return $seqid1, $start1, $end1, $seqid2, $start2, $end2;
}


sub rand_kmer_from_collection {
  # Pick a kmer at random amongst all possible kmers in the collection
  my ($self, $kmer_arr, $kmer_cdf) = @_;
  my $kmers = defined $kmer_arr ? $kmer_arr : $self->{chimera_kmer_arr};
  my $cdf   = defined $kmer_cdf ? $kmer_cdf : $self->{chimera_kmer_cdf};
  my $kmer  = $$kmers[rand_weighted($cdf)];
  return $kmer;
}


sub rand_seq_with_kmer {
   # Pick a random sequence ID that contains the given kmer. An optional sequence
   # ID to exclude can be provided.
   my ($self, $kmer, $excl) = @_;
   my $source;
   my ($sources, $freqs) = $self->{chimera_kmer_col}->sources($kmer, $excl, 1);

   my $num_sources = scalar @$sources;
   if ($num_sources > 0) {
     my $cdf = $self->proba_cumul($freqs);
     $source = $$sources[rand_weighted($cdf)];
   }

   return $source;
}


sub rand_kmer_of_seq {
  # Pick a kmer amongst the possible kmers of the given sequence
  my ($self, $seqid) = @_;
  my $kmer;
  my ($kmers, $freqs) = $self->{chimera_kmer_col}->kmers($seqid, 1);
  if (scalar @$kmers > 0) {
    my $cdf = $self->proba_cumul($freqs);
    $kmer = $$kmers[rand_weighted($cdf)];
  }
  return $kmer;
}


sub rand_kmer_start {
  # Pick a kmer starting position at random for the given kmer and sequence ID.
  # An optional minimum start position can be given.
  my ($self, $kmer, $source, $min_start) = @_;
  my $start;
  $min_start ||= 1;
  my $kmer_col = $self->{chimera_kmer_col};
  my $kmer_starts = $kmer_col->positions($kmer, $source);

  # Find index of first index min_idx where position respects min_start
  my $min_idx;
  for (my $i = 0; $i < scalar @$kmer_starts; $i++) {
    my $start = $kmer_starts->[$i];
    if ($start >= $min_start) {
      $min_idx = $i;
      last;
    }
  }

  if (defined $min_idx) {
    # Get a random index between min_idx and the end of the array
    my $rand_idx = $min_idx + int rand (scalar @$kmer_starts - $min_idx);
    # Get the value for this random index
    $start = $kmer_starts->[ $rand_idx ];
  }

  return $start;
}


sub rand_chimera_fragments {
  # Pick which sequences and breakpoints to use to form a chimera
  my ($self, $m, $sequence, $positions, $oids) = @_;

  # Pick random sequences
  my @seqs = ($sequence);
  my $min_len = $sequence->length;
  for (my $i = 2; $i <= $m; $i++) {
    my $prev_seq = $seqs[-1];
    my $seq;
    do {
      $seq = $self->rand_seq($positions, $oids);
    } while ($seq->seq->id eq $prev_seq->seq->id);
    push @seqs, $seq;
    my $seq_len = $seq->length;
    if ( (not defined $min_len) || ($seq_len < $min_len) ) {
      $min_len = $seq_len;
    }
  }

  # Pick random breakpoints
  my $nof_breaks = $m - 1;
  my %breaks = ();
  while ( scalar keys %breaks < $nof_breaks ) {
    # pick a random break
    my $rand_pos = 1 + int( rand($min_len - 1) );
    $breaks{$rand_pos} = undef;
  }
  my @breaks = (1, sort {$a <=> $b} (keys %breaks));
  undef %breaks;

  # Assemble the positional array
  my @pos;
  for (my $i = 1; $i <= $m; $i++) {
    my $seq   = $seqs[$i-1];
    my $start = shift @breaks;
    my $end   = $breaks[0] || $seq->length;
    $breaks[0]++;
    push @pos, ($seq, $start, $end);
  }

  return @pos;
}


sub assemble_chimera {
  # Create a chimera sequence object based on positional information:
  #   seq1, start1, end1, seq2, start2, end2, ...
  my (@pos) = @_;

  # Create the ID, sequence and split location
  my ($chimera_id, $chimera_seq);
  my $chimera_loc = Bio::Location::Split->new();
  while ( my ($seq, $start, $end) = splice @pos, 0, 3 ) {

    # Add amplicon position
    $chimera_loc->add_sub_Location( $seq->location );

    # Add amplicon ID
    if (defined $chimera_id) {
      $chimera_id .= ',';
    }

    # Add subsequence
    my $chimera = $seq->seq;
    $chimera_id .= $chimera->id;
    $chimera_seq .= $chimera->subseq($start, $end);

  }

  # Create a sequence object
  my $chimera = Bio::SeqFeature::SubSeq->new(
    -seq => Bio::PrimarySeq->new( -id => $chimera_id, -seq => $chimera_seq ),
  );

  # Save split location object (a bit hackish)
  $chimera->{_chimera} = $chimera_loc;

  return $chimera;
}


sub rand_seq_orientation {
  # Return a random read orientation: 1 for uncomplemented, or -1 for complemented
  return int(rand()+0.5) ? 1 : -1;
}


sub rand_seq_errors {
  # Introduce sequencing errors (point mutations, homopolymers) in a sequence
  # based on error models
  my ($self, $seq) = @_;
  my $seq_str = $seq->seq();
  my $error_specs = {}; # Error specifications

  # First, specify errors in homopolymeric stretches
  $error_specs = $self->rand_homopolymer_errors($seq_str, $error_specs)
    if $self->{homopolymer_dist};

  # Then, specify point sequencing errors: substitutions, insertions, deletions
  $error_specs = $self->rand_point_errors($seq_str, $error_specs) 
    if $self->{mutation_para1};

  # Finally, actually implement the errors as per the specifications
  $seq->errors($error_specs) if (scalar keys %$error_specs > 0);

  return $seq;
}


sub rand_homopolymer_errors {
  # Specify sequencing errors in a sequence's homopolymeric stretches
  my ($self, $seq_str, $error_specs) = @_;
  while ( $seq_str =~ m/(.)(\1+)/g ) {

    # Found a homopolymer
    my $res = $1;                       # residue in homopolymer
    my $len = length($2) + 1;           # length of the homopolymer
    my $pos = pos($seq_str) - $len + 1; # start of the homopolymer (residue no.)

    # Apply homopolymer model based on normal distribution N(mean, standard deviation)
    #   Balzer:    N(n, 0.03494 + n * 0.06856)  Balzer et al. 2010
    #   Richter:   N(n, 0.15 * sqrt(n))         Richter et al. 2008
    #   Margulies: N(n, 0.15 * n)               Margulies et al. 2005
    my ($stddev, $new_len, $diff) = (0, 0, 0);
    if ( $self->{homopolymer_dist} eq 'balzer' ) {
      $stddev = 0.03494 + $len * 0.06856;
    } elsif ($self->{homopolymer_dist} eq 'richter') {
      $stddev = 0.15 * sqrt($len);
    } elsif ($self->{homopolymer_dist} eq 'margulies') {
      $stddev = 0.15 * $len;
    } else {
      die "Error: Unknown homopolymer distribution '".$self->{homopolymer_dist}."'\n";
    }
    $new_len = int( $len + $stddev * randn() + 0.5 );
    $new_len = 0 if $new_len < 0;
    # We're done if no error was introduced
    $diff = $new_len - $len;
    next unless $diff;
    # Otherwise, track the error generated
    if ($diff > 0) { 
      # Homopolymer extension
      push @{$$error_specs{$pos}{'+'}}, ($res) x $diff;
    } elsif ($diff < 0) {
      # Homopolymer shrinkage
      for my $offset ( 0 .. abs($diff)-1 ) {
        push @{$$error_specs{$pos+$offset}{'-'}}, undef;
      }
    }

  }

  return $error_specs;
}


sub rand_point_errors {
  # Do some random point sequencing errors on a sequence based on a model
  my ($self, $seq_str, $error_specs) = @_;

  # Mutation cumulative density functions (cdf) for this sequence length
  my $seq_len = length $seq_str;
  if ( not defined $self->{mutation_cdf}->{$seq_len} ) {
    my $mut_pdf  = []; # probability density function
    my $mut_freq =  0; # average
    my $mut_sum  =  0;

    if ($self->{mutation_model} eq 'uniform') {
      # Uniform error model
      # para1 is the average mutation frequency
      my $proba = 1 / $seq_len;
      $mut_pdf  = [ map { $proba } (1 .. $seq_len) ];
      $mut_freq = $self->{mutation_para1};
      $mut_sum  = 1;

    } elsif ($self->{mutation_model} eq 'linear') {
      # Linear error model
      # para 1 is the error rate at the 5' end of the read
      # para 2 is the error rate at the 3' end
      $mut_freq = abs( $self->{mutation_para2} + $self->{mutation_para1} ) / 2;
      if ($seq_len == 1) {
        $$mut_pdf[0] = $mut_freq;
        $mut_sum = $mut_freq
      } elsif ($seq_len > 1) {
        my $slope = ($self->{mutation_para2} - $self->{mutation_para1}) / ($seq_len-1);
        for my $i (0 .. $seq_len-1) {
          my $val = $self->{mutation_para1} + $i * $slope;
          $mut_sum += $val;
          $$mut_pdf[$i] = $val;
        }
      }
      
    } elsif ($self->{mutation_model} eq 'poly4') {
      # Fourth degree polynomial error model: e = para1 + para2 * i**4
      for my $i (0 .. $seq_len-1) {
        my $val = $self->{mutation_para1} + $self->{mutation_para2} * ($i+1)**4;
        $mut_sum += $val;
        $$mut_pdf[$i] = $val;
      }
      $mut_freq = $mut_sum / $seq_len;

    } else {
      die "Error: '".$self->{mutation_model}."' is not a supported error distribution\n";
    }

    # Normalize to 1 if needed
    if ($mut_sum != 1) {
      $mut_pdf = normalize($mut_pdf, $mut_sum);
    }

    # TODO: Could have sanity checks so that mut_pdf should have no values < 0 or > 100

    $self->{mutation_cdf}->{$seq_len} = $self->proba_cumul($mut_pdf);
    $self->{mutation_avg}->{$seq_len} = $mut_freq;
  }

  my $mut_cdf = $self->{mutation_cdf}->{$seq_len};
  my $mut_avg = $self->{mutation_avg}->{$seq_len};

  # Number of mutations to make in this sequence is assumed to follow a Normal
  # distribution N( mutation_freq, 0.3 * mutation_freq )
  my $read_mutation_freq = $mut_avg + 0.3 * $mut_avg * randn();
  my $nof_mutations = $seq_len * $read_mutation_freq / 100;
  my $int_part = int $nof_mutations;
  my $dec_part = rand(1) < ($nof_mutations - $int_part);
  $nof_mutations = $int_part + $dec_part;

  # Exit without doing anything if there are no mutations to do
  return $error_specs if $nof_mutations == 0;

  # Make as many mutations in read as needed based on model
  my $subst_frac = $self->{mutation_ratio}->[0] / 100;
  for ( 1 .. $nof_mutations ) {

    # Position to mutate
    my $idx = rand_weighted( $mut_cdf );

    # Do a substitution or indel
    if ( rand() <= $subst_frac ) {

      # Substitute at given position by a random replacement nucleotide
      push @{$$error_specs{$idx+1}{'%'}}, $self->rand_res( substr($seq_str, $idx, 1) );
    } else {

      # Equiprobably insert or delete
      if ( rand() < 0.5 ) {
        # Insertion after given position
        push @{$$error_specs{$idx+1}{'+'}}, $self->rand_res(); 
      } else {
        # Make a deletion at given position
        next if length($seq_str) == 1; # skip this deletion to avoid a 0 length
        push @{$$error_specs{$idx+1}{'-'}}, undef;
      }

    }

  }

  return $error_specs;
}


sub rand_res {
  # Pick a residue at random from the stored alphabet (dna, rna or protein).
  # An optional residue to exclude can be specified.
  my ($self, $not_nuc) = @_;
  my $cdf;
  my @res;
  if (not defined $not_nuc) {
    # Use complete alphabet
    @res = keys %{$self->{alphabet_hash}};
    $cdf = $self->{alphabet_complete_cdf};
  } else {
    # Remove non-desired residue from alphabet
    my %res = %{$self->{alphabet_hash}};
    delete $res{uc($not_nuc)};
    @res = keys %res;
    $cdf = $self->{alphabet_truncated_cdf};
  }
  my $res = $res[rand_weighted($cdf)];
  return $res; 
}


sub rand_seq_length {
  # Choose the sequence length following a given probability distribution
  my($avg, $model, $stddev) = @_;
  my $length;
  if (not $model) {
    # No specified distribution: all the sequences have the length of the average
    $length = $avg;
  } else {
    if ($model eq 'uniform') {
      # Uniform distribution: integers uniformly distributed in [min, max]
      my ($min, $max) = ($avg - $stddev, $avg + $stddev);
      $length = $min + int( rand( $max - $min + 1 ) );
    } elsif ($model eq 'normal') {
      # Gaussian distribution: decimal number normally distribution in N(avg,stddev)
      $length = $avg + $stddev * randn();
      $length = int( $length + 0.5 );
    } else {
      die "Error: '$model' is not a supported read or insert length distribution\n";
    }
  }
  $length = 1 if ($length < 1);
  return $length;
}


sub rand_seq_pos {
  # Pick the coordinates (start and end) of an amplicon or random shotgun read.
  # Coordinate system: the first base is 1 and the number is inclusive, ie 1-2
  # are the first two bases of the sequence
  my ($seq_obj, $read_length, $amplicon, $mid) = @_;
  # Read length includes the MID
  my $length = $read_length - length($mid);
  # Pick starting position
  my $start;
  if (defined $amplicon) {
    # Amplicon always start at first position of amplicon
    $start = 1;
  } else {
    # Shotgun reads start at a random position in genome
    $start = int( rand($seq_obj->length - $length + 1) ) + 1;
  }
  # End position
  my $end = $start + $length - 1;
  return $start, $end;
}


sub randn {
  # Normally distributed random value (mean 0 and standard deviation 1) using
  # the Box-Mueller transformation method, adapted from the Perl Cookbook
  my ($g1, $g2, $w);
  do {
    $g1 = 2 * rand() - 1; # uniformly distributed
    $g2 = 2 * rand() - 1;
    $w = $g1**2 + $g2**2; # variance
  } while ( $w >= 1 );
  $w = sqrt( (-2 * log($w)) / $w ); # weight
  $g1 *= $w; # gaussian-distributed
  if ( wantarray ) {
    $g2 *= $w;
    return ($g1, $g2);
  } else {
    return $g1;
  }
}


sub randig {
   # Random value sampled from the inverse gaussian (a.k.a. Wald) distribution,
   # using the method at http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
   my ($mu, $lambda) = @_;
   my $y = randn()**2;
   my $x = $mu + ($mu**2 * $y)/(2 * $lambda) - $mu / (2 * $lambda)
           * sqrt(4 * $mu * $lambda * $y + $mu**2 * $y**2);
   if ( rand() <= $mu / ($mu + $x) ) {
      $y = $x;
   } else {
      $y = $mu**2 / $x;
   }
   return $y;
}


sub randomize {
  # Randomize an array using the Fisher-Yates shuffle described in the Perl
  # cookbook.
  my ($array) = @_;
  my $i;
  for ($i = @$array; --$i; ) {
   my $j = int rand($i+1);
   next if $i == $j;
   @$array[$i,$j] = @$array[$j,$i];
  }
  return $array;
}


sub database_create {
  # Read and import sequences
  # Parameters:
  #   * FASTA file containing the sequences or '-' for stdin. REQUIRED
  #   * Sequencing unidirectionally? 0: no, 1: yes forward, -1: yes reverse
  #   * Amplicon PCR primers (optional): Should be provided in a FASTA file and
  #     use the IUPAC convention. If a primer sequence is given, any sequence
  #     that does not contain the primer (or its reverse complement for the
  #     reverse primer) is skipped, while any sequence that matches is trimmed
  #     so that it is flush with the primer sequence
  #   * Abundance file (optional): To avoid registering sequences in the database
  #     unless they are needed
  #   * Delete chars (optional): Characters to delete form the sequences.
  #   * Minimum sequence size: Skip sequences smaller than that
  my ($self, $fasta_file, $unidirectional, $forward_reverse_primers,
    $abundance_file, $delete_chars, $min_len) = @_;
  $min_len = 1 if not defined $min_len;
  # Input filehandle
  if (not defined $fasta_file) {
    die "Error: No reference sequences provided\n";
  }
  my $in;
  if ($fasta_file eq '-') {
    $in = Bio::SeqIO->newFh(
      -fh     => \*STDIN,
      -format => 'fasta',
    );
  } else {
    $in = Bio::SeqIO->newFh(
      -file   => $fasta_file,
      -format => 'fasta',
    );
  }

  # Get list of all IDs with a manually-specified abundance
  my %ids_to_keep;
  if ($abundance_file) {
    my ($ids) = community_read_abundances($abundance_file);
    for my $comm_num (0 .. $#$ids) {
      for my $gen_num ( 0 .. scalar @{$$ids[$comm_num]} - 1 ) {
        my $id = $$ids[$comm_num][$gen_num];
        $ids_to_keep{$id} = undef;
      }
    }
  }

  # Initialize search for amplicons
  my $amplicon_search;
  if (defined $forward_reverse_primers) {
    $amplicon_search = Bio::Tools::AmpliconSearch->new(
      -primer_file => $forward_reverse_primers,
    );
  }

  # Process database sequences
  my %seq_db;      # hash of BioPerl sequence objects (all amplicons)
  my %seq_ids;     # hash of reference sequence IDs and IDs of their amplicons
  my %mol_types;   # hash of count of molecule types (dna, rna, protein)
  while ( my $ref_seq = <$in> ) {
    # Skip empty sequences
    next if not $ref_seq->seq;
    # Record molecule type
    $mol_types{$ref_seq->alphabet}++;
    # Skip unwanted sequences
    my $ref_seq_id = $ref_seq->id;
    next if (scalar keys %ids_to_keep > 0) && (not exists $ids_to_keep{$ref_seq_id});
    # If we are sequencing from the reverse strand, reverse complement now
    if ($unidirectional == -1) {
      $ref_seq = $ref_seq->revcom;
    }

    # Extract amplicons if needed
    my $amp_seqs;
    if (defined $amplicon_search) {
      $amplicon_search->template($ref_seq);
      while (my $amp_seq = $amplicon_search->next_amplicon) {
        push @$amp_seqs, $amp_seq;
      }
      next if not defined $amp_seqs;
    } else {
      $amp_seqs = [ Bio::SeqFeature::SubSeq->new( -start    => 1,
                                                  -end      => $ref_seq->length,
                                                  -template => $ref_seq,    ) ];
    }

    for my $amp_seq (@$amp_seqs) {
      # Remove forbidden chars
      if ( (defined $delete_chars) && (not $delete_chars eq '') ) {

        ### TODO: Use Bio::Location::Split here as well?
        my $clean_seq = $amp_seq->seq;
        my $clean_seqstr = $clean_seq->seq;
        my $dirty_length = length $clean_seqstr;

        $clean_seqstr =~ s/[$delete_chars]//gi;
        my $num_dels = $dirty_length - length $clean_seqstr;
        if ($num_dels > 0) {
          # Update sequence with cleaned sequence string
          $clean_seq->seq($clean_seqstr);
          $amp_seq->seq($clean_seq);
          # Adjust (decrease) end of feature
          $amp_seq->end( $amp_seq->end - $num_dels );
        }

      }
      # Skip the sequence if it is too small
      next if $amp_seq->length < $min_len;

      # Save amplicon sequence and create a barcode that identifies it
      my $amp_bc = create_amp_barcode($amp_seq, $ref_seq_id);
      $seq_db{$amp_bc} = $amp_seq;
      $seq_ids{$ref_seq_id}{$amp_bc} = undef;

    }

  }
  undef $in; # close the filehandle (maybe?!)

  # Error if no usable sequences in the database
  if (scalar keys %seq_ids == 0) {
    die "Error: No genome sequences could be used. If you specified a file of".
      " abundances for the genome sequences, make sure that their ID match the".
      " ID in the FASTA file. If you specified amplicon primers, verify that ".
      "they match some genome sequences.\n";
  }

  # Determine database type: dna, rna, protein
  my $db_alphabet = $self->database_get_mol_type(\%mol_types);
  $self->{alphabet} = $db_alphabet;

  # Error if using amplicon on protein database
  if ( ($db_alphabet eq 'protein') && (defined $forward_reverse_primers) ) {
    die "Error: Cannot use amplicon primers with proteic reference sequences\n";
  }

  # Error if using wrong direction on protein database
  if ( ($db_alphabet eq 'protein') && ($unidirectional != 1) ) {
    die "Error: Got <unidirectional> = $unidirectional but can only use ".
      "<unidirectional> = 1 with proteic reference sequences\n";
  }

  my $database = { 'db' => \%seq_db, 'ids' => \%seq_ids };
  return $database;
}


sub create_amp_barcode {
  # Create a barcode that is unique for each amplicon, store it and return it
  my ($amp_sf, $ref_seq_id) = @_;
  my $sep = '/';
  my @elems = ($ref_seq_id, $amp_sf->start, $amp_sf->end, $amp_sf->strand || 1);

  #### TODO: follow the spec: id:start..end/strand

  my $barcode = join $sep, @elems;
  $amp_sf->{_barcode} = $barcode;
  return $barcode;
}

sub get_amp_barcode {
  # Get the amplicon barcode
  my ($amp_sf) = @_;
  return $amp_sf->{_barcode};
}


sub database_get_mol_type {
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
    die "Error: Cannot determine what type of molecules the reference sequences".
        " are. Got $max_count sequences of type '$max_type' and $other_count ".
        "others.\n";
  }
  if ( (not $max_type eq 'dna') && (not $max_type eq 'rna') && (not $max_type eq 'protein') ) {
    die "Error: Reference sequences are in an unknown alphabet '$max_type'\n";
  }
  return $max_type;
}


sub database_get_all_oids {
  # Retrieve all object IDs from the database. These OIDs match the output of
  # the database_get_all_seqs method.
  my ($self) = @_;
  my @oids;
  while ( my ($oid, undef) = each %{$self->{database}->{db}} ) {
    push @oids, $oid;
  }
  return \@oids;
}


sub database_get_all_seqs {
  # Retrieve all sequence objects from the database. These sequence objects match
  # the output of the database_get_all_oids method.
  my ($self)  = @_;
  my @seqs;
  while ( my (undef, $seq) = each %{$self->{database}->{db}} ) {
    push @seqs, $seq;
  }
  return \@seqs;
}


sub database_get_seq {
  # Retrieve a sequence object from the database based on its object ID
  my ($self, $oid)  = @_;
  my $db = $self->{database}->{db};
  my $seq_obj;
  if (not exists $$db{$oid}) {
    warn "Warning: Could not find sequence with object ID '$oid' in the database\n";
  }
  $seq_obj = $$db{$oid};
  return $seq_obj;
}


sub database_get_children_seq {
  # Retrieve all the sequences object made from a reference sequence based on the
  # ID of the reference sequence
  my ($self, $refseqid)  = @_;
  my @children;
  for my $child_oid ( keys %{$self->{database}->{ids}->{$refseqid}} ) {
    push @children, $self->database_get_seq($child_oid);
  }
  return \@children;
}

sub database_get_parent_id {
  # Based on a sequence object ID, retrieve the ID of the reference sequence it
  # came from
  my ($self, $oid) = @_;
  my $seq_id = $self->database_get_seq($oid)->seq->id;
  return $seq_id;
}


sub iupac_to_regexp {
  # Create a regular expression to match a nucleotide sequence that contain
  # degeneracies (in IUPAC standard)
  my ($seq) = @_;
  # Basic IUPAC code
  #my %iupac = (
  #  'A' => ['A'],
  #  'C' => ['C'],
  #  'G' => ['G'],
  #  'T' => ['T'],
  #  'U' => ['U'],
  #  'R' => ['G', 'A'],
  #  'Y' => ['T', 'C'],
  #  'K' => ['G', 'T'],
  #  'M' => ['A', 'C'],
  #  'S' => ['G', 'C'],
  #  'W' => ['A', 'T'],
  #  'B' => ['G', 'T', 'C'],
  #  'D' => ['G', 'A', 'T'],
  #  'H' => ['A', 'C', 'T'],
  #  'V' => ['G', 'C', 'A'],
  #  'N' => ['A', 'G', 'C', 'T'],
  #);
  # IUPAC code
  #   + degenerate primer residues matching ambiguous template residues
  #   + degenerate primer residues matching uracil U
  my %iupac = (
    'A' => ['A'],
    'C' => ['C'],
    'G' => ['G'],
    'T' => ['T'],
    'U' => ['U'],
    'R' => ['G', 'A', 'R'],
    'Y' => ['T', 'U', 'C', 'Y'],
    'K' => ['G', 'T', 'U', 'K'],
    'M' => ['A', 'C', 'M'],
    'S' => ['G', 'C', 'S'],
    'W' => ['A', 'T', 'U', 'W'],
    'B' => ['G', 'T', 'U', 'C', 'Y', 'K', 'S', 'B'],
    'D' => ['G', 'A', 'T', 'U', 'R', 'K', 'W', 'D'],
    'H' => ['A', 'C', 'T', 'U', 'Y', 'M', 'W', 'H'],
    'V' => ['G', 'C', 'A', 'R', 'M', 'S', 'V'],
    'N' => ['A', 'G', 'C', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'],
  );
  # Regular expression to catch this sequence
  my $regexp;
  for my $pos (0 .. length($seq)-1) {
    my $res = substr $seq, $pos, 1;
    my $iupacs = $iupac{$res};
    if (not defined $iupacs) {
      die "Error: Primer sequence '$seq' is not a valid IUPAC sequence. ".
        "Offending character is '$res'.\n";
    }
    if (scalar @$iupacs > 1) {
      $regexp .= '['.join('',@$iupacs).']';
    } else {
      $regexp .= $$iupacs[0];
    }
  }
  $regexp = qr/$regexp/i;
  return $regexp;
}


sub lib_coverage {
  # Calculate number of sequences needed to reach a given coverage. If the
  # number of sequences is provided, calculate the coverage
  my ($self, $c_struct) = @_;
  my $coverage    = $self->{coverage_fold};
  my $nof_seqs    = $self->{total_reads};
  my $read_length = $self->{read_length};
  # 1/ Calculate library length and size
  my $ref_ids    = $c_struct->{'ids'};
  my $diversity  = scalar @$ref_ids;
  my $lib_length = 0;
  for my $ref_id (@$ref_ids) {
    my $seqobj = $self->database_get_seq($ref_id);
    my $seqlen = $seqobj->length;
    $lib_length += $seqlen;
  }
  # 2/ Calculate number of sequences to generate based on desired coverage. If
  # both number of reads and coverage fold were given, coverage has precedence.
  if ($coverage) {
    $nof_seqs = ($coverage * $lib_length) / $read_length;
    if ( int($nof_seqs) < $nof_seqs ){
      $nof_seqs = int($nof_seqs + 1); # ceiling
    }
  }
  $coverage = ($nof_seqs * $read_length) / $lib_length;
  # 3/ Sanity check

  # TODO: Warn only if diversity was explicitely specified on the command line

  if ( $nof_seqs < $diversity) {
    warn "Warning: The number of reads to produce is lower than the required ".
      "diversity. Increase the coverage or number of reads to achieve this ".
      "diversity.\n";
    $self->{diversity}->[$self->{cur_lib}-1] = $nof_seqs;
  }
  return $nof_seqs, $coverage;
}


sub new_subseq {
  # Create a new sequence object as a subsequence of another one and name it so
  # we can trace back where it came from
  my ($fragnum, $seq_feat, $unidirectional, $orientation, $start, $end, $mid,
    $mate_number, $lib_number, $tracking, $qual_levels) = @_;

  # If the length is too short for this read, no choice but to decrease it.
  $start = 1 if $start < 1;
  $end   = $seq_feat->length if $end > $seq_feat->length;

  # Build the sequence ID
  my $name_sep  = '_';
  my $field_sep = ' ';
  my $mate_sep  = '/'; # mate pair indicator, by convention
  my $newid = $fragnum;
  if (defined $lib_number) {
    $newid = $lib_number.$name_sep.$newid;
  }
  if (defined $mate_number) {
    $newid .= $mate_sep.$mate_number;
  }

  # Create a new simulated read object
  my $newseq = Bio::Seq::SimulatedRead->new(
     -id          => $newid,
     -reference   => $seq_feat->seq,
     -start       => $start,
     -end         => $end,
     -strand      => $orientation,
     -mid         => $mid,
     -track       => $tracking,
     -coord_style => 'genbank',
     -qual_levels => $qual_levels,
  );

  # Record location of amplicon on reference sequence in the sequence description
  if ( $seq_feat->isa('Bio::SeqFeature::Amplicon') || exists($seq_feat->{_chimera}) ) {
    my $amplicon_desc = gen_subseq_desc($seq_feat);
    my $desc = $newseq->desc;
    $desc =~ s/(reference=\S+)/$1 $amplicon_desc/;
    $newseq->desc($desc);
  }

  # Database sequences were already reverse-complemented if reverse sequencing
  # was requested
  if ($unidirectional == -1) {
    $orientation *= -1;
    $newseq = set_read_orientation($newseq, $orientation);
  }

  return $newseq;
}


sub gen_subseq_desc {
  my ($seq_feat) = @_;

  # Chimeras have several locations (a Bio::Location::Split object)
  my @locations;
  if (exists $seq_feat->{_chimera}) {
    @locations = $seq_feat->{_chimera}->sub_Location();
  } else {
    @locations = ( $seq_feat->location );
  } 

  for (my $i = 0; $i <= scalar @locations - 1; $i++) {
    my $location = $locations[$i];
    my $strand = $location->strand || 1;
    if ($strand == 1) {
      $location = $location->start.'..'.$location->end;
    } elsif ($strand == -1) {
      $location = 'complement('.$location->start.'..'.$location->end.')';
    } else {
      die "Error: Strand should be -1 or 1, but got '".$location."'\n";
    }
    $locations[$i] = $location;
  }

  my $desc = 'amplicon='.join(',', @locations);
  return $desc;
}


sub set_read_orientation {
  # Set read orientation and change its description accordingly
  my ($seq, $new_orientation) = @_;
  $seq->strand($new_orientation);
  my $desc = $seq->desc;
  $desc =~ s/position=(complement\()?(\d+)\.\.(\d+)(\))?/position=/;
  my ($start, $end) = ($2, $3);
  if ($new_orientation == -1) {
    $desc =~ s/position=/position=complement($start\.\.$end)/;
  } else {
    $desc =~ s/position=/position=$start\.\.$end/;
  }
  $seq->desc( $desc );
  return $seq;
}


sub two_array_sort {
  # Sort 2 arrays by taking the numeric sort of the first one and keeping the 
  # element of the second one match those of the first one
  my ($l1, $l2) = @_;
  my @ids = map { [ $$l1[$_], $$l2[$_] ] } (0..$#$l1);
  @ids = sort { $a->[0] <=> $b->[0] } @ids;
  my @k1;
  my @k2;
  for (my $i = 0; $i < scalar @ids; $i++) {
    $k1[$i] = $ids[$i][0];
    $k2[$i] = $ids[$i][1];
  }
  return \@k1, \@k2;
}


sub normalize {
   # Normalize an arrayref to 1.
   my ($arr, $total) = @_;
   if (not $total) { # total undef or 0
      die "Error: Need to provide a valid total\n";
   }
   $arr = [ map {$_ / $total} @$arr ];
   return $arr;
}


1;

