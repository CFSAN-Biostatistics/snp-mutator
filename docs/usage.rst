========
Usage
========

SNP Mutator is a single script called ``snpmutator``.  After installation,
it will be on your path, ready to use.

Quick Start
-----------

Step 1 - Create a work directory::

    $ mkdir workdir
    $ cd workdir


Step 2 - Gather a reference fasta file::

    # For this example we will use a fasta file from our sister project, snp-pipeline
    $ wget https://raw.githubusercontent.com/CFSAN-Biostatistics/snp-pipeline/master/snppipeline/data/agonaInputs/reference/NC_011149.fasta

Step 3 - Generate the mutated sequences::

    # The mutated sequence files are generated in the current working directory
    # -r 1, set random seed
    # -n 1000, generate 1000 mutated sequences
    # -s 900, nine hundred substitutions in each mutated sequence
    # -i 50, fifty insertion in each mutated sequence
    # -d 50, fifty deletions in each mutated sequence
    # -o summary.tsv, generate a mutation summary file called summary.tsv
    # -v variants.vcf, generate a VCF file of mutations
    # -p 100000, choose mutations from a pool of 100000 positions
    # -g 100, partition the 1000 replicates into 10 groups of 100 replicates, with each group having a separate pool of positions
    # -m, create monomorphic alleles within each pool
    $ snpmutator -r 1 -n 1000 -s 900 -i 50 -d 50 -o summary.tsv -v variants.vcf -p 100000 -g 100 -m NC_011149.fasta

Step 4 - Examine the results::

    $ ls NC_011149_mutated_*.fasta
    $ less summary.tsv


Input Files
-----------
The only input file is the reference fasta file.

Note: Multi-fasta files may produce unexpected results.  All the sequences are concatenated
into a single sequence.  All description lines after the first are discarded.


Output Files
------------

Summary
~~~~~~~
An optional summary file in tab-delimited format lists the positions of the mutations for
each of the replicates with the original base and the resulting mutation at each position.

VCF
~~~
An optional VCF file lists the mutations in Variant Call Format.

Replicate Files
~~~~~~~~~~~~~~~
Multiple mutated replicate files are generated in the current working directory.  Files are
named with the basename of the original reference file, suffixed with ``_mutated_#.fasta``.

For example, if the reference file name is ``NC_011149.fasta``, the first two replicate files
are named ``NC_011149_mutated_1.fasta`` and ``NC_011149_mutated_2.fasta``.

The defline (description) of the generated fasta files is copied from the original reference
fasta file, but with a suffix describing the mutations.  For example, the defline suffix
``(mutated s=2 i=1 d=0)`` indicates there are two substitutions, one insertion, and zero deletions.

Pooling
-------
The ``--pool`` option increases the likelihood of multiple replicates having mutations at the
same positions by limiting the number of positions along the genome where mutations will be
introduced.  The positions in the pool are choosen randomly and uniformly from the positions
in the genome.

Grouping
--------
The ``--group`` option partitions the replicates into groups with each group having a different pool
of eligible positions.  This has the effect of creating more closely related replicates within
groups and more distant replicates between groups.

Monomorphic Alleles
-------------------
The ``--mono`` option ensures that when multiple replicates have a mutation at the same position,
the mutation will be identical in each replicate.  However, when used with the ``--group`` option, the
monomorphic mutations are only within the group.  Different groups of replicates may have polymorphic
alleles with respect to other groups of replicates.


Command Reference
-----------------

::

  usage: snpmutator [-h] [-o FILE] [-n INT] [-s INT] [-i INT] [-d INT] [-r INT]
                    [-p INT] [-g INT] [-m] [-v FILE] [--version]
                    input_fasta_file

  Generate mutated sequence files from a reference genome. Takes a fasta file
  and creates a specified number of randomly generated base substitutions,
  insertions, and deletions. Outputs the mutated genomes, and optionally, a
  summary file listing the mutations by position.

  positional arguments:
    input_fasta_file      Input fasta file.

  optional arguments:
    -h, --help            show this help message and exit
    -o FILE, --summary FILE
                          Output positional summary file. (default: None)
    -n INT, --num-simulations INT
                          Number of mutated sequences to generate. (default:
                          100)
    -s INT, --num-substitutions INT
                          Number of substitutions. (default: 500)
    -i INT, --num-insertions INT
                          Number of insertions. (default: 20)
    -d INT, --num-deletions INT
                          Number of deletions. (default: 20)
    -r INT, --random-seed INT
                          Random number seed; if not set, the results are not
                          reproducible. (default: None)
    -p INT, --pool INT    Choose variants from a pool of eligible positions of
                          the specified size (default: 0)
    -g INT, --group INT   Group size. When greater than zero, this parameter
                          chooses a new pool of positions for each group of
                          replicates. (default: None)
    -m, --mono            Create monomorphic alleles (default: False)
    -v FILE, --vcf FILE   Output VCF file. (default: None)
    --version             show program's version number and exit
