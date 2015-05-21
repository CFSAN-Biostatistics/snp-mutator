========
Usage
========

SNP Mutator is a single python script called ``snpmutator.py``.  After installation,
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
    # -n 2, generate two mutated sequences
    # -s 2, two substitutions in each mutated sequence
    # -i 1, one insertion in each mutated sequence
    # -d 0, zero deletions in each mutated sequence
    # -o summary.tsv, generate a mutation summary file called summary.tsv
    $ snpmutator.py -r 1 -n 2 -s 2 -i 1 -d 0 -o summary.tsv NC_011149.fasta

Step 4 - Examine the results::

    $ head NC_011149_mutated_*.fasta
    $ cat summary.tsv


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

Replicate Files
~~~~~~~~~~~~~~~
Multiple mutated replicate files are generated in the current working directory.  Files are 
named with the basename of the original reference file, suffixed with ``_mutated_#.fasta``.

For example, if the reference file name is ``NC_011149.fasta``, the first two replicate files
are named ``NC_011149_mutated_1.fasta`` and ``NC_011149_mutated_2.fasta``.

The defline (description) of the generated fasta files is copied from the original reference
fasta file, but with a suffix describing the mutations.  For example, the defline suffix 
``(mutated s=2 i=1 d=0)`` indicates there are two substitutions, one insertion, and zero deletions.


Command Reference
-----------------

::

  usage: snpmutator.py [-h] [-o FILE] [-n INT] [-s INT] [-i INT] [-d INT]
                       [-r INT]
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
