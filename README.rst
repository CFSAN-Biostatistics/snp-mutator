===============================
SNP Mutator
===============================


.. Image showing the PyPI version badge - links to PyPI
.. image:: https://img.shields.io/pypi/v/snp-mutator.svg
        :target: https://pypi.python.org/pypi/snp-mutator

.. Image showing the PyPi download per month count  - links to PyPI
.. .. image:: https://img.shields.io/pypi/dm/snp-mutator.svg
..        :target: https://pypi.python.org/pypi/snp-mutator

.. Image showing the Travis Continuous Integration test status, commented out for now
.. .. image:: https://img.shields.io/travis/CFSAN-Biostatistics/snp-mutator.svg
..        :target: https://travis-ci.org/CFSAN-Biostatistics/snp-mutator



Generate mutated sequence files from a reference genome.

SNP Mutator was developed by the United States Food
and Drug Administration, Center for Food Safety and Applied Nutrition.

* Free software
* Documentation: https://snp-mutator.readthedocs.io/en/latest/readme.html
* Source Code: https://github.com/CFSAN-Biostatistics/snp-mutator
* PyPI Distribution: https://pypi.python.org/pypi/snp-mutator


Features
--------

* Reads a fasta file and generates any number of mutated fasta replicate files.
* Mutations can be any number of single-base substitutions, insertions, and deletions at randomly
  chosen positions, uniformly distributed across the genome.
* Mutations can be chosen from a subset (pool) of all possible positions.
* Replicates can be partitioned into multiple groups with each group sharing a pool of eligible positions.
* Generates a summary file listing the original base and the mutation for all mutated positions.
* Mutations can be either monomorphic or polymorphic.
* VCF file generation.
* Various metrics can be saved to an output file.

Citing SNP Mutator
--------------------------------------

To cite SNP Mutator, please cite the publication below:

    `Davis S, Pettengill JB, Luo Y, Payne J, Shpuntoff A, Rand H, Strain E. (2015)
    CFSAN SNP Pipeline: an automated method for constructing SNP matrices from next-generation sequence data.
    PeerJ Computer Science 1:e20   https://doi.org/10.7717/peerj-cs.20 <https://doi.org/10.7717/peerj-cs.20>`_

License
-------

See the LICENSE file included in the SNP Mutator distribution.
