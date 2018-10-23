.. :changelog:

History
-------

1.1.0 (2018-10-23)
---------------------
* Added the capability to specify the fasta defline sequence id in the output files.
* The unmutated, but concatenated sequence reference can be written to a separate output file.
  All the replicates will be mutations of this file.

1.0.0 (2018-10-12)
---------------------
* Mutations can be chosen from a subset (pool) of all possible positions.
* Replicates can be partitioned into multiple groups with each group sharing a pool of eligible positions.
* Mutations can be either monomorphic or polymorphic.
* Add the capability to generate VCF output files.
* Add the capability to generate a metrics output file.
* The installer creates an executable script called ``snpmutator``.
* Insertions are now placed after the original reference position, not before, as in prior versions of
  snpmutator.  This means it will not be possible to exactly reproduce the results of prior versions
  of this software.

0.2.0 (2016-01-15)
---------------------

* Allow lowercase bases in the input fasta file.
* Do not mutate gaps or ambiguous positions.
* Add a command line switch to show the program version.

0.1.1 (2015-06-17)
---------------------

* Remove spaces from the summary file column headings.  This will simplify downstream
  analysis in some scripting languages.


0.1.0 (2015-05-21)
---------------------

* First release on GitHub, Read the Docs, and PyPI.
