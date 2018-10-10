.. :changelog:

History
-------

1.0.0 (2018-10-10)
---------------------
* Mutations can be chosen from a subset (pool) of all possible positions.
* Replicates can be partitioned into multiple groups with each group sharing a pool of eligible positions.
* Mutations can be either monomorphic or polymorphic.
* Add the capability to generate VCF output files.
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
