============
Installation
============

.. highlight:: bash

At the command line::

    $ pip install --user snp-mutator

Update your .bashrc file with the path to user-installed python packages::

    export PATH=~/.local/bin:$PATH

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv snp-mutator
    $ pip install snp-mutator


Upgrading SNP Mutator
-----------------------------------------

If you previously installed with pip, you can upgrade to the newest version from the command line::

    $ pip install --user --upgrade snp-mutator


Uninstalling SNP Mutator
--------------------------------------------

If you installed with pip, you can uninstall from the command line::

    $ pip uninstall snp-mutator


Requirements
------------

SNP Mutator requires requires python version 2.6 or 2.7. It has not been tested on other python versions.
