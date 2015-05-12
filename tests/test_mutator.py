#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_mutator
----------------------------------

Tests for `mutator` module.
"""

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import unittest
from testfixtures import TempDirectory
import random
import argparse
from itertools import izip_longest
from mutator import mutator



def make_random_dna_string(length):
    """
    Create random DNA sequence string of the specified length.
    
    Parameters
    ----------
    length : int
        The length of the DNA string being generated.
        
    Returns
    -------
    str
        Random DNA string
    """
    char_list = [random.choice("GATC") for _ in range(length)]
    return "".join(char_list)
    

def write_fasta(seq_string, directory, file_name, ident="", name="", description=""):
    """
    Write a sequence to a fasta file.
    
    Parameters
    ----------
    seq_string : str
        Sequence of nucleotide characters to write to the fasta file
    directory : str
        Directory of the fasta file to write
    file_name : str
        File name of the fasta file to write
    ident : str, optional
        Identifier such as a locus tag (string)
    name : str, optional
        Sequence name
    description : str, optional
        Additional text
        
    Returns
    -------
    str
        Path of fasta file
    """
    seq = Seq(seq_string)
    record = SeqRecord(seq, ident, name, description)
    file_path = os.path.join(directory, file_name)
    SeqIO.write([record], file_path, "fasta")
    return file_path

    
def write_random_dna_fasta(directory, file_name, length):
    """
    Write a random dna string of a specified length to a fasta file.
    The fasta file will also contain a random ID, name, and description.

    Parameters
    ----------
    directory : str
        Directory of the fasta file to write
    file_name : str
        File name of the fasta file to write
    length : int
        Length of nucleotide sequence generate and write to the fasta file
        
    Returns
    -------
    file_path : str
        Path of fasta file
    dna : str
        Generated DNA string
    """
    dna = make_random_dna_string(length)
    random_num_str = str(random.choice(range(10, 100))) # two digit number
    ident = "Id" + random_num_str
    name = "Name" + random_num_str
    description = "Description" + random_num_str
    file_path = write_fasta(dna, directory, file_name, ident, name, description)
    return (file_path, dna)


def compare_fasta_files(file_path1, file_path2):
    """
    Determine if two fasta files have equivalent contents.  Not necessarily 
    byte-by-byte identical.

    Parameters
    ----------
    file_path1 : str
        Path to the first file.
    file_path2 : str
        Path to the second file.
        
    Returns
    -------
    bool
        True if the two fasta files are equivalent, False otherwise.
    """
    handle1 = open(file_path1)
    handle2 = open(file_path2)
    iter1 = SeqIO.parse(handle1, 'fasta')
    iter2 = SeqIO.parse(handle2, 'fasta')
    fillvalue = object()
    paired_records = izip_longest(iter1, iter2, fillvalue=fillvalue)
    return all(r1.__dict__ == r2.__dict__ for r1, r2 in paired_records)
    

class TestMutator(unittest.TestCase):

    def setUp(self):
        pass

    def test_zero_changes(self):
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 50)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.num_subs = 0
        args.num_insertions = 0
        args.num_deletions = 0
        args.random_seed = 1
        args.summary_file = None
        mutator.main(args)
        no_change = compare_fasta_files(original_file_path, "original_mutated_0.fasta")
        assert(no_change)
    

    def tearDown(self):
        #TempDirectory.cleanup_all()
        pass

if __name__ == '__main__':
    unittest.main()