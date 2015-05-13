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


def read_fasta_seq_record(file_path):
    """
    Reads a fasta file and returns the first sequence record. Any subsequent sequences are ignored.

    Parameters
    ----------
    file_path : str
        Path to the fasta file
        
    Returns
    -------
    SeqRecord
        Biopython SeqRecord
    """
    with open(file_path) as handle:
        iter = SeqIO.parse(handle, 'fasta')
        return iter.next()


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
    with open(file_path1) as handle1:
        with open(file_path2) as handle2:
            iter1 = SeqIO.parse(handle1, 'fasta')
            iter2 = SeqIO.parse(handle2, 'fasta')
            fillvalue = object()
            paired_records = izip_longest(iter1, iter2, fillvalue=fillvalue)
            return all(r1.__dict__ == r2.__dict__ for r1, r2 in paired_records)


class TestMutator(unittest.TestCase):

    def setUp(self):
        random.seed(1)

    def test_zero_changes(self):
        """Verify the output fasta file matches the input fasta file when zero substitions, insertions, and deletions are requested.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
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
        self.assertTrue(no_change, "Generated fasta file does not match original fasta file")

    def test_snp_changes(self):
        """Test various numbers of substitutions.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.num_insertions = 0
        args.num_deletions = 0
        args.random_seed = 1
        args.summary_file = None

        args.num_subs = 1
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCTAAATCGG", "SNP 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 5
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "ACTATAGCGC", "SNP 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 10
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "AGTCTCGTAC", "SNP 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_insert_changes(self):
        """Test various numbers of insertions.
        """
        
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.num_subs = 0
        args.num_deletions = 0
        args.random_seed = 1
        args.summary_file = None

        args.num_insertions = 1
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCGCAAATCGG", "Insert 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_insertions = 5
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "CGCGCATAAATCGCG", "Insert 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_insertions = 10
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "CGACGCTATATAATTCCGCG", "Insert 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_delete_changes(self):
        """Test various numbers of deletions.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.num_subs = 0
        args.num_insertions = 0
        args.random_seed = 1
        args.summary_file = None

        args.num_deletions = 1
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCAAATCGG", "Delete 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_deletions = 5
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "CAACG", "Delete 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_deletions = 10
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "", "Delete 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_mutate_mix_changes(self):
        """Test a mix of substitutions, inserts, and deletes.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.random_seed = 1
        args.summary_file = None

        args.num_subs = 1
        args.num_insertions = 1
        args.num_deletions = 1
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCTAAAATCG", "Mutate mix 1,1,1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "TGCTCAACGC", "Mutate mix 2,2,2 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 3
        args.num_insertions = 4
        args.num_deletions = 3
        mutator.main(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_0.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "CCTTAGTCAGC", "Mutate mix 3,4,3 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))
            
    def test_too_many_mutations(self):
        """Deliberately exceed the maximum allowed number of mutations.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 1
        args.random_seed = 1
        args.summary_file = None

        args.num_subs = 4
        args.num_insertions = 4
        args.num_deletions = 3
        with self.assertRaises(SystemExit):  # Verify exit if number of mutations exceeds sequence length
            mutator.main(args)
            
    def test_not_all_same(self):
        """Verify Mutator creates different mutated fasta files when generating more than one.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = argparse.Namespace()
        args.input_fasta_file = original_file_path
        args.num_sims = 3
        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        args.random_seed = 1
        args.summary_file = None
        mutator.main(args)
        mutated_seq_record0 = read_fasta_seq_record("original_mutated_0.fasta")
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertNotEqual(str(mutated_seq_record0.seq), str(mutated_seq_record1.seq), "Generated sequences 0 and 1 should be different.")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Generated sequences 1 and 2 should be different.")
        self.assertNotEqual(str(mutated_seq_record0.seq), str(mutated_seq_record2.seq), "Generated sequences 0 and 2 should be different.")

    def tearDown(self):
        """
        Delete all the temporary directories and files created during this 
        testing session.
        """
        if os.path.exists("original_mutated_0.fasta"):
            os.remove("original_mutated_0.fasta")
        if os.path.exists("original_mutated_1.fasta"):
            os.remove("original_mutated_1.fasta")
        if os.path.exists("original_mutated_2.fasta"):
            os.remove("original_mutated_2.fasta")
        if os.path.exists("original_snpListMutated.txt"):
            os.remove("original_snpListMutated.txt")
        TempDirectory.cleanup_all()

if __name__ == '__main__':
    unittest.main()
