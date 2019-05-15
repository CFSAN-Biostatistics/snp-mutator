#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_snpmutator
----------------------------------

Tests for `snpmutator` module.
"""

from __future__ import print_function

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import unittest
from testfixtures import TempDirectory
import random
import argparse
try:
    # Python 3
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest

from snpmutator import script as snpmutator


def make_default_args(original_fasta_file_path):
    """Create a default set of arguments that does nothing.
    """
    args = argparse.Namespace()
    args.input_fasta_file = original_fasta_file_path
    args.num_sims = 0
    args.num_subs = 0
    args.num_insertions = 0
    args.num_deletions = 0
    args.random_seed = 0
    args.subset_len = 0
    args.group_size = None
    args.mono = False
    args.concat_ref_file = None
    args.summary_file = None
    args.seq_id = None
    args.vcf_file = None
    args.metrics_file = None
    args.fasta_output_dir = "."

    return args


def make_random_dna_string(length, allowed_bases):
    """
    Create random DNA sequence string of the specified length.
    Parameters
    ----------
    length : int
        The length of the DNA string being generated.
    allowed_bases : str
        A string of possible bases.
    Returns
    -------
    str
        Random DNA string
    """
    char_list = [random.choice(allowed_bases) for _ in range(length)]
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


def write_fixed_dna_fasta(dna, directory, file_name):
    """
    Write a dna string to a fasta file.
    The fasta file will contain a random ID, name, and description.

    Parameters
    ----------
    dna : str
        String of DNA to write to fasta file.
    directory : str
        Directory of the fasta file to write
    file_name : str
        File name of the fasta file to write

    Returns
    -------
    file_path : str
        Path of fasta file
    """
    random_num_str = str(random.choice(range(10, 100))) # two digit number
    ident = "Id" + random_num_str
    name = "Name" + random_num_str
    description = "Description" + random_num_str
    file_path = write_fasta(dna, directory, file_name, ident, name, description)
    return file_path


def write_random_dna_fasta(directory, file_name, length, allowed_bases="GATC"):
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
    allowed_bases : str, optional
        A string of possible bases.  Defaults to 'GATC'

    Returns
    -------
    file_path : str
        Path of fasta file
    dna : str
        Generated DNA string
    """
    dna = make_random_dna_string(length, allowed_bases)
    file_path = write_fixed_dna_fasta(dna, directory, file_name)
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
        return next(iter)


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
            paired_records = zip_longest(iter1, iter2, fillvalue=fillvalue)
            return all(r1.__dict__ == r2.__dict__ for r1, r2 in paired_records)


def compare_mutated_fasta_files(original_file_path, mutated_file_path):
    """
    Determine if two fasta files have equivalent contents after running a
    zero count mutation on the original file.  The files will not be identical.
    The fasta defline will differ in the description suffix and the sequence
    lines may have different lengths.

    Parameters
    ----------
    original_file_path : str
        Path to the original fasta file.
    mutated_file_path : str
        Path to the mutated fasta file.

    Returns
    -------
    bool
        True if the two fasta files are equivalent, False otherwise.
    """
    with open(original_file_path) as handle1:
        with open(mutated_file_path) as handle2:
            iter1 = SeqIO.parse(handle1, 'fasta')
            iter2 = SeqIO.parse(handle2, 'fasta')
            fillvalue = object()
            paired_records = zip_longest(iter1, iter2, fillvalue=fillvalue)
            for r1, r2 in paired_records:
                if type(r1) != type(r2):
                    return False
                if r1.id != r2.id:
                    return False
                if r1.name != r2.name:
                    return False
                if r1.seq != r2.seq:
                    return False
                if not r2.description.startswith(r1.description):
                    return False
    return True


class TestSnpmutator(unittest.TestCase):

    def setUp(self):
        random.seed(1)

    def test_invalid_args(self):
        """Verify the argument parser rejects invalid arguments
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 100)

        sys_arguments = ("-s -1 " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-i -1 " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-d -1 " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-s aa " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-i aa " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-d aa " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-n -1 " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

        sys_arguments = ("-n aa " + original_file_path).split()
        self.assertRaises(SystemExit, snpmutator.parse_arguments, sys_arguments)

    def test_zero_changes(self):
        """Verify the output fasta file matches the input fasta file when zero substitions, insertions, and deletions are requested.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1
        snpmutator.run_from_args(args)
        no_change = compare_mutated_fasta_files(original_file_path, "original_mutated_1.fasta")
        self.assertTrue(no_change, "Generated fasta file does not match original fasta file")

    def test_snp_changes(self):
        """Test various numbers of substitutions.
        """
        directory = TempDirectory()
        dna = "GCCAAATCGG"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_subs = 1
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCTAAATCGG", "SNP 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 5
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "ACTATAGCGC", "SNP 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 10
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "AGTCTCGTAC", "SNP 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_eligible_snp_changes(self):
        """Test substitutions where some positions are ineligible
        """
        directory = TempDirectory()
        dna = "12345aaaaaAAAAA12345"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_subs = 10
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "12345GGTCTCGTGC12345", "Eligible SNP test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_insert_changes(self):
        """Test various numbers of insertions.
        """
        directory = TempDirectory()
        dna = "TTTTAATTTT"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_insertions = 1
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "TTTGTAATTTT", "Insert 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_insertions = 5
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "TCTTGTATATATTTC", "Insert 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_insertions = 10
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "TCTATGTTATATTATTTCTC", "Insert 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_delete_changes(self):
        """Test various numbers of deletions.
        """
        directory = TempDirectory()
        dna = "GCCAAATCGG"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_deletions = 1
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCAAATCGG", "Delete 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_deletions = 5
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "CAACG", "Delete 5 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_deletions = 10
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "", "Delete 10 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_mutate_mix_changes(self):
        """Test a mix of substitutions, inserts, and deletes.
        """
        directory = TempDirectory()
        dna = "GGGGGGGGGG"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_subs = 1
        args.num_insertions = 1
        args.num_deletions = 1
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GGTGGGGAGG", "Mutate mix 1,1,1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GTGTGCGGGC", "Mutate mix 2,2,2 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

        args.num_subs = 3
        args.num_insertions = 4
        args.num_deletions = 3
        snpmutator.run_from_args(args)
        mutated_seq_record = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(str(mutated_seq_record.seq), "GCTGTAGTGAC", "Mutate mix 3,4,3 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record.seq)))

    def test_too_many_mutations(self):
        """Deliberately exceed the maximum allowed number of mutations.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_subs = 4
        args.num_insertions = 4
        args.num_deletions = 3
        # Verify exit if number of mutations exceeds sequence length
        self.assertRaises(SystemExit, snpmutator.run_from_args, args)

    def test_too_many_mutations_pool(self):
        """Deliberately exceed the maximum allowed number of mutations.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1
        args.subset_len = 9

        args.num_subs = 4
        args.num_insertions = 3
        args.num_deletions = 3
        # Verify exit if number of mutations exceeds sequence length
        self.assertRaises(SystemExit, snpmutator.run_from_args, args)

    def test_too_many_mutations_ineligible(self):
        """Deliberately exceed the maximum allowed number of mutations when
        some of the reference positions are not ATGC.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 10, "ACGTSNWR")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        args.num_subs = 10
        args.num_insertions = 0
        args.num_deletions = 0
        # Verify exit if number of mutations exceeds eligible sequence length
        self.assertRaises(SystemExit, snpmutator.run_from_args, args)

    def test_not_all_same(self):
        """Verify Mutator creates different mutated fasta files when generating more than one.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 3
        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2

        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        mutated_seq_record3 = read_fasta_seq_record("original_mutated_3.fasta")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Generated sequences 1 and 2 should be different.")
        self.assertNotEqual(str(mutated_seq_record2.seq), str(mutated_seq_record3.seq), "Generated sequences 2 and 3 should be different.")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record3.seq), "Generated sequences 1 and 3 should be different.")

    def test_not_all_same_pool(self):
        """Verify Mutator creates different mutated fasta files when generating more than one.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 3
        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        args.subset_len = 500

        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        mutated_seq_record3 = read_fasta_seq_record("original_mutated_3.fasta")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Generated sequences 1 and 2 should be different.")
        self.assertNotEqual(str(mutated_seq_record2.seq), str(mutated_seq_record3.seq), "Generated sequences 2 and 3 should be different.")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record3.seq), "Generated sequences 1 and 3 should be different.")

    def test_not_all_same_mono(self):
        """Verify Mutator creates different mutated fasta files when generating more than one.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 3
        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        args.mono = True

        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        mutated_seq_record3 = read_fasta_seq_record("original_mutated_3.fasta")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Generated sequences 1 and 2 should be different.")
        self.assertNotEqual(str(mutated_seq_record2.seq), str(mutated_seq_record3.seq), "Generated sequences 2 and 3 should be different.")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record3.seq), "Generated sequences 1 and 3 should be different.")

    def test_not_all_same_pool_mono(self):
        """Verify Mutator creates different mutated fasta files when generating more than one.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 3
        args.num_subs = 2
        args.num_insertions = 2
        args.num_deletions = 2
        args.subset_len = 500
        args.mono = True

        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        mutated_seq_record3 = read_fasta_seq_record("original_mutated_3.fasta")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Generated sequences 1 and 2 should be different.")
        self.assertNotEqual(str(mutated_seq_record2.seq), str(mutated_seq_record3.seq), "Generated sequences 2 and 3 should be different.")
        self.assertNotEqual(str(mutated_seq_record1.seq), str(mutated_seq_record3.seq), "Generated sequences 1 and 3 should be different.")

    def test_concat_ref(self):
        """Verify the concat ref file is created if and only if requested.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        concat_ref_file = os.path.join(directory.path, "concat_ref.fasta")

        args.concat_ref_file = None
        snpmutator.run_from_args(args)
        file_exists = os.path.exists(concat_ref_file)
        self.assertFalse(file_exists, "The concat ref file should not exist when not explicitly requested")

        args.concat_ref_file = concat_ref_file
        snpmutator.run_from_args(args)
        file_exists = os.path.exists(concat_ref_file)
        self.assertTrue(file_exists, "The concat ref file is missing when requested.")


    def test_summary_creation(self):
        """Verify the summary file is created if and only if requested.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        summary_file_path = os.path.join(directory.path, "original_snpListMutated.txt")

        args.summary_file = None
        snpmutator.run_from_args(args)
        summary_file_exists = os.path.exists(summary_file_path)
        self.assertFalse(summary_file_exists, "The summary file should not exist when not explicitly requested")

        args.summary_file = summary_file_path
        snpmutator.run_from_args(args)
        summary_file_exists = os.path.exists(summary_file_path)
        self.assertTrue(summary_file_exists, "The summary file is missing when requested.")

    def test_vcf_creation(self):
        """Verify the VCF file is created if and only if requested.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 1000)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1

        vcf_file_path = os.path.join(directory.path, "test.vcf")

        args.vcf_file = None
        snpmutator.run_from_args(args)
        vcf_file_exists = os.path.exists(vcf_file_path)
        self.assertFalse(vcf_file_exists, "The VCF file should not exist when not explicitly requested")

        args.vcf_file = vcf_file_path
        snpmutator.run_from_args(args)
        vcf_file_exists = os.path.exists(vcf_file_path)
        self.assertTrue(vcf_file_exists, "The VCF file is missing when requested.")

    def test_pooling(self):
        """Verify that pooling places mutations at the same location in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 2
        args.subset_len = 1

        args.num_subs = 1
        args.num_insertions = 0
        args.num_deletions = 0
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), 'AATAAAAAAA', "Pooling SNP replicate 1 test failed, dna=%s mutated seq1=%s" % (dna, str(mutated_seq_record1.seq)))
        self.assertEqual(str(mutated_seq_record2.seq), 'AACAAAAAAA', "Pooling SNP replicate 2 test failed, dna=%s mutated seq2=%s" % (dna, str(mutated_seq_record2.seq)))

        args.num_subs = 0
        args.num_insertions = 1
        args.num_deletions = 0
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), 'AAAGAAAAAAA', "Pooling INS replicate 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record1.seq)))
        self.assertEqual(str(mutated_seq_record2.seq), 'AAACAAAAAAA', "Pooling INS replicate 2 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record2.seq)))

        args.num_subs = 0
        args.num_insertions = 0
        args.num_deletions = 1
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), 'AAAAAAAAA', "Pooling DEL replicate 1 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record1.seq)))
        self.assertEqual(str(mutated_seq_record2.seq), 'AAAAAAAAA', "Pooling DEL replicate 2 test failed, dna=%s mutated seq=%s" % (dna, str(mutated_seq_record1.seq)))

    def test_mono_snps(self):
        """Verify that Monomorphic mutations are the same in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.subset_len = 8
        args.mono = True

        args.num_sims = 2
        args.num_subs = 8
        args.num_insertions = 0
        args.num_deletions = 0
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Monomorphic SNPs do not match, mutated seq 1=%s mutated seq 2=%s" % (str(mutated_seq_record1.seq), str(mutated_seq_record2.seq)))

    def test_mono_insertions(self):
        """Verify that Monomorphic mutations are the same in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.subset_len = 8
        args.mono = True

        args.num_sims = 2
        args.num_subs = 0
        args.num_insertions = 8
        args.num_deletions = 0
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Monomorphic insertions do not match, mutated seq 1=%s mutated seq 2=%s" % (str(mutated_seq_record1.seq), str(mutated_seq_record2.seq)))

    def test_mono_deletions(self):
        """Verify that Monomorphic mutations are the same in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.subset_len = 8
        args.mono = True

        args.num_sims = 2
        args.num_subs = 0
        args.num_insertions = 0
        args.num_deletions = 8
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Monomorphic deletions do not match, mutated seq 1=%s mutated seq 2=%s" % (str(mutated_seq_record1.seq), str(mutated_seq_record2.seq)))

    def test_mono_mix(self):
        """Verify that Monomorphic mutations are the same in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.subset_len = 9
        args.mono = True

        args.num_sims = 2
        args.num_subs = 3
        args.num_insertions = 3
        args.num_deletions = 3
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Monomorphic mix of mutations do not match, mutated seq 1=%s mutated seq 2=%s" % (str(mutated_seq_record1.seq), str(mutated_seq_record2.seq)))

    def test_mono_mix_no_del(self):
        """Verify that Monomorphic mutations are the same in all replicates.
        """
        directory = TempDirectory()
        dna = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        original_file_path = write_fixed_dna_fasta(dna, directory.path, "original.fasta")
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.subset_len = 6
        args.mono = True

        args.num_sims = 2
        args.num_subs = 3
        args.num_deletions = 0
        args.num_insertions = 3
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        mutated_seq_record2 = read_fasta_seq_record("original_mutated_2.fasta")
        self.assertEqual(str(mutated_seq_record1.seq), str(mutated_seq_record2.seq), "Monomorphic mix of mutations do not match, mutated seq 1=%s mutated seq 2=%s" % (str(mutated_seq_record1.seq), str(mutated_seq_record2.seq)))

    def test_metrics(self):
        """Verify that metrics are correct.
        """
        seq_str = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
        all_eligible_positions = [i for i in range(len(seq_str))]
        base_file_name = "original"
        seq_name = ""
        num_sims = 2
        num_subs = 3
        num_deletions = 0
        num_insertions = 0
        pool_size = 0
        group_size = num_sims
        mono = True
        metrics = snpmutator.run_simulations(seq_str, all_eligible_positions, base_file_name, seq_name, num_sims, num_subs, num_insertions, num_deletions, pool_size, group_size, mono)
        self.assertEqual(metrics.replicates, 2, "Incorrect metrics.replicates")
        self.assertEqual(metrics.variant_positions, 6, "Incorrect metrics.variant_positions")
        self.assertEqual(metrics.monomorphic_variant_positions, 6, "Incorrect metrics.monomorphic_variant_positions")
        self.assertEqual(metrics.polymorphic_variant_positions, 0, "Incorrect metrics.polymorphic_variant_positions")
        self.assertEqual(metrics.single_replicate_positions, 6, "Incorrect metrics.single_replicate_positions")
        self.assertEqual(metrics.multiple_replicate_positions, 0, "Incorrect metrics.multiple_replicate_positions")

        seq_str = "AAA"
        all_eligible_positions = [i for i in range(len(seq_str))]
        num_sims = 5
        pool_size = 3
        mono = False
        metrics = snpmutator.run_simulations(seq_str, all_eligible_positions, base_file_name, seq_name, num_sims, num_subs, num_insertions, num_deletions, pool_size, group_size, mono)
        self.assertEqual(metrics.replicates, 5, "Incorrect metrics.replicates")
        self.assertEqual(metrics.variant_positions, 3, "Incorrect metrics.variant_positions")
        self.assertEqual(metrics.monomorphic_variant_positions, 0, "Incorrect metrics.monomorphic_variant_positions")
        self.assertEqual(metrics.polymorphic_variant_positions, 3, "Incorrect metrics.polymorphic_variant_positions")
        self.assertEqual(metrics.single_replicate_positions, 0, "Incorrect metrics.single_replicate_positions")
        self.assertEqual(metrics.multiple_replicate_positions, 3, "Incorrect metrics.multiple_replicate_positions")

    def test_seqid_override(self):
        """Verify output fasta files can have the defline seqid overridden.
        """
        directory = TempDirectory()
        original_file_path, dna = write_random_dna_fasta(directory.path, "original.fasta", 50)
        args = make_default_args(original_file_path)
        args.random_seed = 1
        args.num_sims = 1
        args.num_subs = 3
        original_seq_record = read_fasta_seq_record(original_file_path)

        args.seq_id = None
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(original_seq_record.id, mutated_seq_record1.id, "Defline seq id should not change when not requested.")

        args.seq_id = "test_override_seqid"
        snpmutator.run_from_args(args)
        mutated_seq_record1 = read_fasta_seq_record("original_mutated_1.fasta")
        self.assertEqual(mutated_seq_record1.id, args.seq_id, 'Overridden defline seq id "%s" does not match expected value "%s"' % (mutated_seq_record1.id, args.seq_id))


    def tearDown(self):
        """
        Delete all the temporary directories and files created during this
        testing session.
        """
        for replicate in range(1, 6):
            file_name = "original_mutated_%d.fasta" % replicate
            if os.path.exists(file_name):
                os.remove(file_name)
        if os.path.exists("original_snpListMutated.txt"):
            os.remove("original_snpListMutated.txt")
        TempDirectory.cleanup_all()


if __name__ == '__main__':
    unittest.main()
