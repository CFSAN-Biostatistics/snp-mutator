#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import os.path
from numpy import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def read_fasta_sequence(fasta_file_path):
    """
    Read a fasta file and return the sequence as a string.

    Parameters
    ----------
    fasta_file_path : str
        Path to the fasta file to read.

    Returns
    -------
    seq_name : str
        The fasta description of the first sequence in the file only.
    sequence : str
        Sequence string.
        Note: if the fast file contains multiple sequences, all the sequences
        are combined into a single string.
    """
    seq_string = []
    with open(fasta_file_path, "r") as infile:
        for idx, line in enumerate(infile):
            if idx == 0:
                seq_name = line.strip(">").strip("\n")
            elif ">" in line:
                pass
            else:
                seq_string.append(line.strip("\n"))

    seq_str = "".join(seq_string)
    return (seq_name, seq_str)


def write_fasta_sequence(seq_id, file_path, sequence_list, mutations):
    """
    Write a mutated sequence to a fasta file.  The fasta defline is suffixed
    with a mutation description of the form (mutated s=100 i=20 d=20)
    describing the number of substitutions, insertions, and deletions.

    Parameters
    ----------
    seq_id : str
        ID of the sequence which will be written to the fasta description line.
    file_path : str
        Path to the fasta file to write.  An existing file will be overwritten.
    sequence_list : list of str
        List of sequence fragments which may be 0 or more bases.  The fragments
        are concatenated before writing to the fast a file.
    mutations : tuple
        3-tuple of number of substitutions, insertions, and deletions in that
        order.  The mutation counts are appended to the fasta defline.
    """
    seq_str = "".join(sequence_list)
    seq = Seq(seq_str)
    description = "(mutated s=%i i=%i d=%i)" % mutations
    record = SeqRecord(seq, id=seq_id, description=description)
    SeqIO.write([record], file_path, "fasta")


def build_mutated_seq(seq_str, num_subs, num_insertions, num_deletions):
    """
    Copy a sequence and randomly introduce the specified numbers of
    substitutions, insertions, and deletions.

    Parameters
    ----------
    seq_str : str
        Sequence string.
    num_subs : int
        Number of base substitutions to make.
    num_insertions : int
        Number of base insertions to make.  Insertions are placed before the
        original base at positions having insertions.
    num_deletions : int
        Number of base deletions to make.

    Returns
    -------
    new_indexed_seq : list of str
        List indexed by original position containing strings of bases.  Each
        string can be 0 - 2 bases long, where a zero length string indicates
        a deletion at the position, a string containing a single base indicates
        no mutation at the position, and a string containing two bases
        indicates an insertion prior to the original base at the position.
    subs_positions : array of int
        Array of positions where substitions are introduced.
    insertion_positions : array of int
        Array of positions where insertions are introduced.
    deletion_positions : array of int
        Array of positions where bases are deleted.
    """
    substitution_choices = {"A" : ["C", "T", "G"],
                            "C" : ["A", "T", "G"],
                            "T" : ["C", "A", "G"],
                            "G" : ["C", "T", "A"],
                           }

    seq_length = len(seq_str)
    num_mutations = num_subs + num_deletions + num_insertions
    positions = random.choice(seq_length, size=num_mutations, replace=False)
    subs_positions = positions[: num_subs]
    deletion_positions = positions[num_subs : (num_subs + num_deletions)]
    insertion_positions = positions[(num_subs + num_deletions) : len(positions)]

    # Copy the original sequence in a way that easily allows mutations
    # while preserving position information
    new_indexed_seq = list(seq_str)

    for i in subs_positions:
        original_base = seq_str[i]
        new_base = random.choice(substitution_choices[original_base], size=1)[0]
        new_indexed_seq[i] = new_base

    for i in deletion_positions:
        new_indexed_seq[i] = ""

    for i in insertion_positions:
        original_base = seq_str[i]
        insert_base = random.choice(["A", "C", "T", "G"], size=1)[0]
        new_indexed_seq[i] = insert_base + original_base

    return (new_indexed_seq, subs_positions, insertion_positions, deletion_positions, )


def run_simulations(seq_str, base_file_name, seq_name, num_sims, num_subs, num_insertions, num_deletions, summary_file_path=None):
    """
    Generate multiple random mutations of a reference sequence, repeatedly 
    calling build_mutated_seq() to create each of the mutated sequences.

    Parameters
    ----------
    seq_str : str
        Original sequence string.
    base_file_name : str
        The base file name of the original reference with the extension
        removed.  This will be the file name prefix of the generated mutated
        files.
    seq_name : str
        ID of the sequence which will be written to the fasta description line.
    num_sims : int
        Number of mutated sequenced to generate.
    num_subs : int
        Number of base substitutions to make.
    num_insertions : int
        Number of base insertions to make.  Insertions are placed before the
        original base at positions having insertions.
    num_deletions : int
        Number of base deletions to make.
    summary_file_path : str, optional
        Path to summary file where a list of positions and corresponding
        mutations will be written.
    """
    if summary_file_path:
        snp_list_file = open(summary_file_path, "w")

    try:
        if summary_file_path:
            snp_list_file.write("Replicate\tPosition\tOriginalBase\tNewBase\n")

        for replicate in range(1, num_sims + 1):
            print("Creating replicate %i" % replicate)
            replicate_name = base_file_name + "_mutated_" + str(replicate)

            new_indexed_seq, subs_positions, insertion_positions, deletion_positions = \
                build_mutated_seq(seq_str, num_subs, num_insertions, num_deletions)

            mutations = (num_subs, num_insertions, num_deletions)
            write_fasta_sequence(seq_name, replicate_name + ".fasta", new_indexed_seq, mutations)

            if summary_file_path:
                summary_list = list()
                summary_list.extend([(pos, new_indexed_seq[pos]) for pos in subs_positions])
                summary_list.extend([(pos, new_indexed_seq[pos] + "_insertion") for pos in insertion_positions])
                summary_list.extend([(pos, "_deletion") for pos in deletion_positions])
                for pos, change in sorted(summary_list):
                    snp_list_file.write("%s\t%i\t%s\t%s\n" % (replicate_name, pos + 1, seq_str[pos], change))

    finally:
        if summary_file_path:
            snp_list_file.close()


def parse_arguments(system_args):
    """
    Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    usage = """Generate mutated sequence files from a reference genome.  Takes a fasta file and creates
               a specified number of randomly generated base substitutions, insertions, and deletions.
               Outputs the mutated genomes, and optionally, a summary file listing the mutations by
               position."""
               
    def non_negative_int(value):
        try:
            ivalue = int(value)
        except:
            raise argparse.ArgumentTypeError("Must be a number >= 0")
        if ivalue < 0:
            raise argparse.ArgumentTypeError("Must be >= 0")
        return ivalue
               
    def positive_int(value):
        try:
            ivalue = int(value)
        except:
            raise argparse.ArgumentTypeError("Must be a number greater than 0")
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("Must be greater than 0")
        return ivalue

    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                                             dest="input_fasta_file", type=str,                            help="Input fasta file.")
    parser.add_argument("-o", "--summary",           metavar="FILE", dest="summary_file",     type=str,              default=None, help="Output positional summary file.")
    parser.add_argument("-n", "--num-simulations",   metavar="INT",  dest="num_sims",         type=positive_int,     default=100,  help="Number of mutated sequences to generate.")
    parser.add_argument("-s", "--num-substitutions", metavar="INT",  dest="num_subs",         type=non_negative_int, default=500,  help="Number of substitutions.")
    parser.add_argument("-i", "--num-insertions",    metavar="INT",  dest="num_insertions",   type=non_negative_int, default=20,   help="Number of insertions.")
    parser.add_argument("-d", "--num-deletions",     metavar="INT",  dest="num_deletions",    type=non_negative_int, default=20,   help="Number of deletions.")
    parser.add_argument("-r", "--random-seed",       metavar="INT",  dest="random_seed",      type=int,              default=None, help="Random number seed; if not set, the results are not reproducible.")

    args = parser.parse_args(system_args)
    return args


def main(args):
    """
    Generate multiple random mutations of a reference sequence.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.

    See Also
    --------
    parse_arguments()
    """
    # Input file arg
    in_file_path = args.input_fasta_file
    in_file_name = os.path.basename(in_file_path)
    base_file_name, in_file_ext = os.path.splitext(in_file_name)

    # Random seed option
    random.seed(args.random_seed)

    # Read the reference and generate mutations
    seq_name, seq_str = read_fasta_sequence(in_file_path)
    seq_length = len(seq_str)

    num_mutations = args.num_subs + args.num_insertions + args.num_deletions
    if num_mutations > seq_length:
        print("ERROR: You have specified a number of substitutions that is greater than the length of the sequence", file=sys.stderr)
        sys.exit()

    run_simulations(seq_str, base_file_name, seq_name, args.num_sims, args.num_subs,
                    args.num_insertions, args.num_deletions, args.summary_file)


if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)


