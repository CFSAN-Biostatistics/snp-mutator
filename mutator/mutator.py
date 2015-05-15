#!/usr/bin/env python
 
from __future__ import print_function
import argparse
import sys
import os.path
from numpy import random
from timeit import Timer


TIMING_RUNS = 10
ENABLE_TIMING_TEST = False

### Read in reference and build dictionary - I originally used SeqIO.parse but went this route to account for a references with more than one contig

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
    with open(fasta_file_path, "r") as file:
        for idx, line in enumerate(file):
            if idx == 0:
                seq_name = line.strip(">").strip("\n")
            elif ">" in line:
                pass
            else:
                seq_string.append(line.strip("\n"))

    seq_str = "".join(seq_string)
    return (seq_name, seq_str)


def write_fasta_sequence(seq_name, file_path, sequence_list):
    with open(file_path, "w") as file:
        file.write(">" + seq_name + "\n")
        file.write("".join(sequence_list))


### Run simulations to get mutated genome
def runSimulations(in_fileBase, seq_name, numSims, numSubs, numDels, numInsertions, seq_str, ranSeed):
    seq_length = len(seq_str)
    with open(in_fileBase + "_snpListMutated.txt", "w") as snpList:
        snpList.write("replicate\tposition\toriginalBase\tnewBase\n")

        for replicate in range(0, numSims):  
            positions = random.choice(range(0, seq_length), size=int(numSubs + numDels + numInsertions), replace=False)
            subPositions = positions[:numSubs]
            deletionPositions = positions[numSubs:(numSubs + numDels)]
            insertionPositions = positions[(numSubs + numDels):len(positions)]

            print("Creating replicate %i" % replicate)

            def buildNewSeq():
                """
                Create a copy of the sequence and mutate sites that were 
                identified in the previous step
                """
                substitution_choices = {"A" : ["C", "T", "G"],
                                        "C" : ["A", "T", "G"],
                                        "T" : ["C", "A", "G"],
                                        "G" : ["C", "T", "A"],
                                        }
                new_indexed_seq = list(seq_str)
                for i in subPositions:
                    state = seq_str[i]
                    newState = random.choice(substitution_choices[state], size = 1)[0]
                    new_indexed_seq[i] = newState
                    snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 

                for i in deletionPositions:
                    state = seq_str[i]
                    newState = ""
                    new_indexed_seq[i] = newState
                    snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_deletion\n") 

                for i in insertionPositions:
                    state = seq_str[i]
                    newState = str(random.choice(["A", "C", "T", "G"], size = 1)[0]) + state
                    new_indexed_seq[i] = newState
                    snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_insertion\n") 
                    
                return new_indexed_seq

            if ENABLE_TIMING_TEST:
                t = Timer(lambda: buildNewSeq())
                print("Min buildNewSeq Time = %f" % ( min(t.repeat(repeat=TIMING_RUNS, number=1))))

            new_indexed_seq = buildNewSeq()
            
            
            if ENABLE_TIMING_TEST:
                t = Timer(lambda: write_fasta_sequence(seq_name, in_fileBase + "_mutated_" + str(replicate) + ".fasta", new_indexed_seq))
                print("Min Write Time = %f" % ( min(t.repeat(repeat=TIMING_RUNS, number=1))))
            
            write_fasta_sequence(seq_name, in_fileBase + "_mutated_" + str(replicate) + ".fasta", new_indexed_seq)
            
            
def parse_arguments(system_args):
    """
    """
    usage = """Generate mutated sequence files from a reference genome.  Takes a fasta file and creates 
               a specified number of randomly generated base substitutions, insertions, and deletions.  
               Outputs the mutated genomes, and optionally, a summary file listing the mutations by 
               position."""

    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                                             dest="input_fasta_file", type=str,               help="Input fasta file.")
    parser.add_argument("-o", "--summary",           metavar="FILE", dest="summary_file",     type=str, default=None, help="Output positional summary file.")
    parser.add_argument("-n", "--num-simulations",   metavar="INT",  dest="num_sims",         type=int, default=100,  help="Number of mutated sequences to generate.")
    parser.add_argument("-s", "--num-substitutions", metavar="INT",  dest="num_subs",         type=int, default=500,  help="Number of substitutions.")
    parser.add_argument("-i", "--num-insertions",    metavar="INT",  dest="num_insertions",   type=int, default=20,   help="Number of insertions.")
    parser.add_argument("-d", "--num-deletions",     metavar="INT",  dest="num_deletions",    type=int, default=20,   help="Number of deletions.")
    parser.add_argument("-r", "--random-seed",       metavar="INT",  dest="random_seed",      type=int, default=None, help="Random number seed making the results reproducible.")
    
    args = parser.parse_args()
    return args
            

def main(args):
    """
    """
    # Input file arg
    in_file = args.input_fasta_file
    in_file_name = os.path.basename(in_file)
    in_fileBase, in_fileExt = os.path.splitext(in_file_name)
    
    # Summary arg
    if args.summary_file:
        summary_file = os.path.abspath(args.summary_file)
    
    # Get number of simulations, substitutions and indels
    numSims = args.num_sims
    numSubs = args.num_subs
    numInsertions = args.num_insertions
    numDels = args.num_deletions
    
    # Random seed option
    ranSeed = args.random_seed
    random.seed(ranSeed)
    
    # Read the reference and generate mutations
    if ENABLE_TIMING_TEST:
        t = Timer(lambda: read_fasta_sequence(in_file))
        print("Min Read Time = %f" % ( min(t.repeat(repeat=TIMING_RUNS, number=1))))
    seq_name, indexed_sequence = read_fasta_sequence(in_file)
    seq_length = len(indexed_sequence)
    
    num_mutations = args.num_subs + args.num_insertions + args.num_deletions
    if num_mutations > seq_length:
        print("ERROR: You have specified a number of substitutions that is greater than the length of the sequence", file=sys.stderr)
        sys.exit()
    
    runSimulations(in_fileBase, seq_name, numSims, numSubs, numDels, numInsertions, indexed_sequence, ranSeed)


if __name__ == '__main__':
    args = parse_arguments(sys.argv)
    main(args)
    

