#!/usr/bin/env python
 
from __future__ import print_function
import argparse
import sys
import os.path
from numpy import random


### Read in reference and build dictionary - I originally used SeqIO.parse but went this route to account for a references with more than one contig
seqString = []
seqDict = {}
seqName = []
seqLength = []

def getSequence(in_file):
    with open(in_file, "r") as inFile:
        for idx, line in enumerate(inFile):
            if idx == 0:
                seqName.append(line.strip(">").strip("\n"))
            elif ">" in line:
                pass
            else:
                seqString.append(line.strip("\n"))


def buildSeqDict(seqString):
    seq = "".join(seqString)
    seqLength.append(len(seq))
    for idx, i in enumerate(seq):
        seqDict[idx]= i



### Run simulations to get mutated genome
def runSimulations(in_fileBase, numSims, numSubs, numDels, numInsertions, seqDict, ranSeed):
    with open(in_fileBase + "_snpListMutated.txt", "w") as snpList:
        snpList.write("replicate\tposition\toriginalBase\tnewBase\n")

        for replicate in range(0, numSims):  
            positions = random.choice(range(0, seqLength[0]), size=int(numSubs + numDels + numInsertions), replace=False)
            subPositions = positions[:numSubs]
            deletionPositions = positions[numSubs:(numSubs + numDels)]
            insertionPositions = positions[(numSubs + numDels):len(positions)]

### Create a copy of the sequence dictionary and mutate sites that were identified in the previous step
            print("Creating replicate %i" % replicate)

            def buildNewSeq():
                newSeqDict = dict(seqDict)  
                for i in subPositions:
                    state = seqDict[i]
                    if state == "A":
                        newState = str(random.choice(["C", "T", "G"], size = 1)[0])
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                    if state == "C":
                        newState = str(random.choice(["A", "T", "G"], size = 1)[0])
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                    if state == "T":
                        newState = str(random.choice(["C", "A", "G"], size = 1)[0])
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                    if state == "G":
                        newState = str(random.choice(["C", "T", "A"], size = 1)[0])
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 

                for i in deletionPositions:
                    state = seqDict[i]
                    newState = ""
                    newSeqDict[i] = newState
                    snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_deletion\n") 

                for i in insertionPositions:
                    state = seqDict[i]
                    newState = str(random.choice(["A", "C", "T", "G"], size = 1)[0]) + state
                    newSeqDict[i] = newState
                    snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_insertion\n") 

                def writeNewSequence():
                    newSequence = []
                    with open(in_fileBase + "_mutated_" + str(replicate) + ".fasta", "w") as newFasta:
                        newFasta.write(">" + seqName[0] + "\n")
                        for key, value in newSeqDict.items():
                            newSequence.append(value)   
                        newFasta.write("".join(newSequence))
                writeNewSequence()
            buildNewSeq()
            
            
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
    # TODO: remove these globals
    global seqLength
    global seqDict
    global seqName
    global seqString
    seqLength = []
    seqDict = {}
    seqName = []
    seqString = []
    
    
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
    getSequence(in_file)
    buildSeqDict(seqString)
    
    num_mutations = args.num_subs + args.num_insertions + args.num_deletions
    if num_mutations > seqLength[0]:
        print("ERROR: You have specified a number of substitutions that is greater than the length of the sequence", file=sys.stderr)
        sys.exit()
    
    runSimulations(in_fileBase, numSims, numSubs, numDels, numInsertions, seqDict, ranSeed)


if __name__ == '__main__':
    args = parse_arguments(sys.argv)
    main(args)
    

