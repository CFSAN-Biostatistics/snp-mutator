#!/usr/bin/env python
 
import argparse
import sys
import os.path
from numpy.random import choice,seed


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
                
### Deal with random seed and generate the positions where the numSubs will occur   
            try:
                seed(ranSeed)
            except NameError:
                pass
                
            try:
                positions = choice(range(0, seqLength[0] - 1), size=int(numSubs + numDels + numInsertions), replace=False)
            except ValueError:
                print "ERROR: You have specified a number of substitutions that is greater than the length of the sequence"
                sys.exit()
            subPositions = positions[:numSubs]
            deletionPositions = positions[numSubs:(numSubs + numDels)]
            insertionPositions = positions[(numSubs + numDels):len(positions)]

### Create a copy of the sequence dictionary and mutate sites that were identified in the previous step
            print "Creating replicate ", str(replicate)

            def buildNewSeq():
                newSeqDict = dict(seqDict)  
                for i in subPositions:
                    try:
                        seqDict[i]
                        state = seqDict[i]
                        if state == "A":
                            newState = str(choice(["C", "T", "G"], size = 1)[0])
                            newSeqDict[i] = newState
                            snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                        if state == "C":
                            newState = str(choice(["A", "T", "G"], size = 1)[0])
                            newSeqDict[i] = newState
                            snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                        if state == "T":
                            newState = str(choice(["C", "A", "G"], size = 1)[0])
                            newSeqDict[i] = newState
                            snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                        if state == "G":
                            newState = str(choice(["C", "T", "A"], size = 1)[0])
                            newSeqDict[i] = newState
                            snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "\n") 
                    except:
                        pass

                for i in deletionPositions:
                    try:
                        seqDict[i]
                        state = seqDict[i]
                        newState = ""
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_deletion\n") 
                    except:
                        pass

                for i in insertionPositions:
                    try:
                        seqDict[i]
                        state = seqDict[i]
                        newState = str(choice(["A", "C", "T", "G"], size = 1)[0]) + state
                        newSeqDict[i] = newState
                        snpList.write(str(replicate) + "\t" + str(i + 1) + "\t" + state + "\t" + newState + "_insertion\n") 
                    except:
                        pass
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
    # Input file arg
    in_file = args.input_fasta_file
    in_filePath, in_fileWholeName = os.path.split(in_file)
    in_fileBase, in_fileExt = os.path.splitext(in_fileWholeName)
    
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
    
    # Read the reference and generate mutations
    getSequence(in_file)
    buildSeqDict(seqString)
    runSimulations(in_fileBase, numSims, numSubs, numDels, numInsertions, seqDict, ranSeed)


if __name__ == '__main__':
    args = parse_arguments(sys.argv)
    main(args)
    

