#!/usr/bin/env python
 
import argparse
import sys
import os.path
from numpy.random import choice,seed
import subprocess
import shutil


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
def runSimulations(seqDict):
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
			


if __name__ == '__main__':
	usage="""Takes a fasta file and creates N base substitutions and N insertion/deletions. Assumes you have the simulation program installed in your path and outputs a summary file, mutated genomes, and fastq files form art.
	usage: %prog -i FILE (args)"""


	parser = argparse.ArgumentParser(description= usage)


	parser.add_argument("-i", "--in-file", metavar="FILE", dest="in_file", default=None,
					help="Specify the input fasta FILE")

	parser.add_argument("-s", "--summary", metavar="FILE", dest="summary_file", default=None,
					help="Specify the name of the file to output positional information to.")

	parser.add_argument("-n", "--number-of-simulations", metavar="INT", dest="num_sims", default=None,
					help="Specify the number of simulations to run. Default is 100.")

	parser.add_argument("-m", "--number-of-substitutions", metavar="INT", dest="num_subs", default=None,
					help="Specify the number of substitutions. Default is 500.")

	parser.add_argument("-t", "--num-of-insertions", metavar="INT", dest="num_insertions", default=None,
					help="Specify the number of insertions. Default is 20.")

	parser.add_argument("-d", "--num-of-deletions", metavar="INT", dest="num_deletions", default=None,
					help="Specify the number of deletions. Default is 20.")

	parser.add_argument("-e", "--random-number-seed", metavar="INT", dest="ran_seed", default=None,
					help="Specify the random number seed - this will cause all mutated genomes to be identical. Default is to get one from the computer.")

	args = parser.parse_args()


	### Input file args
	if args.in_file:
		in_file = os.path.abspath(args.in_file)
	else:	
		print "You must specify an in_file. Use '-h' for help."
		sys.exit()

	(in_filePath, in_fileWholeName) = os.path.split(in_file)
	(in_fileBase, in_fileExt) = os.path.splitext(in_fileWholeName)

	###	 Summary args
	if args.summary_file:
		contigs_file = os.path.abspath(args.summary_file)
	else:
		summary_file = os.path.join(in_filePath, in_fileBase + "_mutatedSNPlist.txt")

	### Define number of simulations, substitutions and indels
	if args.num_sims == None:
		numSims = 100
	else:
		numSims = int(args.num_sims)

	if args.num_subs == None:
		numSubs = 500
	else:
		numSubs = int(args.num_subs)

	if args.num_deletions == None:
		numDels = 20
	else:
		numDels = int(args.num_deletions)

	if args.num_insertions == None:
		numInsertions = 20
	else:
		numInsertions = int(args.num_insertions)

	### Parse random seed options
	if args.ran_seed == None:
		pass
	else:
		ranSeed = int(args.ran_seed)

	getSequence(in_file)	
	buildSeqDict(seqString)	
	runSimulations(seqDict)
