#!usr/bin/env python3
import sys, os, re, gzip, time

########################################################################################
#### @Authors:	Josephine Strange & Jan Turowski									####
#### @Version:	1.0																	####
#### @Date:		01-Dec-2019															####
####																				####
#### 			Takes a database (FASTA-file) of resistance genes to 				####
####			scan reads for presence of antimicrobial resistance (AMR)			#### 
####			in metagenomic samples.												####
####																				####
########################################################################################
#### File prompt ####
# Define fasta extension tuple
fastaExt = (".fasta", ".fsa", ".faa", ".fna", ".ffn", ".frn")

# Prompt for a resistance gene database and an arbitrary number of read files
if len(sys.argv) == 1:
	resFile = input("Please enter path/name of the resistance gene \"database\" (fasta file): ")
	readFiles = input("Please enter path/name of the sample (gzip) files (separated by space): ").split()
	# Sanity check that all read files are actually gzipped
	for file in readFiles:
		if not file.endswith(".gz"):
			sys.stderr.write("Error! Non-gzipped file detected ({})... Terminating!\n".format(os.path.basename(file)))
			sys.exit(1)

elif len(sys.argv) > 1: 
	if sys.argv[1].endswith(fastaExt):
		resFile = sys.argv[1]
		readFiles = sys.argv[2:]
		# Sanity check that all read files are actually gzipped
		for file in readFiles:
			if not file.endswith(".gz"):
				sys.stderr.write("Error! Identified non-gzipped file ({})... Terminating!\n".format(os.path.basename(file)))
				sys.exit(1)
	else:
		sys.stderr.write("Usage: abres.py <resistance database>.fasta <read file1>.gz ... <read fileN>.gz\n") 
		sys.exit(1)

sys.stdout.write("Initializing...")
########################################################################################
#### Define functions ####

def yield_kmer(seq, kmerLen):
	"""Generates kmers from input sequence and specified kmer size"""
	for i in range(0, len(seq) - kmerLen + 1):
		yield seq[i:i+kmerLen]

def reverse_complement(seq):
	"""Reverse-complements an input kmer"""
	transTable = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
	seq = seq.translate(transTable)[::-1]
	return seq

def check_kmers(readSeq, kmerLen):
	"""Scans resistance dict for generated kmers: returns true if equal to or above certain threshold"""
	i = 0 
	kmerThreshold = 3
	obsKmerCounter = 0
	while i < len(readSeq)-kmerLen:
		kmer = readSeq[i:i+kmerLen]
		if kmer in kmerDict:
			obsKmerCounter += 1
			if obsKmerCounter >= kmerThreshold:
				return True
		i += kmerLen
	return False

########################################################################################
#### Create kmer resistance dictionary ####
# Define kmer size controller 
kmerLen = 19				

# Do operations
# Program progress text including start time for elapsed time for step
sys.stdout.write("\nIndexing resistance database...")
start_resdb_time = time.time()

# File opening and error handling
try: 
	resdb = open(resFile, "r")
except IOError as err:
	sys.stderr.write("Error when loading file: {}\n".format(str(err)))
	sys.exit(1)

# Extract headers and DNA sequences
resLine = "nil"

# The line containing resistance sequence
while resLine != "" and resLine[0] != ">":
	resLine = resdb.readline()

sequences, headers = [], []
while resLine != "":
	# Extract headers
	header = None
	if resLine != "" and resLine[0] == ">":
		header = resLine.rstrip()
	
	# Extract DNA sequences
	seq = ""
	resLine = resdb.readline()
	while resLine != "" and resLine[0] != ">":
		resLine = re.sub("\s", "", resLine).upper()

		if re.search(r"[^ATCGN]", resLine) is not None:
			sys.stderr.write("\nIdentified non-nucleotide sequence: {}\n".format(line))
			sys.exit(1)

		seq += resLine	
		resLine = resdb.readline()

	# Store sequence in list and its header in another list (same indices)
	sequences += [seq]
	headers += [header]

	# Check that the lists are equal in length, otherwise terminate!
	if len(sequences) != len(headers):
		sys.stderr.write("\nUnequal number of headers and sequences.... Terminating!\n")
		sys.exit(1)
	
resdb.close()

# Generates a resistance database of kmers (key) linked to a list of associated indices (values) in sequences and headers
kmerDict = dict()			# {sequence:header} dictionary 
for i,seq in enumerate(sequences): 
	kmerPos = 0
	for kmer in yield_kmer(seq, kmerLen):
		if kmer not in kmerDict:
			kmerDict[kmer] = [(i,kmerPos)]
		else:
			kmerDict[kmer] += [(i,kmerPos)]
		kmerPos += 1

sys.stdout.write("\nResistance database indexing took {:.3f} seconds.".format((time.time() - start_resdb_time)))
########################################################################################
####  Work with reads (gzipped files) and find depth at each position in sequence ####

# Define final depth dictionary storing depth of each base in resistance sequence 
depthDict = dict()

# Do operations
# Program progess text including start time of elapsed time for step
sys.stdout.write("\nObtaining reads and scanning for resistances...")
start_scan_time = time.time()

# Iterate over each file provided
for filename in readFiles:
	try:
		with gzip.open(filename, "rt") as readFile:
			
			readSeq = ""	# Read sequence string 
			lineNum = 3		# Reads on every 4th line
			for line in readFile:
				if lineNum == 4:
					readSeq = line.rstrip()
					lineNum = 1	

					# Check if kmers in readSeq
					found = check_kmers(readSeq, kmerLen)
					if not found:
						readSeq = reverse_complement(readSeq)
						found = check_kmers(readSeq, kmerLen)

					# If kmers are found in readSeq increment depth counter of temporary dict
					if found:
						temp_depthDict = dict()
						for kmer in yield_kmer(readSeq, kmerLen):
							if kmer in kmerDict:
								# "tup" is tuple from resistance dict where tup[0] = gene index and tup[1] = kmer start position
								for tup in kmerDict[kmer]:
									if tup[0] not in temp_depthDict:
										# Initialize an index key linked to a list which has length(sequence) and is filled with 0s.  
										temp_depthDict[tup[0]] = [0]*len(sequences[tup[0]])
									# Add 1 depth at observed position (tup[1] = kmer start pos)
									for i in range(tup[1], tup[1]+kmerLen):
										temp_depthDict[tup[0]][i] = 1


						# Add contents of temporary depth dict to final depth dict
						for index in temp_depthDict:
							if index not in depthDict:
								depthDict[index] = temp_depthDict[index]
							elif index in depthDict:
								for i,depth in enumerate(depthDict[index]):
									depthDict[index][i] += temp_depthDict[index][i]

				else: 
					lineNum += 1

	except IOError as err: 
		sys.stderr.write("\nLoad file error: {}\n".format(str(err)))
		sys.exit(1)

sys.stdout.write("\nObtaining reads and resistance scan took {:.3f} seconds.".format(time.time()-start_scan_time))

########################################################################################
#### Determine average depth and coverage, and write to file ####

# Define variables
depthThreshold = 10 # depth threshold controller
covThreshold = 0.95	# coverage threshold controller (max 1.0)
numHits = 0

# Program progress text and start time for elapsed time of step
sys.stdout.write("\nWriting to file...")
start_write_time = time.time()

# Do operations 
try:
	outfile = open("resistances.fsa", "w")
	covDict = dict()
	avgDepthDict = dict()

	for index in depthDict:
		sequence = sequences[index]
		
		# Calculate average depth and coverage
		sumDepth = 0
		depthHit = 0

		for depth in depthDict[index]:
			sumDepth += depth
			if depth >= depthThreshold:
				depthHit += 1
		avgDepth = sumDepth // len(sequence) 
		coverage = depthHit / len(sequence)

		avgDepthDict[index] = avgDepth

		# Write positive hits above a certain threshold to file
		if coverage >= covThreshold:
			numHits += 1
			covDict[index] = coverage

	valuesToOutput = []
	for index in covDict:

		header = headers[index]
		sequence = sequences[index]

		avgDepth = avgDepthDict[index]
		coverage = covDict[index]

		valuesToOutput += [{
			'header': header,
			'sequence': sequence,
			'coverage': coverage,
			'avgDepth': avgDepth
		}]

	valuesToOutput = sorted(valuesToOutput, key = lambda i: (i['coverage'], i['avgDepth']), reverse = True)

	for toWrite in valuesToOutput:

		gene = toWrite['header'].split("_")[0].replace(">", "")
		res = re.search(r"\s(.+):", toWrite['header']).group(1)
		
		outfile.write(">{} | Gene: {} | cov: {:.3f} | avg depth: {}\n".format(res, gene, toWrite['coverage'], toWrite['avgDepth']))
		for i in range(0, len(toWrite['sequence'])-1,60):
			outfile.write(toWrite['sequence'][i:i+60] + "\n")	

				
except IOError as err:
	sys.stderr.write("\nError while writing to file: {}\n".format(str(err)))
	sys.exit(1)

except KeyError as err:
	sys.stderr.write("\nError! {} not found in dictionary\n".format(str(err)))
	sys.exit(1)

outfile.close()
sys.stdout.write("\nWriting to file took {:.3f} seconds.\nIdentified {} resistance gene(s)!\nResistance sequences logged in \"resistances.fsa\"\n".format(time.time()-start_write_time, numHits))
########################################################################################