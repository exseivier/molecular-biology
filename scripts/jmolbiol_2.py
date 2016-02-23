#!/usr/bin/env python


import sys
import subprocess

def extract_seq(datafilename, spath, outpufile):
	"""
	
	"""
	DFH = open(datafilename, "r")
	OFH = open("%s%s" % (spath, outpufile), "w+")
	for line in DFH:
		tag = None
		line = line.strip()
		items = line.split("\t")
		SFH = open("%s%s" % (spath, items[4]), "r")
		for line1 in SFH:
			if tag == False:
				break
			else:
				pass
			line1 = line1.strip()
			if line1[0] == ">" and items[1] in line1:
				tag = True
				continue
			else:
				pass
			if tag == True:
				seq = ""
				for i in xrange(len(line1)):
					if int(items[2]) <= i and int(items[3]) >= i:
						seq = seq + line1[i]
					elif int(items[3]) < i:
						OFH.write(">%s\n%s\n" % (items[0], seq))
						tag = False
						break
					else:
						pass
			else:
				pass
	
	return 1

def get_gc(gcfilename):
	"""
	Returns a data structure holding the genetic code (dictionary variable).
	key -> codon : value -> amino acid letter.
	Requires:
		o gcfilename: the name of a file which contains the genetic code
	"""
	gc = {}
	IN = open(gcfilename, "r")
	# Set variables.
	for line in IN:
		if "AAs" in line:
			aas = line.strip().split("=")[1][1:]
		elif "Base1" in line:
			b1 = line.strip().split("=")[1][1:]
		elif "Base2" in line:
			b2 = line.strip().split("=")[1][1:]
		elif "Base3" in line:
			b3 = line.strip().split("=")[1][1:]
		else:
			pass
	# Building dictionary.
	for i in xrange(len(aas)):
		gc[b1[i]+b2[i]+b3[i]] = aas[i]
	
	IN.close()
	
	return gc

def translate(seq, gc):
	"""
	Returns a string variable which holds the translated sequence of
	amino acids residues of an open reading frame (ORF).
	Requires:
		o seq: the ORF.
		o gc: genetic code data structure (dictionary).
	"""
	pseq = ""
	for i in xrange(0, len(seq), 3):
		# Adding the amino acid based on gc dictionary data.
		pseq = pseq + gc[seq[i:i+3].upper()]
	
	return pseq

def find_orfs(seqfilename, gcfilename, spath, sizep, scodon, ecodon):
	"""
	Finds ORFs over a sequence.
	Requires:
		o seqfilename: sequences file name (str).
		o gcfilename: genetic code file name (str).
		o spath: path to sequences file (str).
		o sizep: maximal protein size to search (int).
		o scodon: start codons (array:str).
		o ecodon: stop codons (array:str).
	"""
	gc = get_gc(gcfilename)
	the_last_word = ""
	SFH = open("%s%s" % (spath, seqfilename), "r")
	NOFH = open("%s%s_nuc.fasta" % (spath, seqfilename.split(".")[0]), "w+")
	POFH = open("%s%s_prot.fasta" % (spath, seqfilename.split(".")[0]), "w+")
	FOFH = open("%s%s_feat.txt" % (spath, seqfilename.split(".")[0]), "w+")
	FOFH.write("ORFname\tChrName\tStart\tend\tfilename\n")
	orf_counter = 0
	pivot = 1000
	name = ""
	end = 0
	# Iterating over sequences on sequences file.
	for line in SFH:
		line = line.strip()
		if line[0] == ">":
			name = line.split(".")[0]
			chrname = line
		else:
			i = 0
			# Iterating by nucleotide over the sequence.
			while i <= len(line)-3:  
				if line[i:i+3].upper() in scodon:
					start = i
					seq = ""
					pseq = ""
					# A start codon is found.
					# Aditional sliding window. Iterating by nucleotides
					# over the sequence until a stop codon is found
					for j in xrange(i, len(line)-3, 3):
						seq = seq + line[j:j+3]
						# If protein size is greater than sizep * 3
						# break
						if sizep != None and j > i+(sizep * 3):
							i = i + 1
							break
						else:
							pass
						# A stop codon is found
						if line[j:j+3].upper() in ecodon:
							end = j+3
							# If this new ORF ends like the last ORF
							# break.
							if the_last_word == line[end-20:end]:
								i = i + 1
								break
							# If protein size is lower than 300
							# break
							if len(seq) < 300:
								i = i + 1
								break
							# Here, the_last_word is updated
							the_last_word = line[end-20:end]
							# Here orf_counter is updated
							orf_counter = orf_counter + 1
							# Here every output file is updated
							NOFH.write("%s_ORF%s|len_%s\n%s\n" % (name, \
							orf_counter, len(seq), seq))
							pseq = translate(seq, gc)
							POFH.write("%s_ORF%s|len_%s\n%s\n" % (name, \
							orf_counter, len(pseq), pseq))
							FOFH.write("%s_ORF%s|len_%s\t%s\t%s\t%s\t%s\n" %\
							 (name, orf_counter, len(seq), chrname, start, end,\
							 seqfilename))
							# Printing results to screen
							if orf_counter == pivot:
								pivot = pivot + 1000
								print "%s ORFs have been found!" % orf_counter
							i = i + 1
							break
						else:
							i = i + 1
					i = i + 1
				else:
					i = i + 1
	SFH.close()
	NOFH.close()
	POFH.close()
	FOFH.close()
	print "Total ORFs found were %s" % orf_counter
	return "%s%s_prot.fasta" % (spath, seqfilename.split(".")[0])


def format_seq(seqfilename, spath):
	"""
	Re-shapes the sequences format. Removes the \n return special character
	from the end of the found sequence.
	Requires:
		o seqfilename: sequences file name
		o spath: path to sequences file
	"""
	IN = open("%s%s" % (spath, seqfilename), "r")
	TMP = open("%s_tmp" % spath, "w+")
	TMP.write("%s" % IN.readline())
	# Formating the sequences.
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			# Storing them in a temporary file.
			TMP.write("\n%s\n" % line)
		else:
			TMP.write(line)
	IN.close()
	TMP.close()
	# Transferring the formated sequences to the original file.
	IN = open("%s%s" % (spath, seqfilename), "w+")
	TMP = open("%s_tmp" % spath, "r")
	for line in TMP:
		IN.write(line)
	IN.close()
	TMP.close()

	return 1

def blasting(protfile, protdb):
	"""
	Blasting the hypothetical proteome to a database of proteins.
	Requires:
		o protfile: the hypothetical proteome (translated ORFs found in genome) file name
		o protdb: the protein database file name
		ERROR:
			This function sends the process to the shell and it does not wait for process end,
			it returns to the main function and resumes all the python script.
			If a blast output is not found, as is in this case, the script asks for it and
			an error raises crashing the script.
	"""
	cmd = "blastall -p blastp -d protdb -i protfile -e 1E-20 -m 9 -o %s.blast" % protfile.strip().split(".")[0]
	print "the output file is: %s.blast" % protfile.strip().split(".")[0]
	print "the proteins file is %s: " % protfile
	print "the database is: %s" % protdb
	print "Blast in process..."
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	print p.wait()
	return "%s.blast" % protfile.strip().split(".")[0]

def parsingBlastResult(filename, cutoff):
	"""
	
	"""
	IN = open(filename, "r")
	OUT = open("%s.filter" % filename.strip()[:-6], "w+")
	bit_score = 0
	for line in IN:
		line = line.strip()
		if line[0] != "#":
			items = line.split("\t")
			if bit_score == 0:
				OUT.write("%s\t%s\n" % (items[0], items[1]))
				bit_score = int(items[11])
				continue
			elif bit_score > 0 and ((int(items[11]) * 100) / bit_score) >= cutoff:
				OUT.write("%s\t%s\n" % (items[0], items[1]))
				bit_score = int(items[11])
				continue
			else:
				continue
		else:
			continue
	IN.close()
	OUT.close()
	return "%s.filter" % filename.strip()[:-6]


def doc_help():
	"""

	"""
	return """
	Under construction
	Molecular Biology Library
	Dedicated to sequence data analysis
	"""

if __name__ == "__main__":
	print doc_help()
