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

	"""
	gc = {}
	IN = open(gcfilename, "r")
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
	for i in xrange(len(aas)):
		gc[b1[i]+b2[i]+b3[i]] = aas[i]
	
	IN.close()
	
	return gc

def translate(seq, gc):
	"""
	
	"""
	pseq = ""
	for i in xrange(0, len(seq)-3, 3):
		pseq = pseq + gc[seq[i:i+3].upper()]
	
	return pseq

def find_orfs(seqfilename, gcfilename, spath, sizep, scodon, ecodon):
	"""
	
	"""
	gc = get_gc(gcfilename)
	SFH = open("%s%s" % (spath, seqfilename), "r")
	NOFH = open("%s%s_nuc.fasta" % (spath, seqfilename.split(".")[0]), "w+")
	POFH = open("%s%s_prot.fasta" % (spath, seqfilename.split(".")[0]), "w+")
	FOFH = open("%s%s_feat.txt" % (spath, seqfilename.split(".")[0]), "w+")
	FOFH.write("ORFname\tChrName\tStart\tend\tfilename\n")
	orf_counter = 0
	pivot = 1000
	name = ""
	end = 0
	for line in SFH:
		line = line.strip()
		if line[0] == ">":
			name = line.split(".")[0]
			chrname = line
		else:
			i = 0
			while i <= len(line)-3:  
				if line[i:i+3].upper() in scodon:
					seq = ""
					pseq = ""
					for j in xrange(i, len(line)-3, 3):
						seq = seq + line[j:j+3]
						if sizep != None and j > i+(sizep * 3):
							i = i + 1
							break
						else:
							pass
						if line[j:j+3].upper() in ecodon:
							start = i
							end = j+3
							seq = seq + line[j:j+3]
							if len(seq) < 300:
								i = i + 1
								break
							orf_counter = orf_counter + 1
							i = i + (((end + 1) - i) / 4)
							NOFH.write("%s_ORF%s|len_%s\n%s\n" % (name, \
							orf_counter, len(seq), seq))
							pseq = translate(seq, gc)
							POFH.write("%s_ORF%s|len_%s\n%s\n" % (name, \
							orf_counter, len(pseq), pseq))
							FOFH.write("%s_ORF%s|len_%s\t%s\t%s\t%s\t%s\n" %\
							 (name, orf_counter, len(seq), chrname, start, end,\
							 seqfilename))
							if orf_counter == pivot:
								pivot = pivot + 1000
								print "%s ORFs have been found!" % orf_counter
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

	"""
	IN = open("%s%s" % (spath, seqfilename), "r")
	TMP = open("%s_tmp" % spath, "w+")
	TMP.write("%s" % IN.readline())
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			TMP.write("\n%s\n" % line)
		else:
			TMP.write(line)
	IN.close()
	TMP.close()

	IN = open("%s%s" % (spath, seqfilename), "w+")
	TMP = open("%s_tmp" % spath, "r")
	for line in TMP:
		IN.write(line)
	IN.close()
	TMP.close()

	return 1

def blasting(protfile, protdb):
	"""

	"""
	cmd = "blastall -p blastp -d protdb -i protfile -e 1E-20 -m 9 -o %s.blast -n T" % protfile.strip().split(".")[0]
	print "the output file is: %s.blast" % protfile.strip().split(".")[0]
	print "the proteins file is %s: " % protfile
	print "the database is: %s" % protdb
	print "Blast in process..."
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p.wait()
	return "%s.blast" % protfile.strip().split(".")[0]

def 

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
