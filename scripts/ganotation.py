#!/usr/bin/env python

import sys
from jmolbiol_2 import format_seq, \
					get_gc, \
					translate, \
					find_orfs, \
					blasting

test = format_seq(sys.argv[2], sys.argv[1])
if test != 1:
	sys.exit("Error in format_seq function, it does not return 1")
else:
	print "Successful! The fasta file %s was correctly formated" % sys.argv[2]

scodon = sys.argv[5].split("|")
print "Start codons: ",
print scodon
ecodon = sys.argv[6].split("|")
print "End codons: ",
print ecodon

print "ORFs' quest..."
protfile = find_orfs(sys.argv[2], sys.argv[3], sys.argv[1], int(sys.argv[4]), \
scodon, ecodon)
print protfile

