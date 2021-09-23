#!/usr/bin/env python

from sys import argv
from Bio import SeqIO

cds_lines = []
pep_lines = []

for rec in SeqIO.parse(open(argv[1]), "fasta"):
  for i in range(3):
    remainder = len(rec.seq[i:])%3
    if remainder == -0:
      for_cds = rec.seq[i:]
    else:
      for_cds = rec.seq[i:-remainder]

    for_str = str(for_cds.translate(table=1))
    cds_lines.append(">%s_for_%s\n%s\n" % (rec.name, i, for_cds))
    pep_lines.append(">%s_for_%s\n%s\n" % (rec.name, i, for_str.replace("*", "X")))

    rev_remainder = len(rec.seq[:-(i+1)])%3
    rev_cds = rec[rev_remainder:-(i+1)].seq.reverse_complement()
    rev_str = str(rev_cds.translate(table=1))

    cds_lines.append(">%s_rev_%s\n%s\n" % (rec.name, i, rev_cds))
    pep_lines.append(">%s_rev_%s\n%s\n" % (rec.name, i, rev_str.replace("*", "X")))

outfh_cds = open(argv[2], 'w+')
outfh_cds.writelines(cds_lines)
outfh_cds.close()

outfh_pep = open(argv[3], 'w+')
outfh_pep.writelines(pep_lines)
outfh_pep.close()


# vim: set sts=2 ts=2 sw=2:
