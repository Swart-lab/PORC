#!/usr/bin/env python

from Bio import SeqIO, Seq
from collections import OrderedDict
from numpy import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
  '--cds', '-c', help='Fasta file with CDS sequences')
parser.add_argument(
  '--hmmer', '-a', help='Alignments produced by hmmsearch, output file from -o option')
parser.add_argument(
  '--matrix', '-m', help='Output path to write matrix for Weblogo')
parser.add_argument(
  '--codon_threshold', '-t', type=float, default=0.05,
  help='Cut-off for including codon to be displayed; fraction of the mean codon frequency of other non-stop codons')
parser.add_argument(
  '--eval_threshold', '-v', type=float, default=1e-20,
  help='Conditional evalue threshold for HMMER domains to use')
args = parser.parse_args()

if not args.cds or not args.hmmer or not args.matrix:
    print('Options --cds, --hmmer, and --matrix are mandatory')
    parser.print_help()
    exit()


codon_l = ( # standard genetic code
  ('TTT', 'F'),
  ('TTC', 'F'),
  ('TTA', 'L'),
  ('TTG', 'L'),
  ('TCT', 'S'),
  ('TCC', 'S'),
  ('TCA', 'S'),
  ('TCG', 'S'),
  ('TAT', 'Y'),
  ('TAC', 'Y'),
  ('TAA', '*'),
  ('TAG', '*'),
  ('TGT', 'C'),
  ('TGC', 'C'),
  ('TGA', '*'),
  ('TGG', 'W'),
  ('CTT', 'L'),
  ('CTC', 'L'),
  ('CTA', 'L'),
  ('CTG', 'L'),
  ('CCT', 'P'),
  ('CCC', 'P'),
  ('CCA', 'P'),
  ('CCG', 'P'),
  ('CAT', 'H'),
  ('CAC', 'H'),
  ('CAA', 'Q'),
  ('CAG', 'Q'), 
  ('CGT', 'R'),
  ('CGC', 'R'),
  ('CGA', 'R'),
  ('CGG', 'R'),
  ('ATT', 'I'),
  ('ATC', 'I'),
  ('ATA', 'I'),
  ('ATG', 'M'),
  ('ACT', 'T'),
  ('ACC', 'T'),
  ('ACA', 'T'),
  ('ACG', 'T'),
  ('AAT', 'N'),
  ('AAC', 'N'),
  ('AAA', 'K'),
  ('AAG', 'K'),
  ('AGT', 'S'),
  ('AGC', 'S'),
  ('AGA', 'R'),
  ('AGG', 'R'),
  ('GTT', 'V'),
  ('GTC', 'V'),
  ('GTA', 'V'),
  ('GTG', 'V'),
  ('GCT', 'A'),
  ('GCC', 'A'),
  ('GCA', 'A'),
  ('GCG', 'A'),
  ('GAT', 'D'),
  ('GAC', 'D'),
  ('GAA', 'E'),
  ('GAG', 'E'),
  ('GGT', 'G'),
  ('GGC', 'G'),
  ('GGA', 'G'),
  ('GGG', 'G'))
codon_d = OrderedDict(codon_l)

def aln_cds_to_pep_aln(cds, pep, pep_start=145):
  """Align CDS relative to amino acid sequence with gaps

  With inputs: ungapped CDS sequence and gapped amino acid sequence, introduce
  gaps to the CDS sequence so that it is aligned to the amino acids.
  """
  new_cds_l = []
  j = pep_start
  cds_str = str(cds)
  for i in range(len(pep)):
    if pep[i] == '-':
      j -= 1
      new_cds_l.append("---")
    else:
      new_cds_l.append(cds_str[j*3:(j+1)*3])
    j += 1
  return(new_cds_l)


def ok_cod(acod):
  """Check if a codon is valid

  Returns True if not a gap, has length 3, and only contains A, T, C or G,
  False otherwise.
  """
  chars = ['A', 'C', 'G', 'T']
  if acod == '---':
    return(False)
  elif len(acod) == 3 and acod[0] in chars and acod[1] in chars and acod[2] in chars:
    return(True)
  else: 
    return(False)


# Read CDS file to dict, keyed by sequence name
cds_d = {}
for rec in SeqIO.parse(open(args.cds), 'fasta'): # CDS file
  cds_d[rec.name] = rec.seq

cod_d = {}
acod_d = {}

text = open(args.hmmer, 'rt').read() # HMMer output

chunks = text.split("//")
found_both_l = []
found_uar_only_l = []
found_uga_only_l = []

cod_count_d = dict([(cod, 0) for cod in codon_d])
used_cods = 0

for chunk in chunks:
  if "Domain annotation for each sequence" in chunk:
    hit_chunks = chunk.split(">>")[1:]
    for hit_chunk in hit_chunks:
      hit_lines = hit_chunk.split("\n")
      if hit_lines[1][:4] == '   #':
        aln_part = hit_chunk.split("Alignments for each domain:\n")[1]
        for dom_part in aln_part.split("score:")[1:]:

          cond_eval = float(dom_part.split("conditional E-value: ")[1].split("\n")[0])
          if cond_eval < args.eval_threshold:

            aln_lines = dom_part.split("\n")
            hmm_aln = []
            query_aln = []

            #if " RF\n" not in hit_chunk:
            ok_case = False
            if aln_lines[1][-3:] == " CS" and aln_lines[5][-3:] == ' PP':#  == " CS" and aln_lines[2][:-3] == " RF": 
              query_n = 4
              hmm_n = 2
              chunk_aln_size = 6
              ok_case = True
            elif aln_lines[1][-3:] == " CS" and aln_lines[2][-3:] == " RF": #" CS\n" in hit_chunk:
              query_n = 5
              hmm_n = 3
              chunk_aln_size = 7
              ok_case = True
            elif aln_lines[4][-3:] == ' PP':
              query_n = 3
              hmm_n = 1
              chunk_aln_size = 5
              ok_case = True
          
            if ok_case:
              query_aln_line = aln_lines[query_n]
              hmm_aln_line = aln_lines[hmm_n]
              
              #print "QUERY", query_aln_line
              #print "HMM", hmm_aln_line
              #print 

              query_aln_line_atoms = aln_lines[query_n].split()
              hmm_aln_line_atoms = aln_lines[hmm_n].split()

              query_id = query_aln_line_atoms[0]
              
              hmm_id = hmm_aln_line_atoms[0]
              query_aln_start = int(query_aln_line_atoms[1])
              hmm_aln_start = int(hmm_aln_line_atoms[1])

              for i in range(0, len(aln_lines) - chunk_aln_size, chunk_aln_size):
                if len(aln_lines[i+hmm_n]) < 10: break

                hmm_aln_line_atoms = aln_lines[i+hmm_n].split()
                #try:
                hmm_aln.append(hmm_aln_line_atoms[2])
                #except IndexError:
                #  print "XX", aln_lines[i+hmm_n], "XX"
                #  print hit_chunk
                #  break

                query_aln_line_atoms = aln_lines[i+query_n].split()
                query_aln.append(query_aln_line_atoms[2])

              query_aln_seq = "".join(query_aln)
              hmmer_cons = "".join(hmm_aln)

              cds_aln_l = aln_cds_to_pep_aln(cds_d[query_id], query_aln_seq, pep_start=query_aln_start-1)
              
              found_uar = False
              found_uga = False

              for h, q, c in zip(hmmer_cons, query_aln_seq, cds_aln_l):
                if ok_cod(c):
                  cod_count_d[c] += 1

                #h = h.upper() #uncomment to examine all positions
                if q == 'X':
                  #cod_d.setdefault(c, []).append(h.upper())
                  if h.isupper():
                    used_cods += 1
                    cod_d.setdefault(c, []).append(h)

                    if c == 'TGA': found_uga = True
                    elif c == 'TAA' or c == 'TAG': found_uar = True
                h_upper = h.upper()
                if h.isupper():
                  used_cods += 1
                  acod_d.setdefault(c, []).append(h_upper)

              if found_uar and found_uga:
                found_both_l.append((hmm_id, query_id))
              elif found_uar and not found_uga:
                found_uar_only_l.append((hmm_id, query_id))
              elif not found_uar and found_uga:
                found_uga_only_l.append((hmm_id, query_id))

for hmm_id, ref_id in found_both_l:
  print("found both UGA and UAR codons", hmm_id, ref_id)

print() 

for hmm_id, ref_id in found_uar_only_l:
  print("found UAR only", hmm_id, ref_id)

print

for hmm_id, ref_id in found_uga_only_l:
  print("found UGA only", hmm_id, ref_id)

aa_list = "A C D E F G H I K L M N P Q R S T V W Y".split()

stop_hist = {'TGA': dict([(aa, 0) for aa in aa_list]), 'TAA': dict([(aa, 0) for aa in aa_list]), 'TAG': dict([(aa, 0) for aa in aa_list])}
for cod in cod_d:
  for aa in cod_d[cod]:
    stop_hist.setdefault(cod, {})
    stop_hist[cod].setdefault(aa, 0)
    stop_hist[cod][aa] += 1

all_hist = OrderedDict([(cod, dict([(aa, 0) for aa in aa_list])) for cod in codon_d])
aa_count_d = dict([(aa, 0) for aa in aa_list])
for cod in acod_d:
  for aa in acod_d[cod]:
    all_hist.setdefault(cod, {})
    all_hist[cod].setdefault(aa, 0)
    all_hist[cod][aa] += 1
    aa_count_d[aa] += 1

aa_count_tot = sum(aa_count_d.values())
aa_freq_d = dict([(aa, float(aa_count_d[aa])) for aa in aa_list])


stop_cods = ('TAA', 'TAG', 'TGA')

outlines = []
outlines.append("PO" + "\t" + "\t".join(aa_list) + "\n")


stop_tot = {'TGA': 0, 'TAA': 0, 'TAG': 0}
for cod in stop_hist:
  newl = [(v, k) for k, v in stop_hist[cod].items()]
  newl.sort(reverse=True)

  if cod in stop_cods:
    stop_tot[cod] = sum([v for v in stop_hist[cod].values()])

other_tot = dict([(cod[0], 0) for cod in codon_l])
for cod, aa in codon_l:
  other_tot[cod] = sum([v for v in all_hist[cod].values()])
    

cod_tot = sum(list(cod_count_d.values()))

print ()
print("UAA: %s UAG: %s UGA: %s" % (stop_tot['TAA'], stop_tot['TAG'], stop_tot['TGA']))
print("other codons:", mean(list(other_tot.values())), median(list(other_tot.values())))
print('used codons', used_cods)

print()

# Summary statistics on codon frequencies in alignments
for cod in codon_d: # codon sequence
  print("%s\t" % cod, end='')
print()
for cod in codon_d: # amino acid in the standard genetic code (for reference)
  print("%s\t" % codon_d[cod], end='')
print()
for cod in codon_d: # counts of this codon in HMMer alignments
  print("%s\t" % cod_count_d[cod], end= '')
print()
for cod in codon_d: # as percentage of total codon counts
  print("%.1f\t" % (round(100*float(cod_count_d[cod])/cod_tot, 1)), end='')

# Write frequency of observed amino acid per codon to matrix for Weblogo to plot
i = 0
for cod in codon_d:
  i += 1
  if cod in stop_cods:
    # if codon frequency is over the specified threshold, report in count
    # matrix for Weblogo plotting
    if stop_tot[cod] > args.codon_threshold * mean(list(other_tot.values())):
      outlines.append(
        "\t".join([str(i)] + [str(float(all_hist[cod][aa])/aa_freq_d[aa]) for aa in aa_list]) + "\n")
      #outlines.append("%s" % i + "\t" + "\t".join([str(float(all_hist[cod][aa])) for aa in aa_list]) + "\n")
    # else append zeroes (column not displayed by Weblogo)
    else:
      outlines.append(
        "\t".join([str(i)] + ["0" for aa in aa_list]) + "\n")
  else:
    outlines.append(
      "\t".join([str(i)] + [str(float(all_hist[cod][aa])/(aa_freq_d[aa] + 0.00001)) for aa in aa_list]) + "\n")
    #outlines.append("%s" % i + "\t" + "\t".join([str(float(all_hist[cod][aa]) + 0.00001) for aa in aa_list]) + "\n")

outfh = open(args.matrix, 'w+') # Weblogo matrix file
outfh.writelines(outlines)
outfh.close()

  
# vim: set sts=2 ts=2 sw=2:
