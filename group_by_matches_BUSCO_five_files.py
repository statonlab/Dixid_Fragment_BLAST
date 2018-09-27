#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

##---------------------------------
## FUNTIONS
##---------------------------------

##---------------------------------
def create_match_files():
	match_seqs = SeqIO.to_dict(SeqIO.parse(open(match_fasta_file), "fasta"))
	for match in common_matches:
		with open(match + ".fasta", "w") as output_handle:
			if match in match_seqs:
				SeqIO.write(match_seqs[match], output_handle, "fasta")
			else:
				print("Missing sequence for match: " + match)
	return

##---------------------------------
def add_query(blast_file, query_fasta_file):

	# load fasta file into memory
	query_seqs = SeqIO.to_dict(SeqIO.parse(open(query_fasta_file), "fasta"))

	# parse through blast file and add to fasta files
	with open(blast_file) as f:
		for line in f:
			fields = line.split('\t')
			query = fields[0]
			match = fields[1]
			if match in common_matches:
				with open(match + ".fasta", "a") as output_handle:
					if query in query_seqs:
						SeqIO.write(query_seqs[query], output_handle, "fasta")
					else:
						print("Missing sequence for query: " + query)

	return

##---------------------------------
##---------------------------------
##---------------------------------
## This script builds individual fasta files containing the sequences associated with every 
## query sequence observed in an input file
## python group_by_matches.py blast_file hits_fasta_file matches_fasta_file

# Best hits only
match_fasta_file = str(sys.argv[1])

blast1 = "DCORN_augustus_ORFs_2_blastx.tsv.best"
query_fasta1 = "DCORN_augustus_ORFs_2.codingseq"

blast2 = "DgermMidge_86M_ORFs_blastx.tsv.best"
query_fasta2 = "DgermMidge_86M_ORFs.codingseq"

blast3 = "Dpuget_52M_ORFs_blastx.tsv.best"
query_fasta3 = "Dpuget_52M_ORFs.codingseq"

blast4 = "K46i6_Dtern_ORFs_blastx.tsv.best"
query_fasta4 = "K46i6_Dtern_ORFs.codingseq"

blast5 = "K47i6_Dalic_ORFs_blastx.tsv.best"
query_fasta5 = "K47i6_Dalic_ORFs.codingseq"

blast6 = "DBLAX_augustus_ORFS_blastx.tsv.best"
query_fasta6 = "DBLAX_augustus_ORFS.codingseq"

# Find the matches shared in all blast files
matches1 = [line.split("\t")[1] for line in open(blast1).readlines()]
matches2 = [line.split("\t")[1] for line in open(blast2).readlines()]
matches3 = [line.split("\t")[1] for line in open(blast3).readlines()]
matches4 = [line.split("\t")[1] for line in open(blast4).readlines()]
matches5 = [line.split("\t")[1] for line in open(blast5).readlines()]
matches6 = [line.split("\t")[1] for line in open(blast6).readlines()]

common_matches = set.intersection(set(matches1), set(matches2), set(matches3), set(matches4), set(matches5), set(matches6))
print(len(common_matches))


## Create the match fasta files, starting with just the match seq
create_match_files()

## Add the query seqs to these fasta files
add_query(blast1, query_fasta1)
add_query(blast2, query_fasta2)
add_query(blast3, query_fasta3)
add_query(blast4, query_fasta4)
add_query(blast5, query_fasta5)
add_query(blast6, query_fasta6)


