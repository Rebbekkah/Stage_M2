import numpy as np 
import pandas as pd

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split


def read_seq(*fasta):

	files = list(fasta)
	#print(files)
	#print(type(files))

	dico = {}

	for file in files :
		for seq_record in SeqIO.parse(file, "fasta") :
			record = SeqIO.read(file, "fasta")
			print(record.id)
			#print(record.seq)

			seq_name = record.id
			#dico[seq_name] = ""

			my_seq = ProteinAnalysis(str(seq_record.seq))
			freq = my_seq.get_amino_acids_percent()
			dico[seq_name] = freq
	print(dico)
	print("----------------------")

	return dico
		
def tsne(dico) :
	

if __name__ == '__main__' :

	dico_freq = read_seq("Q9SNB7.fasta", "F4HXU3.fasta")





