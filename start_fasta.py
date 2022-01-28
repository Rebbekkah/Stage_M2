import os
import sys
import glob
import pandas as pd
import numpy as np 


path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/"

def Score_aa() :
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']
	Z1 = [-2.49, 2.18, 0.07, -4.19, 1.96, -4.44, -1.22, 2.84, 2.23, -2.69, 2.88, 3.08, -4.92, 3.64, 0.71, 0.92, \
			3.22, -4.75, -1.39, 2.41]
	Z2 = [-0.27, 0.53, -1.73, -1.03, -1.63, -1.68, 0.88, 1.41, -5.36, -2.53, 2.52, 0.39, 1.30, 1.13, -0.97, -2.09, \
			1.45, 3.65, 2.32, 1.74]
	Z3 = [-0.41, -1.14, 0.09, -0.98, 0.57, -1.03, 2.23, -3.14, 0.30, -1.29, -3.44, -0.07, 0.45, 2.36, 4.13, -1.40, \
			0.84, 0.85, 0.01, 1.11]
	df_Zscore = pd.DataFrame({'AA' : list_of_aa, 'Z1' : Z1, 'Z2' : Z2, 'Z3' : Z3})
	return df_Zscore

######### SUR 1 FICHIER FASTA

def read_fasta(fichier) :
	#fich = list(fichiers)
	#for f in fich :

	with open(fichier, 'r') as filin :

		dico = {}
		seq = ""
		for line in filin :
			if line.startswith('>') :
				prot_id = line.split()[0]
				dico[prot_id] = ""
			else :
				dico[prot_id] += line.strip()
		
	seq += dico[prot_id]
	#seq.split()
	print(seq)
	print("Longueur de la séquence : ", len(seq))
	print("--------------------")

	return(dico) ######### Peut être return juste le dico pour la suite


def freq_aa(sequence) : 
	dico = {}
	for letter in sequence :
		if letter not in dico.keys() :
			dico[letter] = 1
		else :
			dico[letter] += 1

	total_aa = sum(dico.values()) #ou len(sequence)

	freq_dico = {}
	for cle, value in dico.items() :
		if cle not in freq_dico.keys() :
			freq_dico[cle] = value/total_aa


	print("Amino acid and their occurence : ", dico)
	print("And their frequency : ", freq_dico)
	print(len(freq_dico.keys()))
	print("----------freq_aa----------")
	return freq_dico


######### SUR PLUSIEURS FICHIERS FASTAS

def listing(path) :
	fich = glob.glob(path+'*.fasta')

	reads = []
	for fasta in fich :
		reads.append(read_fasta(str(fasta)))
	#print(reads)
	return reads


def specific_occurence(dico) :
	#print(path)
	#fich = os.listdir(path)
	#fich = os.listdir(path+'*.fasta')
	#fich = glob.glob('path+*.fasta')
	#fich = glob.glob(path+'*.fasta')
	#print("path et fasta ---> ", 'path'+'*.fasta')
	#print(fich)
	#print(type(fich))
	#print("----------specific_occ----------")
	reads = listing(path)
	#reads = []
	'''
	for fasta in fich :
		#print("ERREUR")
		print(type(fasta))
		print(fasta)
		reads.append(read_fasta(str(fasta)))
	'''
	print(reads)
	print(type(reads))

	frequencies = []
	for dicct in reads :
		#print(type(dicct))
		for value in dicct.values() :
			#print(value)
			frequencies.append(freq_aa(value))
			#[dicct.keys()]
	#print(frequencies)

	print("----------DF FREQ-----------")
	list_df = []
	for dicct in frequencies :
		dicct = pd.DataFrame([dicct])
		#print(type(dicct))
		print(dicct)
		list_df.append(dicct)

	print(list_df, type(list_df))

	return list_df


def Z_aa(Z_score) :
	dico = listing(path)
	print(dico)

	
	for dicct in dico :
		for seq in dicct.values() :
			#print(len(seq))
			for aa in seq :
				

def Auto_cross_variance() :
	pass
	#ACCz_same = sum()



######### Biopython

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#def bio(path) :
#def bio(*fasta) :
def bio() :
	'''
	fich = os.listdir(path)
	print(fich)
	print("------------")
	for fasta in fich :
		for seq_record in SeqIO.parse(fasta, "fasta") :
			print(seq_record.id)
			print(seq_record.seq)
			print(len(seq_record))
	'''

	'''
	files = list(fasta)
	print(files)
	frequencies = []
	#files = os.listdir(path)
	#for fasta in files :
	for file in files :
		for seq_record in SeqIO.parse(file, "fasta") :
			my_seq = ProteinAnalysis(str(seq_record.seq))
			#print(my_seq)
			#print(my_seq)
			#print(my_seq.id)
			#print(my_seq.seq)
			#print(len(my_seq))
		
			freq = my_seq.get_amino_acids_percent()
			#count = my_seq.count_amino_acids()['E']
			print(type(freq))

			frequencies.append(freq)
			print(frequencies)
			#return freq

	'''

	for seq_record in SeqIO.parse("Q9SNB7.fasta", "fasta") :
		record = SeqIO.read("Q9SNB7.fasta", "fasta")
		print(record, "\n -----------------------")
		seq = record.seq
		print(seq)
		my_seq = ProteinAnalysis(str(seq_record.seq))
		#print(my_seq)
		#print(my_seq)
		#print(my_seq.id)
		#print(my_seq.seq)
		#print(len(my_seq))
	
		freq = my_seq.get_amino_acids_percent()
		print(freq)


		#print(seq_record.annotations)

		return freq


def Tsne(frequencies) :
	pass





if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)

	df_Score = Score_aa()

	# Ensemble de fichiers
	proportion = specific_occurence(path)
	ACC_score = Z_aa(df_Score)


	# Biopython
	#seq = bio()













