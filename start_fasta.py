path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/"

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
	print("----------freq_aa----------")
	return freq_dico


######### SUR PLUSIEURS FICHIERS FASTAS
import os
import sys
import glob
import pandas as pd
import numpy as np 

def specific_occurence(path) :
	#print(path)
	#fich = os.listdir(path)
	#fich = os.listdir(path+'*.fasta')
	#fich = glob.glob('path+*.fasta')
	fich = glob.glob(path+'*.fasta')
	#print("path et fasta ---> ", 'path'+'*.fasta')
	print(fich)
	print(type(fich))
	print("----------specific_occ----------")

	reads = []
	for fasta in fich :
		#print("ERREUR")
		print(type(fasta))
		print(fasta)
		reads.append(read_fasta(str(fasta)))

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


def Auto_cross_variance() :
	pass



if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)

	# Ensemble de fichiers
	proportion = specific_occurence(path)

	# Biopython
	#seq = bio()













