
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

	return(dico, seq)


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
	print("--------------------")
	return dico, freq_dico


######### SUR PLUSIEURS FICHIERS FASTAS
import os
import sys

def specific_occurence(path) :
	#print(path)
	fich = os.listdir(path)
	#print(fich)
	#print(type(fich))

	reads = []
	for fasta in fich :
		#print(type(fasta))
		print(fasta)
		reads.append(read_fasta(fasta))

	print(reads)

	'''
	for fasta in fich : 
		with open(fasta, 'r') as fastin :
	'''




if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/Q9SNB7.fasta")
	#dico_number, frequency = freq_aa(sequence)

	# Ensemble de fichiers
	proportion = specific_occurence("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta")


