import os
import sys
import glob
import pandas as pd
import numpy as np 
import seaborn as sns
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt

# R
import rpy2.robjects as robjects
r = robjects.r
'''
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind = 1)

# Install packages
packnames = ('protr')
utils.install_packages(StrVector(packnames))
# Load packages 
protr = rpackages.importr('protr')
'''


# important
path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/proteome_Chlamydomonas.fa"
path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes"
list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']


def Score_aa() :
	Z1 = [-2.49, 2.18, 0.07, -4.19, 1.96, -4.44, -1.22, 2.84, 2.23, -2.69, 2.88, 3.08, -4.92, 3.64, 0.71, 0.92, \
			3.22, -4.75, -1.39, 2.41]
	Z2 = [-0.27, 0.53, -1.73, -1.03, -1.63, -1.68, 0.88, 1.41, -5.36, -2.53, 2.52, 0.39, 1.30, 1.13, -0.97, -2.09, \
			1.45, 3.65, 2.32, 1.74]
	Z3 = [-0.41, -1.14, 0.09, -0.98, 0.57, -1.03, 2.23, -3.14, 0.30, -1.29, -3.44, -0.07, 0.45, 2.36, 4.13, -1.40, \
			0.84, 0.85, 0.01, 1.11]
	df_Zscore = pd.DataFrame({'AA' : list_of_aa, 'Z1' : Z1, 'Z2' : Z2, 'Z3' : Z3})
	df_Zscore.set_index('AA', inplace = True)

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
	
	print(seq)
	print("Longueur de la séquence : ", len(seq))
	print("--------------------")

	dico_delete = {}	
	not_aa = []
	id_not_aa = []

	'''
	for idt, seq in dico.items() :
		if seq[-1] != '*' :
			for aa in seq :
				if aa not in list_of_aa :
					not_aa.append(aa)
					id_not_aa.append(idt)
					break
		elif seq[-1] == '*' :
			for aa in seq[:-1] :
				if aa not in list_of_aa :
					not_aa.append(aa)
					id_not_aa.append(idt)
					break

	'''
	for idt, seq in dico.items() :
		if seq[-1] != '*' :
			for aa in seq :
				if aa not in list_of_aa :
					if aa not in not_aa :
						not_aa.append(aa)	
					if idt not in id_not_aa :				
						id_not_aa.append(idt)

		elif seq[-1] == '*' :
			for aa in seq[:-1] :
				if aa not in list_of_aa :
					if aa not in not_aa :
						not_aa.append(aa)
					if idt not in id_not_aa :
						id_not_aa.append(idt)
						
		for ident in id_not_aa :
			dico_delete[ident] = ""
			dico_delete[ident] += seq

	print("----NOT------")
	print(not_aa, len(not_aa))
	print(id_not_aa, len(id_not_aa))
	#print(dico_delete)
	#nb_delete = len(dico_delete.keys())
	nb_delete = len(id_not_aa)
	print("NOMBRE DE SEQ DELETE : ", nb_delete)

	
	with open("Deleted_seq", "w") as filout :
		for idt, seq in dico_delete.items() :
			#filout.write(idt+"\n"+seq+"\n")
			filout.write(idt+"\n")
	
	print("avant :", len(dico.keys()), len(dico.values()))

	dico_nb = {}
	for idt, seq in dico.items() :
		#print("-----")
		dico_nb[idt] = parsing_seq(seq)

	dico_total = {}

	for dic in dico_nb.values() :
		for cle, val in dic.items() :
			if cle not in dico_total.keys() :
				dico_total[cle] = 0
			else :
				dico_total[cle] += val

	print("TOTAL :", dico_total)


	for idt in id_not_aa :
		del dico[idt]

	print("après :", len(dico.keys()), len(dico.values()))
	
	with open("New_proteome", "w") as filout :
		for idt, seq in dico.items() :
			#filout.write(idt+"\n"+seq+"\n")
			filout.write(idt+"\n")

	return(dico) 


def parsing_seq(seq) :
	
	dico = {'codon-stop' : 0, 'pseudo-gène' : 0, 'other' : 0}
	
	if seq[-1] == '*' :
		dico['codon-stop'] += 1
		for aa in seq[:-1] :
			if aa == '*' :
				dico['pseudo-gène'] += 1
			elif aa not in list_of_aa :
				dico['other'] += 1
			elif aa in list_of_aa :
				continue
			else :
				break

	elif seq[-1] != '*' :
		for aa in seq :
			if aa == '*' :
				dico['pseudo-gène'] += 1
			elif aa not in list_of_aa :
				dico['other'] += 1
			elif aa in list_of_aa :
				continue
			elif (aa == '*') or (aa not in list_of_aa) :
				break


	return dico


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


	for aa in list_of_aa :
		if aa not in freq_dico :
			freq_dico[aa] = 0

	#print("Amino acid and their occurence : ", dico)
	#print("And their frequency : ", freq_dico)
	#print(len(freq_dico.keys()))
	#print("----------freq_aa----------")
	

	return freq_dico


######### SUR PLUSIEURS FICHIERS FASTAS

def listing(path) :
	#fich = glob.glob((path+'*.fasta')||(path+'*.faa')||(path+'*.fna'))
	fich = glob.glob(path+'*.fasta')

	reads = []
	for fasta in fich :
		reads.append(read_fasta(str(fasta)))
	#print(reads)
	return reads


def specific_occurence(dico) :

	reads = read_fasta(path)

	frequencies = []
	for idt, seq in reads.items() :
		frequencies.append(freq_aa(seq))
	#print(frequencies, type(frequencies))

	print("----------DF FREQ-----------")
	list_df = []
	for dicct in frequencies :
		dicct = pd.DataFrame([dicct])
		list_df.append(dicct)

	#print(list_df[1:100], type(list_df))
	
	for df in list_df :
		for col, val in df.iteritems() :
			if col not in list_of_aa :
				del df[col]

	for df in list_df :
		for col, val in df.iteritems() :
			if col not in list_of_aa :
				print(col, val)

	#print(list_df[1:100], type(list_df))

	return list_df


def Z_aa(Z_score) :
	#dico = listing(path)
	dico = read_fasta(path)
	#print("-------------DICO", dico)
	#print(Z_score)
	#print(Z_score.loc['H'])
	#print(Z_score.iloc[:, 1:])
	print("Nombre lignes : ", Z_score.shape[0])
	#pos = 0

	dico_Z = {}
	#new_seq = []
	#print(dico)
	for idt, seq in dico.items() :
		new_seq = []
		for aa in seq :
			new_seq.append(list(Z_score.loc[aa]))
			#print(new_seq)
		dico_Z[idt] = new_seq
	#print(new_seq[1:10])

		#print(new_seq)
	print("-----------vdfgegse-----------")	
	print(dico_Z)
	print("-----------vdfgegse-----------")

	return dico_Z




def Auto_cross_variance(Zscores) :

	print("--------ACC----------")

	#Rstudio()
	#ACCs = protr.acc(Zscores)
	#ACCs = robjects.r('acc(Zscores)')

	'''
	print(Zscores)
	ACC = r.acc(Zscores)
	print(ACC)	
	'''

	#robjects.r.source('acc.R') 

	lag_r = 4
	dico_Acc = {}
	for dicct in Zscores :
		for idt, score in dicct.items() :
			N = len(score)
			j = []
			k = []
			l = []
			Acc = []
			newj = []
			#for lag in range(0, N-lag_r, 4) : 
			#for lag in range(0, N, 4) : 
			for Z in score : 
				#print(Z)
				j.append(Z[0])
				k.append(Z[1])
				l.append(Z[2])
			#for i in range(0, len(j), 4) :
			#	newj.append(j[i])

			print("----j-----")
			print(j, len(j))
			print("----k----")
			print(k, len(k))
			print("----l----")
			print(l, len(l))
			#print(len(j))
			#print("-------newj", newj, len(newj))


			#print(sum(j))


			#mult = []
			#res_same_one = []
			#res_diff_one = []

			dico_jkl = {'j' : j, 'k' : k, 'l' : l}
			#print(dico_jkl)

			#mult = []
			#res_same_one = []
			mult = []
			for jkl in dico_jkl.values() :
				#print("---jkl----", jkl)

				##### Caler qql par un for i in range(lagr)
				#mult = []
				res_same_one = []
				res_diff_one = []
				#print(len(jkl))
				for i in range(len(jkl)-1) :
					#print(jkl[i], jkl[i+1])
					#break
					mult.append(jkl[i]*jkl[i+1])
				#print(mult)
				for m in mult :
					res_same_one.append((m**2+lag_r)/(N-lag_r))
				res_same = sum(res_same_one)
				Acc.append(res_same)
			dico_Acc[idt] = Acc

			length = len(j)
			print("---------length : ", length)
			mult2 = []
			mult3 = []
			mult4 = []
			res_diff_one = []
			for i in range(length) :
				mult2.append(j[i]*k[i])
				mult3.append(j[i]*l[i])
				mult4.append(k[i]*l[i])
			for m in mult2 : 
				res_diff_one.append((m+lag_r)/(N-lag_r))
			res_diff = sum(res_diff_one)
			Acc.append(res_diff)
			dico_Acc[idt] += Acc

			res_diff_one = []
			for m in mult3 : 
				res_diff_one.append((m+lag_r)/(N-lag_r))
			res_diff = sum(res_diff_one)
			Acc.append(res_diff)
			dico_Acc[idt] += Acc

			res_diff_one = []
			for m in mult4 : 
				res_diff_one.append((m+lag_r)/(N-lag_r))
			res_diff = sum(res_diff_one)
			Acc.append(res_diff)
			dico_Acc[idt] += Acc

	#print(dico_Acc)
	
	for i in dico_Acc.values() :
		print(i)
		print(len(i))
	
	return dico_Acc


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

	tsne = TSNE(n_components = 2, random_state = 0)
	#print(frequencies)

	fit_list = []
	print("--------")

	df = pd.DataFrame()
	matrix = []

	for mat in frequencies :
		matrix.append(mat.to_numpy()[0])
		if len(mat.columns) != 20 :
			print(len(mat.columns))
	
	for mat in matrix :
		mat = np.asarray(mat)
		#mat = np.reshape(mat, (1, 20))
	#print(matrix)
	
	matrix = np.array(matrix)
	
	for arr in matrix :
		if arr.shape != (20,) :
			print(arr.shape)
			arr = np.reshape(arr, 20,)
	
	print(matrix)
	print(matrix.shape)
	print(matrix[0].dtype)

	print("-------------OK----------")
	X_2d = tsne.fit_transform(matrix)
	df['x'] = X_2d[:,0]
	df['y'] = X_2d[:,1]
	print(df)

	df_data_tsne = df


	sns.scatterplot(x = 'x', y = 'y', data = df)
	plt.title("Tsne", fontsize = 15)
	
	plt.show()

	return df_data_tsne
	


def tsne_proteomes(path_to_proteom) :
	fich = listing(path_proteom)




if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)

	df_Score = Score_aa()

	# Ensemble de fichiers
	proportion = specific_occurence(path)
	#dico_score = Z_aa(df_Score)
	#ACCs = Auto_cross_variance(dico_score)
	#tsne = Tsne(proportion)
	#tsne_all_proteom = tsne_proteomes(path_proteom)
	#ACCs = Auto_cross_variance(dico_score)
	
	# Biopython
	#seq = bio()













