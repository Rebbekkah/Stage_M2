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
path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/"
list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']



def Rstudio() :
	pass
	#robjects.r.source('acc.R') 


	'''
	#import rpy2
	#from rpy2 import *
	import rpy2.robjects as robjects
	import rpy2.robject.packages as rpackages
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


	for aa in list_of_aa :
		if aa not in freq_dico :
			freq_dico[aa] = 0


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
	print("-------------DICO", dico)
	print(Z_score)
	#print(Z_score.loc['H'])
	#print(Z_score.iloc[:, 1:])
	print("Nombre lignes : ", Z_score.shape[0])
	#pos = 0

	for dicct in dico :
		new_seq = []
		#print(dicct)
		for idt, seq in dicct.items() :
			for aa in seq :
				new_seq.append(list(Z_score.loc[aa]))
			dicct[idt] = new_seq

		#print(new_seq)
	print("-----------vdfgegse-----------")	
	print(dico)
	print("-----------vdfgegse-----------")

	return dico




def Auto_cross_variance(Zscores) :
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
	ACCz_zj_lag4 = []
	ACC_zjzk_lag4 = []

	#print(Zscores)


	'''
	for dicct in Zscores :
		#print(dicct)
		#break
		Acc = []
		for idt, score in dicct.items() :
			N = len(score)
			#print(len(score))
			#for i in range(len(score)) : 
			for Z in score : 
				print("------Z-----------", Z)
				for s in Z :
					z = s
					for i in range(0, 3) :
						#print("------I-----------", i)
						res_same = sum((z[i]**2 + lag_r)/(N - lag_r))

	'''
	dico_Acc = {}
	#newj = []
	#Acc = []
	for dicct in Zscores :
		#print(dicct)
		#break
		#Acc = []
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
				##### Caler qql par un for i in range(lagr)
				#mult = []
				res_same_one = []
				res_diff_one = []
				#print(len(jkl))
				for i in range(len(jkl)-1) :
					print(jkl[i], jkl[i+1])
					#break
					mult.append(jkl[i]*jkl[i+1])
				#print(mult)
				for m in mult :
					res_same_one.append((m**2+lag_r)/(N-lag_r))
				res_same = sum(res_same_one)
				Acc.append(res_same)
				dico_Acc[idt] = Acc

			length = len(j)

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
			


	print(dico_Acc)
			
		
'''
			mult = []
			res_same_one = []
			for i in range(len(j)-1) :
				mult.append(j[i]*j[i+1])
			#print(len(mult))
			for m in mult :
				#print(type(m))
				#print(m)
				res_same_one.append((m**2+lag_r)/N-lag_r)
			#print("----res_same1----", res_same_one, len(res_same_one))
			res_same = sum(res_same_one)
			Acc.append(res_same)
			dico_Acc[idt] = Acc
		
			#Acc.append(sum((j**2+lag_r)/N-lag_r))
	print(Acc)
	print(dico_Acc)
			#break
'''
	


'''
def formula_acc(list) :
	dico_Acc = {}
	mult = []
	res_same_one = []
	for i in range(len(j)-1) :
		mult.append(j[i]*j[i+1])
	#print(len(mult))
	for m in mult :
		#print(type(m))
		#print(m)
		res_same_one.append((m**2+lag_r)/N-lag_r)
	#print("----res_same1----", res_same_one, len(res_same_one))
	res_same = sum(res_same_one)
	Acc.append(res_same)
	dico_Acc[idt] = Acc
'''


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
	print(frequencies)

	fit_list = []
	print("--------")
	#for i in range(len(frequencies)-1) :
		#print(i)

	'''
	for mat in frequencies :
		print(mat)
		#break
		arr = mat.to_numpy()
		print(arr)
		print(type(arr))

		arr = np.reshape(arr, (-1, 2))
		print(arr)

		X_2d = tsne.fit_transform(arr)
		fit_list.append(X_2d)
		print("---------")
		print(X_2d)
		#plt.scatter(X_2d, )
		#plt.plot(X_2d, linestyle = 'none', marker = 'o', color = 'red', markersize = 12)
		#plt.show()
		#break
	'''

	df = pd.DataFrame()
	matrix = []

	'''
	for mat in frequencies :
		arr = mat.to_numpy()
		print(arr)
		print(type(arr))

		arr = np.reshape(arr, (-1, 2))
		print(arr)

		X_2d = tsne.fit_transform(arr)
		df['x'] = X_2d[:,0]
		df['y'] = X_2d[:,1]
		print(df)

		sns.scatterplot(x = 'x', y = 'y', data = df)
		plt.show()
	'''

	for mat in frequencies :
		matrix.append(mat.to_numpy())

	matrix = np.asarray(matrix)
	matrix = np.reshape(matrix, (3, -1))
	print(matrix.ndim)
	print(type(matrix))
	X_2d = tsne.fit_transform(matrix)
	df['x'] = X_2d[:,0]
	df['y'] = X_2d[:,1]
	print(df)


	sns.scatterplot(x = 'x', y = 'y', data = df)
	plt.title("Tsne", fontsize = 15)

	'''
	for i, label in enumerate(annotations) :
		plt.annotate()
	'''
	
	plt.show()


	


if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)

	df_Score = Score_aa()

	# Ensemble de fichiers
	proportion = specific_occurence(path)
	dico_score = Z_aa(df_Score)
	ACCs = Auto_cross_variance(dico_score)
	#tsne = Tsne(proportion)

	
	# Biopython
	#seq = bio()













