"""Code that allows to perform Tsne 
   on amino acid frequency and ACC (based on Zscale)
   for a set of proteoms (fasta format)

------------------------------------------------------------------
Rebecca GOULANCOURT
M2 BIOLOGIE - INFORMATIQUE
Université de Paris 2021 - 2022
Stage M2 - supervisor : Ingrid Lafontaine
------------------------------------------------------------------

This code calls another script in R : acc.R.

"""

# Import of the necessary modules
import os
import sys
import glob
import csv
import subprocess
import pandas as pd
import numpy as np 
import seaborn as sns
from os.path import basename
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt

# Variables 
path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/proteome_diatom.faa"
path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/"
path_output = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/"
list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']


def Score_aa() :
	""" Initialization of the Zscales dataframefor each amino acid

	Parameters
	----------
	None

	Returns
	-------
	df_Zscore : dataframe
		Dataframe that contains the (Z1, Z2, Z3) of the 20 amino acids

	"""
	Z1 = [-2.49, 2.18, 0.07, -4.19, 1.96, -4.44, -1.22, 2.84, 2.23, -2.69, 2.88, 3.08, -4.92, 3.64, 0.71, 0.92, \
			3.22, -4.75, -1.39, 2.41]
	Z2 = [-0.27, 0.53, -1.73, -1.03, -1.63, -1.68, 0.88, 1.41, -5.36, -2.53, 2.52, 0.39, 1.30, 1.13, -0.97, -2.09, \
			1.45, 3.65, 2.32, 1.74]
	Z3 = [-0.41, -1.14, 0.09, -0.98, 0.57, -1.03, 2.23, -3.14, 0.30, -1.29, -3.44, -0.07, 0.45, 2.36, 4.13, -1.40, \
			0.84, 0.85, 0.01, 1.11]
	df_Zscore = pd.DataFrame({'AA' : list_of_aa, 'Z1' : Z1, 'Z2' : Z2, 'Z3' : Z3})
	df_Zscore.set_index('AA', inplace = True)

	return df_Zscore


######### These functions below are made for 1 proteom/file or sequence and will be recalled later

def read_fasta(fichier) :
	""" Reading and parsing of a proteom file
		Also writes in a .txt file the results of sequence filtering

	Parameters
	----------
	fichier : str
		Name of the proteom file

	Returns
	-------
	dico : dictionnary
		Dictionnary where keys represents the sequence id and
		values are the corresponding sequence

	"""
	print("File : ", fichier)
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

	dico_delete = {}	
	not_aa = []
	id_not_aa = []

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
			dico_delete[ident] += dico[ident]

	nb_delete = len(id_not_aa)
	print("Number of deleted sequences : ", nb_delete)


	with open("Deleted_seq_"+str(basename(fichier)), "w") as filout :
		for idt, seq in dico_delete.items() :
			filout.write(idt+"\n")
	
	print("Before :", len(dico.keys()), len(dico.values()))

	dico_nb = {}
	for idt, seq in dico.items() :
		dico_nb[idt] = parsing_seq(seq)

	dico_total = {}

	for dic in dico_nb.values() :
		for cle, val in dic.items() :
			if cle not in dico_total.keys() :
				dico_total[cle] = 0
			else :
				dico_total[cle] += val

	for idt in id_not_aa :
		del dico[idt]

	print("After :", len(dico.keys()), len(dico.values()))
	#print("On the proteom we count : ", dico_total)
	'''
	with open("New_proteome_"+str(basename(fichier)), "w") as filout :
		for idt, seq in dico.items() :
			filout.write(idt+"\n")
	'''

	# If you want to try o a limited number of sequence for each proteom
	dico2 = {}
	idt_list = []
	for idt, seq in dico.items() :
		idt_list.append(idt)
	for i in idt_list[:50] :
		dico2[i] = ""
		dico2[i] += dico[i]
				
	#return(dico) 
	return(dico2)


def parsing_seq(seq) :
	""" Function that reads a sequence to verify its nature 
		i.e presence of * or a non amino acid

	Parameters
	----------
	seq : str
		Sequence to analyze

	Returns
	-------
	dico : dictionnary
		Dictionnary where keys represents the nature of the sequence

	"""
	dico = {'stop codon' : 0, 'pseudo-gene' : 0, 'other' : 0}
	
	if seq[-1] == '*' :
		dico['stop codon'] += 1
		for aa in seq[:-1] :
			if aa == '*' :
				dico['pseudo-gene'] += 1
			elif aa not in list_of_aa :
				dico['other'] += 1
			elif aa in list_of_aa :
				continue
			else :
				break

	elif seq[-1] != '*' :
		for aa in seq :
			if aa == '*' :
				dico['pseudo-gene'] += 1
			elif aa not in list_of_aa :
				dico['other'] += 1
			elif aa in list_of_aa :
				continue
			elif (aa == '*') or (aa not in list_of_aa) :
				break


	return dico


def freq_aa(sequence) : 
	""" Calculates for a given sequence its frequency of amino acid

	Parameters
	----------
	sequence : str
		Sequence to analyze

	Returns
	-------
	freq_dico : dictionnary
		Dictionnary of frequencies, where each key is a amino acid
		and the dictionnary value its frequency in the sequence

	"""
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


	return freq_dico


######### These funtions below are written in order to read many proteoms

def listing(path) :
	""" Collects proteoms files in a directory and applies the read_fasta
		function on each proteoms

	Parameters
	----------
	path : str
		path to the proteoms (should finish by "/")

	Returns
	-------
	reads : list
		List of dictionnaries for all proteoms that contains sequence id and the sequence

	"""
	fich = glob.glob(path+'*.f'+'*')
	reads = []
	
	for fasta in fich :
		reads.append(read_fasta(str(fasta)))

	return reads


def specific_occurence(reads) :
	""" Calculates for a given proteom the amino acid frequency of each sequences
		Also take into account the nature of the sequence
		--> frequencies with unknown amino acid or '*' are deleted

	Parameters
	----------
	reads : dictionnary
		Dictionnary of id and sequences

	Returns
	-------
	list_df : list
		List of dataframes that contains amino acid frequency 
		--> each dataframe is a sequence

	"""

	#reads = read_fasta(path) # ACTIVATE IF YOU WANT TO READ ONLY 1 SPCIFIC PROTEOM

	frequencies = []
	for idt, seq in reads.items() :
		frequencies.append(freq_aa(seq))

	list_df = []
	for dicct in frequencies :
		dicct = pd.DataFrame([dicct])
		list_df.append(dicct)
	
	for df in list_df :
		for col, val in df.iteritems() :
			if col not in list_of_aa :
				del df[col]

	for df in list_df :
		for col, val in df.iteritems() :
			if col not in list_of_aa :
				print(col, val)

	return list_df


def Tsne(frequencies) :
	""" Calculates the Tsne on a proteom for its amino acid frequency

	Parameters
	----------
	frequencies : list
		List of frequecies as a dataframe for a proteom and its sequences

	Returns
	-------
	df_data_tsne : dataframe
		Dataframe that contains the Tsne results

	"""
	tsne = TSNE(n_components = 2, random_state = 0)

	df = pd.DataFrame()
	matrix = []


	for mat in frequencies :
		matrix.append(mat.to_numpy()[0]) 
		if len(mat.columns) != 20 :
			print(len(mat.columns))


	for mat in matrix :
		mat = np.asarray(mat)
	
	matrix = np.array(matrix)
	
	for arr in matrix :
		if arr.shape != (20,) :
			print(arr.shape)
			arr = np.reshape(arr, 20,)

	X_2d = tsne.fit_transform(matrix)
	df['x'] = X_2d[:,0]
	df['y'] = X_2d[:,1]

	df_data_tsne = df


	return df_data_tsne
	


def tsne_proteomes(path_to_proteom) :
	""" Perform a tsne for all the selected proteoms

	Parameters
	----------
	path_to_proteom : str
		Path to the directory that contains all the proteoms

	Returns
	-------
	None

	Plot
	----
		A scatterplot for the tsne on frequencies for all the proteoms

	"""
	reads = listing(path_to_proteom)
	occ = []
	
	for dico in reads :
		occ.append(specific_occurence(dico))
	
	tsne = []
	list_df = []

	for proteom in occ :
		tsne.append(Tsne(proteom))
	

	print("--------------TSNE ON FREQUENCY PERFORMING--------------")

	label = proteom_name(path_proteom)

	
	for data in tsne :
		sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	plt.title("Tsne of frequencies", fontsize = 15)
	
	plt.show()


def proteom_name(path_to_proteom) :
	""" Save the proteom name only (not the entire path) 

	Parameters
	----------
	path_to_proteom : str
		Path to the directory that contains all the proteoms

	Returns
	-------
	label : str
		Each proteoms name 

	"""
	fich = glob.glob(path_to_proteom+'*.f'+'*')

	label = []
	for f in fich :
		label.append(basename(f))

	return label

	'''
def Acc_Tsne(path_to_output) :
	""" Calculates via R the acc of each sequence of each proteoms and perform a tsne

	Parameters
	----------
	None

	Returns
	-------
	None

	Plot
	----
		A tsne of ACC and put the pdf image in the working directory

	"""	
	

	#os.system('Rscript --vanilla acc.R') # à la place lire fichier des ACC issu de R
	#print("Script acc.R lancé")

	
	fich = glob.glob(path_to_output+'Acc_output_*')
	names = proteom_name(path_output+'Acc_output_*')
	print("FICH", fich, '\n', "NAMES", names)

	list_df = []
	#essai = fich[0]
	for f in fich :
		df = pd.read_csv(f, sep = "\t", header = None)
		#df = df.head(100)
		list_df.append(df)
	#print(list_df)

	matrix = []
	for df in list_df :
		if len(df.columns) != 36 :
			print(df, len(df.columns))
		#print(df)
		#print(type(df))
		matrix.append(df.to_numpy()[0])
		#matrix.append(df.to_numpy())
	
	print("----------ICI----------")
	#for m in matrix :
	#	print(type(m))
	
	#print(data)
	#for i in matrix :
	#	print(i, type(i), i.shape)

	for mat in matrix :
		mat = np.asarray(mat)
	
	#matrix = np.array(matrix)
	#matrix = np.array(matrix, ndmin = 2)[0]
	matrix = np.array(matrix, ndmin = 2)

	print(matrix)
	for arr in matrix :
		#print(arr)
		print(arr.shape)
		print(len(arr))
		#if arr.shape != (100, 36) :
		#	print("--------DIFF")
		#	print(arr.shape)
		#	arr = np.reshape(arr, 100, 36)
		#print(arr)
		if arr.ndim == 3 :
			print(arr)
			print(arr.dim)
			print("----")
		else :
			print("None")

		#arr.reshape(len(arr), 36)
	

	tsne = TSNE(n_components = 2, random_state = 0)
	#data = []
	#for m in matrix :
	#	break
		#X_2d = tsne.fit_transform(matrix)
	#	data.append(tsne.fit_transform(m))
	print("OK -----------")


	#list_data = []
	#for proteom in label :
	X_2d = tsne.fit_transform(matrix)
	print("ALORSPEUTETRE")

	data = pd.DataFrame()
	data['x'] = X_2d[:,0]
	data['y'] = X_2d[:,1]
		#list_data.append(data)

	
	sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	plt.title("Tsne of frequencies", fontsize = 15)
	
	plt.show()
	'''

	'''
	print("-------------TSNE ON ACC PERFORMING------------")

	
	for data in list_data :
		sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	plt.title("Tsne of frequencies", fontsize = 15)
	
	plt.show()
	
	return data

	'''
def tsne_data(matrix) : 

	tsne = TSNE(n_components = 2, random_state = 0)
	#data = []
	#for m in matrix :
	#	break
		#X_2d = tsne.fit_transform(matrix)
	#	data.append(tsne.fit_transform(m))
	print("OK -----------")


	#list_data = []
	#for proteom in label :
	X_2d = tsne.fit_transform(matrix)
	print("ALORSPEUTETRE")

	data = pd.DataFrame()
	data['x'] = X_2d[:,0]
	data['y'] = X_2d[:,1]

	return data

def ACC_tsne_plot() :


	fich = glob.glob(path_output+'Acc_output_*')
	names = proteom_name(path_output+'Acc_output_*')
	print("FICH", fich, '\n', "NAMES", names)

	list_df = []
	#essai = fich[0]
	for f in fich :
		df = pd.read_csv(f, sep = "\t", header = None)
		df = df.head(10000)
		list_df.append(df)
	#print(list_df)


	matrix = []
	for df in list_df :
		if len(df.columns) != 36 :
			print(df, len(df.columns))
		#print(df)
		#print(type(df))
		matrix.append(df.to_numpy()[0])
		#matrix.append(df.to_numpy())
	
	print("----------ICI----------")
	#for m in matrix :
	#	print(type(m))
	
	#print(data)
	#for i in matrix :
	#	print(i, type(i), i.shape)

	for mat in matrix :
		mat = np.asarray(mat)
	
	#matrix = np.array(matrix)
	#matrix = np.array(matrix, ndmin = 2)[0]
	matrix = np.array(matrix, ndmin = 2)

	print(matrix)
	for arr in matrix :
		#print(arr)
		print(arr.shape)
		print(len(arr))
		#if arr.shape != (100, 36) :
		#	print("--------DIFF")
		#	print(arr.shape)
		#	arr = np.reshape(arr, 100, 36)
		#print(arr)
		if arr.ndim == 3 :
			print(arr)
			print(arr.dim)
			print("----")
		else :
			print("None")

		#arr.reshape(len(arr), 36)

	tsne = []
	#for proteom in matrix :
	for proteom in list_df :
		tsne.append(tsne_data(proteom))

	label = []
	for n in names :
		label.append(n)

	for data in tsne :
		sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	plt.title("Tsne of frequencies", fontsize = 15)
	
	plt.show()



if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)
	#proportion = specific_occurence(path)
	#tsne = Tsne(proportion)
	
	# ALL proteoms
	#tsne_all_proteom = tsne_proteomes(path_proteom)
	ACC_tsne = ACC_tsne_plot()












