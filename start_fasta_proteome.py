import os
import sys
import glob
import pandas as pd
import numpy as np 
import seaborn as sns
from os.path import basename
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
path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/proteome_diatom.faa"
path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/"
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
	print(fichier)
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

	print("----NOT------")
	print(not_aa, len(not_aa))
	print(id_not_aa, len(id_not_aa))
	nb_delete = len(dico_delete.keys())
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


	dico2 = {}
	idt_list = []
	for idt, seq in dico.items() :
		idt_list.append(idt)
	for i in idt_list[:50] :
		dico2[i] = ""
		dico2[i] += dico[i]
				
	print("------------DICO2----------")
	print(dico2)
	print("------------FIN DICO2----------")
	#return(dico) 
	return(dico2)


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
	fich = glob.glob(path+'*.f'+'*')
	#print(fich)
	#fich.append(path+'*.faa')
	#fich.append(path+'*fna')
	#print(fich, type(fich))
	reads = []
	
	for fasta in fich :
		reads.append(read_fasta(str(fasta)))

	return reads


def specific_occurence(reads) :

	#reads = read_fasta(path) # ACTIVATE ONLY FOR 1 PROTEOM

	frequencies = []
	for idt, seq in reads.items() :
		frequencies.append(freq_aa(seq))

	print("----------DF FREQ-----------")
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


def Z_aa(Z_score) :
	#dico = listing(path)
	dico = read_fasta(path)
	#print("-------------DICO", dico, "---------FIN DICO")
	#print(Z_score)
	#print(Z_score.loc['H'])
	#print(Z_score.iloc[:, 1:])
	print("Nombre lignes : ", Z_score.shape[0])
	#pos = 0

	dico_Z = {}
	#new_seq = []
	for idt, seq in dico.items() :
		if '*' in seq[-1] :
			#dico[idt] = seq[:-1]
			seq = seq[:-1]
			#print(dico[idt])
		new_seq = []
		for aa in seq :
			#print(seq)
			new_seq.append(list(Z_score.loc[aa]))
			#print("ok2")
		dico_Z[idt] = new_seq
		#break

	#print("-----------vdfgegse-----------")
	#print(dico_Z)
	#print("-----------vdfgegse-----------")

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
	print(Zscores, type(Zscores))
	print(Zscores['>XP_002176126.1'])
	print(Zscores.keys())


	print("----------lalallalaa")
	for seq in Zscores.values() :
		print(seq, "\n")
	print("----------lalallalaa")

	r = 1
	dico_Acc = {}

	print("---TYPE ZSCORES---", type(Zscores.items()))
	
	for idt, seq in Zscores.items() :
		#print(idt,  seq)
		for aa in seq :
			#print(aa)
			#print(aa[2])
			pass
			#for Z in aa :
				#print(Z)
				#break
	
	#### SAME Z
	N = len(Zscores.keys())
	Z1 = []
	Z2 = []
	Z3 = []
	z_same = []
	# for k in range(1, r) :
	for idt, seq in Zscores.items() :
		for aa in seq :
			Z1.append(aa[0]) # tous les z1 de chaque seq
			Z2.append(aa[1]) # tous les z2 de chaque seq
			Z3.append(aa[2]) # tous les z3 de chaque seq
	#for Z in Z1 :
	for i in range(0, len(Z1), r) :
		#print(i, type(i))
		z_same.append(Z1[i])
	print(z_same, len(z_same))

	res = []
	for i in range(len(Z1)-1) :
		num = (Z1[i]*Z1[i+1])+r
		denum = N-r
		res.append(num/denum)
	z1z1 = sum(res)
	print("z1 -->", z1z1)

	z_same = []
	for i in range(0, len(Z2), r) :
		#print(i, type(i))
		z_same.append(Z2[i])
	print(z_same, len(z_same))

	res = []
	for i in range(len(Z2)-1) :
		num = (Z2[i]*Z2[i+1])+r
		denum = N-r
		res.append(num/denum)
	#print(sum(m1))
	z2z2 = sum(res)
	print("z2 -->", z2z2)

	z_same = []
	for i in range(0, len(Z3), r) :
		#print(i, type(i))
		z_same.append(Z3[i])
	print(z_same, len(z_same))

	res = []
	for i in range(len(Z3)-1) :
		num = (Z3[i]*Z3[i+1])+r
		denum = N-r
		res.append(num/denum)
	#print("RES z3 -->", res)
	z3z3 = sum(res)
	print("z3 -->", z3z3)


	# Z DIFF
	#print(len(Z1), "\n", len(Z2), "\n", len(Z3), "\n")
	
	mz1z2 = []
	mz1z3 = []
	mz2z1 = []
	mz2z3 = []
	mz3z1 = []
	mz3z2 = []
	
	for i in range (len(Z1)) : # ou Z2 ou Z3
		mz1z2.append(Z1[i]*Z2[i])
		mz1z3.append(Z1[i]*Z3[i])
		mz2z1.append(Z2[i]*Z1[i])
		mz2z3.append(Z2[i]*Z3[i])
		mz3z1.append(Z3[i]*Z1[i])
		mz3z2.append(Z3[i]*Z2[i])

	print(len(mz1z2))
	print(mz1z2[1:5])
	print(mz3z1[1:5])
	print(mz1z3[1:5])

	res = []
	for m in range(len(mz1z2)) : # no matter mz actually
		num = (mz1z2[m]+r)
		denum = N-r
		res.append(num/denum)
	z1z2 = sum(res)
	#print(res[0])

	res = []
	for m in range(len(mz1z3)) : # no matter mz actually
		num = (mz1z3[m]+r)
		denum = N-r
		res.append(num/denum)

	z1z3 = sum(res)
	print(z1z3)
	print(res[0])

	res = []
	for m in range(len(mz2z1)) : # no matter mz actually
		num = (mz2z1[m]+r)
		denum = N-r
		res.append(num/denum)
	z2z1 = sum(res)
	print(z2z1)

	res = []
	for m in range(len(mz2z3)) : # no matter mz actually
		num = (mz2z3[m]+r)
		denum = N-r
		res.append(num/denum)
	z2z3 = sum(res)
	print(z2z3)

	res = []
	for m in range(len(mz3z1)) : # no matter mz actually
		num = (mz3z1[m]+r)
		denum = N-r
		res.append(num/denum)
	z3z1 = sum(res)
	print(z3z1)

	res = []
	for m in range(len(mz3z2)) : # no matter mz actually
		num = (mz3z2[m]+r)
		denum = N-r
		res.append(num/denum)
	z3z2 = sum(res)
	print(z3z2)

	dico_Acc = {'z1z1_'+str(r) : z1z1, 'z1z2_'+str(r) : z1z2, 'z1z3_'+str(r) : z1z3, \
		'z2z1_'+str(r) : z2z1, 'z2z2_'+str(r) : z2z2, 'z2z3_'+str(r) : z2z3, \
		'z3z1_'+str(r) : z3z1, 'z3z2_'+str(r) : z3z2, 'z3z3_'+str(r) : z3z3}

	print(dico_Acc)






	'''
	for dicct in Zscores :
		print(dicct.items())
		break
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
	'''
	return dico_Acc


def Tsne(frequencies) :

	tsne = TSNE(n_components = 2, random_state = 0)

	print("--------")

	df = pd.DataFrame()
	matrix = []


	for mat in frequencies :
		print(mat)
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


	print("-------------OK----------")
	X_2d = tsne.fit_transform(matrix)
	df['x'] = X_2d[:,0]
	df['y'] = X_2d[:,1]
	#print(df)

	df_data_tsne = df
	print(df)

	#sns.scatterplot(x = 'x', y = 'y', data = df). # For 1 proteom
	#plt.title("Tsne", fontsize = 15)
	
	#plt.show()

	return df_data_tsne
	


def tsne_proteomes(path_to_proteom) :
	reads = listing(path_to_proteom)
	occ = []
	
	for dico in reads :
		occ.append(specific_occurence(dico))
	print("OCC", occ, type(occ), len(occ))
	print("--------------------------------------------------")
	
	tsne = []
	list_df = []

	for proteom in occ :
		tsne.append(Tsne(proteom))
	

	print("--------------TSNE--------------")
	print(tsne, len(tsne), type(tsne))

	label = proteom_name(path_proteom)

	
	for data in tsne :
		sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	plt.title("Tsne of frequencies", fontsize = 15)
	
	plt.show()


def proteom_name(path_to_proteom) :
	
	fich = glob.glob(path_to_proteom+'*.f'+'*')

	label = []
	for f in fich :
		label.append(basename(f))

	return label




if __name__ == '__main__' :
	# 1 fichier
	#reading, sequence = read_fasta("/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Exemple_fasta/F4HXU3.fasta")
	#dico_number, frequency = freq_aa(sequence)

	df_Score = Score_aa()

	# Ensemble de fichiers
	#proportion = specific_occurence(path)
	dico_score = Z_aa(df_Score)
	ACCs = Auto_cross_variance(dico_score)
	#tsne = Tsne(proportion)
	#ACCs = Auto_cross_variance(dico_score)
	
	# ALL proteom
	#tsne_all_proteom = tsne_proteomes(path_proteom)












