"""Code that parses and analyze the different outputs of softwares
   and returns a dataframe for each proteins and their results

Softwares : 
	- TMHMM
	- ARD2
	- DEEPLOC
	- LOCALIZER
	- RADAR
	- WOLFPSORT

------------------------------------------------------------------
Rebecca GOULANCOURT
M2 BIOLOGIE - INFORMATIQUE
Université de Paris 2021 - 2022
Stage M2 - supervisor : Ingrid Lafontaine & Céline Cattelin
------------------------------------------------------------------

"""

# Modules
import pandas as pd
import numpy as np
#print(np.version.version)
import sys
import os
import statistics
from os.path import basename
import glob
from operator import itemgetter
import seaborn as sns
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
import scipy.stats.distributions as dist
#import umap
#import umap.plot
#from sklearn.decomposition import PCA
from scipy.stats import shapiro
from scipy.stats import f
import scipy.stats as stats



def reading(fichier) :
	''' Function that reads the output of the TMHMM software and decides
		wether the protein id is to keep for the rest of the analysis

	Parameters
	----------
	fichier : str
		Output file

	Writes
	------
	output_Analyzeseq_ToKeep.txt : txt file
		File where each line contain an protein id to keep

	Returns
	-------
	id_to_keep : list
		List of protein identifiers that doesn't have transmembrane domain
		(after the 68 amino acid)

	'''
	print("File :", basename(fichier))
	with open(fichier, "r") as filin :
		id_list_ok = []
		idt_all = []
		dico = {}
		dico_h = {}
		dico_ok = {}
		id_TB = []
		for line in filin :
			if line.startswith('#') : 
				i = line.split(' ')[1]
				if i not in idt_all :
					idt_all.append(i)
			else :	
				if line.split('\t')[2] == 'TMhelix' :
					id_TB.append(line.split('\t')[0])
					dico[line.split('\t')[0]] = ""
					dico[line.split('\t')[0]] = line.split('\t')[3].strip()
				else : 
					idt = line.split('\t')[0]
					id_list_ok.append(idt)

	for cle, val in dico.items() :
		start = int(val.split(' ')[0])
		end = int(val.split(' ')[-1])

		if start > 68 or end > 68 :
			dico_h[cle] = {}
			dico_h[cle] = {'start' : start, 'end' : end}
		elif end < 68 :
			if cle not in id_list_ok :
				id_list_ok.append(cle)
			dico_ok = {}
			dico_ok[cle] = {'start' : start, 'end' : end}

	id_to_keep = []
	for elem in id_list_ok :
		if elem not in id_to_keep and elem not in dico_h.keys():
			id_to_keep.append(elem)

	#print("Total prot : ", idt_all, len(idt_all))
	#print("PROT TO KEEP : ", id_to_keep, len(id_to_keep))

	id_to_delete = []	
	for idt in idt_all :
		if idt in dico_h.keys() :
			id_to_delete.append(idt)
	#print("TO DELETE : ", id_to_delete)
	
	#os.chdir(path_output)
	#with open("output_Analyzeseq_ToKeep_"+basename(fichier)+".txt", "w") as filout :
	#	for idt in id_to_keep :
	#		filout.write(idt+"\n")


	return id_to_keep


def listing(path, pattern) :
	''' Funtion that perform the reading function on many files

	Parameters
	----------
	path : str
		Path to where the tmhmm outputs are located

	pattern : str
		Pattern to recognize the tmhmm files 

	Returns
	-------
	prot : list
		List of the id to keep for each file, 1 element of the list = 1 file

	'''

	fich = glob.glob(path+pattern)
	print(fich)
	prot = []
	
	for f in fich :
		print("----------", basename(f), "----------")
		prot.append(reading(f))

	return prot



def proteome_maker(ToKeep, path_proteom, pattern) :
	''' Production of the new proteom from the id that were selected thru 
		the reading funtion ('meta proteom' that concatenate all ids from the files)

	Parameters
	----------
	ToKeep : list
		List of protein identifiers that we kept

	path_proteom : str
		Path to the initial proteoms that contains id and their sequences

	pattern : str
		pattern to recognize the fasta proteome files

	Returns
	-------
	dico2 : dict
		Dictionnary of the metaproteom (key = id, value = sequence)

	'''
	print("-----------------proteom maker-----------------")
	proteom = glob.glob(path_proteom+pattern)
	print(proteom, len(proteom))

	dico = {}
	dico2 = {}
	for p in proteom :
		print(basename(p))
		with open(p, "r") as filin :
			for line in filin :
				if line.startswith('>') :
					idt = line.split()[0]
					dico[idt] = ""
				else : 
					dico[idt] += line.strip()
	for proteom in ToKeep :
		for ident in proteom :
			ident = ">"+ident
			if ident in dico.keys() :
				dico2[ident] = ""
				dico2[ident] = dico[ident]
				

	idt_del = []
	for idt, seq in dico2.items() :
		if seq[-1] == '*' :
			new_seq = seq[:-1]
			dico2[idt] = new_seq
	
	for seq in dico2.values() :
		for aa in seq :
			if aa not in list_of_aa :
				idt_del.append(idt)

	for idt in idt_del :
		if idt in dico2.keys() :
			#print(idt)
			del dico2[idt]
	print("LEN DICO2------", len(dico2))


	with open("New_Proteom_All.txt", "w") as filout :
		for idt, seq in dico2.items() :
			if len(seq) > 120 :
				filout.write(idt+"\n"+seq+"\n")
	
	#return dico2


def sep(path_proteom, pattern1, pattern2, pattern3) :
	''' Function that separates the meta proteom in the positive and
		negative proteom that we had initially

	Parameters
	----------
	path_proteom : str
		Path to where the proteoms are located

	pattern1 : str
		pattern to recognize fasta proteom files

	pattern2 : str
		pattern to recognize the TMHMM files

	pattern3 : str
		pattern to locate the output results

	Returns
	-------
	None

	Writes
	------
	New_Proteom.txt : txt file
		File of the new proteoms 

	'''
	print("-----------------separateur-----------------")
	proteom = glob.glob(path_proteom+pattern1)
	proteom.sort()
	print(proteom, len(proteom))
	file = glob.glob(path_proteom+pattern2)
	file.sort()
	print(file, len(file))

	idt_ok = []
	for f in file :
		idt_ok.append(reading(f))

	os.chdir(path_output+pattern3)
	idt_all = []
	dico2 = {}
	dico_all = {}
	i = 0
	for p in proteom :
		dico2 = {}
		dico_all = {}
		print(basename(p))
		with open(p, "r") as filin :
			for line in filin :
				if line.startswith('>') :
					elem = line.split()[0]
					dico_all[elem] = ""
				else :
					dico_all[elem] += line.strip()
			print(len(dico_all.keys()))
			#print(dico_all.keys())
		for idt in idt_ok[i] :
			idt = ">"+idt
			#print(idt)
			if idt in dico_all.keys() :
				#print("oui")
				dico2[idt] = ""
				dico2[idt] += dico_all[idt]
		print("LEN DICO2 KEYS", len(dico2.keys()))
		#print("LEN DICO2 VAL", len(dico2.values()))
		with open("New_Proteom_"+basename(p)+".txt", "w") as filout :
			for idt, seq in dico2.items() :
				if len(seq) > 120 :
					filout.write(idt+"\n"+seq+"\n")
		print("--------------------")
		i += 1


def ard2(file, pattern) :
	''' Read and parses ard2 output that we had on our new proteoms
		(proteom of proteins without transmembrane domains)
	
	Parameters
	----------
	file : str
		Output ard2 file to perform the analysis

	Returns
	-------
	dico_f : dict
		Dictionnary of the amino acid linker and its probability ( > 0.10)

	'''

	print(basename(file))

	dico = {}
	with open(file, 'r') as filin :
		for line in filin :
			if line.startswith('>') :
				idt = line.split()[0]
				dico[idt] = ""
				dico_linker = {}
				pos = 0
				window = []
			else :
				pos += 1
				elem = line.split("\t")[0]
				elem = elem.split()[0]
				aa = elem[0]
				proba = float(elem[1:5])
				window.append([aa, proba, pos])

			if len(window) == 6 :
				l = []
				for elem in window :
					if elem[1] > 0.1 :
						l.append(elem[1])

				if not not l :
					p = max(l)

					for elem in window :
						if p == elem[1] :
							if elem[0] in dico_linker.keys() :
								dico_linker[elem[0]].extend([elem[1:]])
							else :
								dico_linker[elem[0]] = [elem[1:]]

				window = []
			dico[idt] = dico_linker

	#print(dico)
	proteoms = glob.glob(path_output+pattern)
	#print(proteoms)

	dico_Arabi = {}
	dico_Chlamy = {}
	dico_pos = {}
	dico_neg = {}
	dico_interm = {}
	for p in proteoms :
		#print(basename(p))
		#print(basename(file))
		#print("----------")
		base = basename(file.split('STDOUT_')[1])
		if base in basename(p) : 
			print(base)
			print(basename(p))
			print(basename(file))
			with open(p, "r") as filin :
				for line in filin :
					if line.startswith('>') :
						idt = line.strip()
						dico_interm[idt] = ""
					else :
						dico_interm[idt] += line.strip()
		
	#print(dico_interm)
	dico_f = {}
	link = []
	
	#elif basename(file) == 'STDOUT_' :
	if 'STDOUT_' in basename(file) :
		for idt in dico_interm.keys() :
			dico_f[idt] = {}

		for linker in dico.values() :
			link.append(linker)

		#print(link)
		for index, key in enumerate(dico_f) :
			#print(dico_f)
			#print(key)
			#print(len(link), len(dico_f.keys()))
			dico_f[key] = link[index]


	'''
	else : 
		#if basename(file) == 'STDOUT' :
		for idt in dico_else.keys() :
			dico_f[idt] = {}

		for linker in dico.values() :
			link.append(linker)
		print(len(link))
		print(len(dico.keys()))
		print(len(dico_f.keys()))
		print(len(dico_else.keys()))
		#print(dico_else.keys())
	
		for index, key in enumerate(dico_f) : 
			dico_f[key] = link[index]
		
	'''

	return dico_f


def wolfpsort(file) :
	''' Read and parses wolfpsort output that we had on our new proteoms
		(proteom of proteins without transmembrane domains)
	
	Parameters
	----------
	file : str
		Output wolfpsort file to perform the analysis

	Returns
	-------
	dico : dict
		Dictionnary of the adressing for each protein

	'''
	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			l_mc = []
			l_tot = []
			if not line.startswith('#') :
				elem = line.split(" ")
				elem[-1] = elem[-1].strip()
				idt = elem[0]
				#print(elem, len(elem), "\n", idt)
				for i in range(len(elem)) :
					if elem[i][-1] == ',' :
						elem[i] = elem[i][:-1]
				#print(elem, len(elem))
				if len(elem)%2 == type(int(1)) :
					print("____________ATTENTION____________")
					print(elem, len(elem))
				

				for i in range(2, len(elem), 2) :
					l_tot.append(float(elem[i]))

				for i in range(len(elem)-1) :
					if elem[i] == 'chlo' or elem[i] == 'mito' :
						#print(elem[i], elem[i+1])
						l_mc.append(float(elem[i+1]))
				#print(l_tot)
				#print(l_mc)
				score_tot = sum(l_tot)
				score_mc = sum(l_mc)
				ratio = score_mc/score_tot
				#print(score_tot, score_mc, ratio)
				#print("---------")

				dico[idt] = 0
				dico[idt] += ratio


	return dico


def targetp2(file) :
	''' Read and parses targetp2 output that we had on our new proteoms
		(proteom of proteins without transmembrane domains)
	
	Parameters
	----------
	file : str
		Output targetp2 file to perform the analysis

	Returns
	-------
	dico : dict
		Dictionnary of the adressing for each protein

	'''
	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if not line.startswith('#') :
				idt = line.split()[0]
				dico[idt] = ""
				dico[idt] += line.split()[1]

	return dico


def deeploc(file) :

	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if not line.startswith('ID') :
				idt = line.split()[0]
				mito = line.split()[6]
				plastid = line.split()[9]
				dico[idt] = {}
				dico[idt] = {'mito' : mito, 'plastid' : plastid}
				#dico[idt] = [mito, plastid]
			#else :
			#	print(line.split())

	return dico



def localizer(file) : 
	''' Read and parses localizer output that we had on our new proteoms
		(proteom of proteins without transmembrane domains)
	
	Parameters
	----------
	file : str
		Output localizer file to perform the analysis

	Returns
	-------
	dico : dict
		Dictionnary of the adressing for each protein

	'''
	dico = {}
	non = ['Over', '#', 'Identifier', '-']
	with open(file, "r") as filin :
		for line in filin :
			first = line.split(" ")[0]
			if first not in non :
				one = line.split("\t")[0]
				if one[0] not in non :
					idt = line.split("\t")[0]
					idt = idt.split()[0]
					elem = line.split("\t")[1]
					if elem.split()[0] != '-' :
						dico[idt] = ""
						dico[idt] += elem.split()[0]
					#else : 
					#	dico[idt] += 'None'
				else :
					break

	return dico		


def radar(file, path_proteom_all) :
	''' Read and parses radar output that we had on our new proteoms
		(proteom of proteins without transmembrane domains)
	
	Parameters
	----------
	file : str
		Output radar file to perform the analysis

	Returns
	-------
	dico : dict
		Dictionnary of the postion, length and proportion of amino acid
		of a found repetition in the sequence for each protein

	'''

	print(basename(file))

	idt_list = []
	dico = {}
	with open(file, "r", encoding = 'ascii', errors = 'ignore') as filin :
		for line in filin :
			if line.startswith('>') or '>' in line :
				repet = []
				l = []
				aa_prop = []
				idt_list.append(line.strip())
				idt = line.strip()
				dico[idt] = {'l_seq', 'l_rep', 'pos', 'rep_prop', 'aa_prop', 'seq'}
			if '- ' in line :
				ligne = line.strip()
				ligne = ligne.split()
				debut = ligne[0].split('-')
				pos = [debut[0], ligne[1]]
				repet.append(pos)


			dico[idt] = {'pos' : repet, 'l_rep' : l, 'aa_prop' : aa_prop, \
			'seq' : ""}

		k = 0
		for idt, dic in dico.items() :
			to_remove = []
			to_change = []
			rep = []
			to_modif = []
			nb = 0
			for repet in dic['pos'] :
				rep.append(repet)
			rep2 = []
			for lis in rep :
				p = []
				for ps in lis :
					p.append(int(ps))
				rep2.append(p)

			rep2 = sorted(rep2, key = itemgetter(0))
			dic['pos'] = sorted(rep2, key = itemgetter(0))
			
			for i in range(len(rep2)-1) : 
				if rep2[i+1][0] > rep2[i][0] and rep2[i+1][1] < rep2[i][1] :
					nb += 1
					#print(rep2, len(rep2))
					#print(rep2[i], rep2[i+1])
					to_remove.append(rep2[i+1])
			
				if rep2[i][0] > rep2[i+1][0] and rep2[i][1] < rep2[i+1][0] :
					nb += 1
					#print(rep2)
					#print(rep2[i], rep2[i+1])
					to_remove.append(rep2[i])

				if rep2[i][0] < rep2[i+1][1] and rep2[i][1] > rep2[i+1][1] \
				and rep2[i+1][0] < rep2[i][0] :
					nb += 1
					#print(rep2)
					#print(rep2[i], rep2[i+1])
					debut = int(rep[i+1][0])
					fin = int(rep[i][1])
					new = [debut, fin]
					to_modif.append(rep2[i], rep2[i+1])
					to_change.append(new)

				if rep2[i][1] < rep2[i+1][1] and rep2[i+1][0] > rep2[i+1][0] \
				and rep2[i][1] > rep2[i+1][0] :
					nb += 1
					#print(rep2)
					#print(rep2[i], rep2[i+1])
					debut = int(rep2[i][0])
					fin = int(rep2[i+1][1])
					new = [debut, fin]
					to_modif.append(rep2[i], rep2[i+1])
					to_change.append(new)

			'''
			if to_remove :
				print("REPT À SUPP ----> ", to_remove, len(to_remove))
			if to_change :
				print("POSITIONS À MODIFIER ----> ", to_change)
			#else :
			#	print("type de chevauchement inexistant")
			'''
			if nb != 0 :
				#print("nb de chevauchement dans la séquence :", nb)
				k += 1
				for elem in to_remove :
				#	print("REMOVING")
				#	print(elem)
					rep2.remove(elem)
				#	print(rep2, len(rep2))
				#print("\n============================\n")

			dic['pos'] = rep2

			'''
			for i in range(len(rep2)) :
				#print(rep2[i])
				if [rep2[i], rep2[i+1]] in to_modif :
					rep2.remove(rep2[i+1])
					rep2.remove(rep2[i])
			rep2.append(to_change)
			rep2 = sorted(rep2, key = itemgetter(0))
			dic['pos'] = sorted(rep2, key = itemgetter(0))
			'''
			#print(to_change)

		print('nb de seq avec chevauchement :', k)
		

	dico_all = Proteom_all(path_proteom_all)

	for idt, dic in dico.items() :
		if idt in dico_all.keys() :
			dic['seq'] += dico_all[idt]

	for idt, dic in dico.items() :
		prop = []
		lgr = []
		for rep in dic['pos'] :
			d = rep[0]
			f = rep[1]
			seq_ = dic['seq'][d:f]
			prop.append(prop_calculator(seq_))
			l = f - d
			lgr.append(l)
		dic['aa_prop'] = prop
		dic['l_rep'] = lgr
		if sum(lgr) != 0 and len(dic['seq']) != 0 :
			dic['rep_prop'] = sum(lgr)/len(dic['seq'])
		else :
			dic['rep_prop'] = 0
		del dic['seq']

	return dico


def Proteom_all(path) :
	''' Function that read the new meta proteom that we'll use in other functions 
	
	Parameters
	----------
	path : str
		Path to the New_Proteom_All.txt file

	Returns
	-------
	dico : dict
		Dictionnary of the metaproteom (key = id, value = sequence)

	'''

	dico = {}
	#with open(path+"TMHMM/files/New_Proteom_All.txt", "r") as filin :
	with open(path+'New_Proteom_All.txt', "r") as filin :
		for line in filin : 
			if line.startswith('>') :
				idt = line.strip()
				dico[idt] = ""
			else :
				dico[idt] += line.strip()

	return dico


def prop_calculator(sequence) :
	''' Calculator of the amino acid proportion in a sequence
	
	Parameters
	----------
	sequence : str
		Sequence on which we want to calculates its amino acid frequency

	Returns
	-------
	freq_dico : dict
		Dictionnary of the frequency (key = amino acid, value = its proportion
		within the sequence)

	'''


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
	

def verif() : 
	''' Verify if all the protein ids has been parsed (some errors in radar files
	may make the parsing incomplete or wrong, linked to the text editor)
	
	Parameters
	----------
	None

	Returns
	-------
	None

	'''
	
	fich = path_output+"output_Analyzeseq_ToKeep_tem_neg.tmhmm.txt"
	idt_radar = radar(path_pb)

	idt_keep = []

	with open(fich, 'r') as filin :
		for line in filin :
			idt_keep.append(line.strip())

	idt_not = []
	for idt in idt_keep :
		idt = '>'+idt
		if idt not in idt_radar :
			idt_not.append(idt)
	print("-----------", idt_not, len(idt_not))


def Data_Create(pattern_ard2, pattern_ard2_2, pattern_radar, pattern_radar_2) :
	''' Function that reads each outputs of each softwares and perform 
		the corresponding function to parse it

	Parameters
	----------
	None

	Returns
	-------
	dico_trgp2, dico_wlf, dico_ard2, dico_loca, dico_dploc, dico_radar : dict
		Dictionnaries of dictionnary (key = negative/positive set, 
		value = results of the parsing (values) for each prot id (keys)

	'''

	file_ard2 = glob.glob(path_ard2+pattern_ard2)
	file_radar = glob.glob(path_radar+pattern_radar)

	dico_ard2 = {}
	for file in file_ard2 :
		dico_ard2[basename(file)] = ard2(file, pattern_ard2_2)	


	dico_radar = {}
	for file in file_radar :
		dico_radar[basename(file)] = radar(file, pattern_radar_2)


	#print(dico_ard2.values())
	#print(dico_ard2)
	#print(dico_radar)
	return dico_ard2, dico_radar




def dataframe_maker(dico_ard2, dico_radar, pattern_prot_all) :
	''' Construction of the dataframe that contains all the results for
		each sequence (index = sequence, columns = results of the parsing
		for a software)
		The 'type' column correspond to the type of the set 
		(0 --> positive samples/ 1 --> negative samples)


	Parameters
	----------
	dico_trgp2, dico_wlf, dico_ard2, dico_loca, dico_dploc, dico_radar : dict
		Same dictionnaries returned in the Data_Create() function

	Returns
	-------
	df : Dataframe
		Dataframe of the parsing results for each sequence

	'''

	dico_all = Proteom_all(pattern_prot_all)
	print(dico_all)
	idt_all = []

	for idt in dico_all.keys() :
		idt_all.append(idt)

	idt_all = np.array(idt_all)

	col = ['type', 'ard2', 'radar']
	df = pd.DataFrame(0, index = idt_all, columns = col)



	for prote, dic in dico_ard2.items() :
		idt_l = []
		#i = 0
		for idt, val in dic.items() :
			idt_l.append(idt)
		#if prote == 'neg' :
		#	for idt in df.index :
		#		idt = str(idt)
		#		if idt in idt_l :
		#			df.loc[idt, 'type'] = 1
		#else : 
		#	for idt in df.index :
		#		df.loc[idt, 'type'] = i
		#		i += 1
				

	for prote, dic in dico_ard2.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'ard2'] = [res]	

	for prote, dic in dico_radar.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'radar'] = [res]


	print("---------df1")
	print(df)
	df.to_csv('dataframe_first.csv', sep = '\t', header = True, index = True)

	#print("------------------------DICO ard2------------------------")
	#print(dico_ard2.values())
	#print("------------------------END DICO ard2------------------------")
	return df


def Modification(dataframe, pattern_prot_all) :
	''' Function that modify the dataframe of the parsing results 
		in order to replace its results by numerics as floats

		- 'type' --> 0 = positive samples / 1 = negative sample / else = 2
		- 'trp2' --> NoTP = 0 / cTP = 1 / mTP = 2 / SP = 3 / LuTP = 4
		- 'wolfpsort' --> chlo = 1 / mito = 2 / else = 0
		- 'Localizer' --> if 'Y' = 1 / else = 0
		- deeploc --> if probability to be adressed to chlo or mito = 1 / else = 0
		- 'ard2' --> if the proteins has at least 3 linkers at a distance of 30 aa = 1 / else = 0
		- 'radar' --> if the sequence has repetitions where 29 < lrep < 46 = 1 / else = 0
		(29 < lrep < 46 can corresponds to PPR, OPR or TPR)

	
	Parameters
	----------
	dataframe : Dataframe
		Dataframe of results that we need to modificate

	Returns
	-------
	dataframe : Dataframe
		Dataframe of the modified parsing results for each sequence 

	'''

	dico_all = Proteom_all(pattern_prot_all)
	
	for col in dataframe :
		print(dataframe[col], type(dataframe[col]))	

	for elem in dataframe['type'] :
		dataframe.replace(elem, float(elem), inplace = True)

	for index, elem in enumerate(dataframe['ard2']) :
		nb = 0
		#print(elem)
		for dic in elem :
			l = []
			#print(dic)
			for aa, link in dic.items() :
				#print(link)
				nb += len(link)
				for i in range(len(link)) :
					l.append(link[i][0])
			if l :
				mini = min(l)
				maxi = max(l)
				med = statistics.median(l)
			else :
				mini = 0
				maxi = 0
				med = 0

			#dataframe['ard2'].iloc[index] = [nb, mini, med, maxi]
			#print(dataframe.index[index])
			#for idt, seq in dico_all.items() :
			#	if idt == dataframe.index[index] :
			#		s = seq

			#prop = nb/len(s)
			dataframe['ard2'].iloc[index] = [nb, mini, med, maxi]
			#dataframe['ard2'].iloc[index] = prop

	'''
	for index, elem in enumerate(dataframe['ard2']) :
		k = 0
		for dic in elem :
			l = []
			for aa, link in dic.items() :
				for lst in link :
					l.append(lst)
			for i in range(len(l)-1) :
				if l[i+1][1] - l[i][1] > 30 :
					k += 1
			if k >= 3 : 
				dataframe['ard2'].iloc[index] = float(1)
			else :
				dataframe['ard2'].iloc[index] = float(0)
	'''


	for index, elem in enumerate(dataframe['radar']) :
		nb = 0
		if type(elem) != type(['liste']) :
			dataframe['radar'].iloc[index] = [elem]
		else :
			for dic in elem :
				if len(dic['l_rep']) == 0 :
					dataframe['radar'].iloc[index] = 0
				else : 
					for item in dic['l_rep'] :
						if item >= 29 and item <= 46 :
							nb += 1
					for idt, seq in dico_all.items() :
						if idt == dataframe['radar'].index[index] :
							s = len(seq)
				prop = nb/s
				dataframe['radar'].iloc[index] = prop
				#if prop != 0 :
				#	dataframe['radar'].iloc[index] = float(1)
				#else : 
				#	dataframe['radar'].iloc[index] = float(0)


	for index, elem in enumerate(dataframe['radar']) :
		if type(elem) != type(float(1)) :
			print(elem, type(elem))
			dataframe['radar'].iloc[index] = float(elem[0])


	dataframe['ard2_mini'] = float(0)
	dataframe['ard2_med'] = float(0)
	dataframe['ard2_max'] = float(0)
	for index, elem in enumerate(dataframe['ard2']) :
		dataframe['ard2'].iloc[index] = float(elem[0])
		dataframe['ard2_mini'].iloc[index] = float(elem[1])
		dataframe['ard2_med'].iloc[index] = float(elem[2])
		dataframe['ard2_max'].iloc[index] = float(elem[3])



	for col in dataframe :
		print(dataframe[col], type(dataframe[col]))

	print("---------df2")
	print(dataframe)

	dataframe.to_csv('dataframe_interm.csv', sep = '\t', header = True, index = True)

	return dataframe

def writing(df) :

	df.to_csv('dataframe_all.csv', sep = '\t', header = True, index = True)


def splitting(df) :
	''' Here we split the dataframe into 2 dataframes : 
	- 1 contains postive samples
	- the other contains negative samples

	Parameters
	----------
	df : Dataframe
		Dataframe to split

	Returns
	-------
	df_pos, df_neg : new dataframes of positive and negative samples


	'''


	df_pos = df[df['type'] == 0]
	df_neg = df[df['type'] == 1]


	print("-----df pos")
	print(df_pos)
	print("-----df neg")
	print(df_neg)

	df_pos.to_csv('dataframe_pos.csv', sep = '\t', header = True, index = True)
	df_neg.to_csv('dataframe_neg.csv', sep = '\t', header = True, index = True)

	return df_pos, df_neg




def Tsne_by_tools(dataframe) :
	''' Perform a tsne on the modified dataframe 

	Parameters
	----------
	dataframe : Dataframe
		Dataframe to perform t-sne on

	Returns
	-------
	None

	Plot
	----
	A t-sne of the dataframe for each columns


	'''
	arr_list = []
	for data in dataframe :
		x = dataframe[data]
		data = np.array(x)
		arr_list.append(data)
	print(arr_list, len(arr_list))


	tsne = []
	i = 0
	for array in arr_list :
		print(array, type(array), array.shape, len(array))
		for i in range(len(array)) :
			if type(array[i]) == type([1]) :
				print(array, array[i])
				array[i] = np.asarray(array[i])

		array = array.reshape(-1, 1)
		tsne.append(tsne_data(array))
		#i += 1



	label = list(dataframe.columns)
	print(label, type(label))
	
	print("--------------TSNE PERFORMING--------------")
	for data in tsne :
		sns.scatterplot(x = 'x', y = 'y', data = data)
	plt.legend(label, prop = {'size' : 5.7})
	if dataframe is df_pos :
		plt.title("Tsne of positive set", fontsize = 15)
	if dataframe is df_neg :
		plt.title("Tsne of negative set", fontsize = 15)
	plt.show()



def tsne_data(to_data) :
	''' Computes the values to perform a tsne

	Parameters
	----------
	to_data : array
		Data to compute a tsne on

	Returns
	-------
	data : array
		Results of the computed values of a tsne

	'''

	tsne = TSNE(n_components = 2, random_state = 0, perplexity = 50) #### OR default perplex = 30

	X_2d = tsne.fit_transform(to_data)

	data = pd.DataFrame()
	data['x'] = X_2d[:,0]
	data['y'] = X_2d[:,1]

	return data


def Tsne_all(pos, neg) :

	print("---------------TSNE POS/NEG---------------")
	arr_list = []
	tsne = []
	list_df = [pos, neg]

	for df in list_df :
		if df is pos :
			print("df pos")
			print(df)
		elif df is neg :
			print("df neg")
			print(df)

		#for index, elem in enumerate(df['ard2']) :
		#	df['ard2'].iloc[index] = elem[to_plot] 

		df = np.array(df)
		print(df, len(df), type(df))

		#for arr in df :
		#	print(arr, len(arr), type(arr))


		tsne.append(tsne_data(df))
		label = list(pos.columns)
		print(label, type(label))

	
		print("--------------TSNE PERFORMING--------------")
		for data in tsne :
			sns.scatterplot(x = 'x', y = 'y', data = data)
	if list_df[0] is pos :
		plt.legend(['pos', 'neg'], prop = {'size' : 5.7})
	elif list_df[0] is neg :
		plt.legend(['neg', 'pos'], prop = {'size' : 5.7})
	plt.title("Tsne on positive and negative sets", fontsize = 15)
	plt.show()
	





def Prop_Test(df1, df2, fold, col, to_plot) :
	''' Calculates and computes a proportion test between two dataframes

	Parameters
	----------
	df1, df2 : Dataframes
		Dataframes to compute the test on

	fold : int
		Fold of pvalue (if computed pvalue > fold --> test is significant)

	col : str (' ')
		Column name of the results to be performed

	'''

	#if col == 'ard2' :
	#	for index, elem in enumerate(df1[col]) :
	#			df1['ard2'].iloc[index] = elem[to_plot]
	#	for index, elem in enumerate(df2[col]) :
	#			df2['ard2'].iloc[index] = elem[to_plot]

	m_df1 = df1[col].mean()
	m_df2 = df2[col].mean()
	total = m_df1 + m_df2
	print("Moyenne prop linker des df :", m_df1, m_df2)
	print("H0 : même proportion\n", "H1 : proportions différentes")

	# Conditions of application
	assert len(df1)*m_df1 > 10, "Condition non validée"
	assert len(df2)*m_df2 > 10, "Condition non validée"
	assert len(df1)*(1-m_df1) > 10, "Condition non validée"
	assert len(df2)*(1-m_df2) > 10, "Condition non validée"

	print("Conditions de validité passées")

	# Computation of the residual standard error
	variance = (total*(1-total))
	#variance = abs(total*(1-total))

	std_error = np.sqrt(variance*(1/len(df1) + 1/len(df2)))
	print("Erreur standard = ", std_error)

	# Calculation of the statistic test
	best = m_df1 - m_df2
	print("Estimation :", best)
	hyp_estimation = 0
	test_stat = (best - hyp_estimation)/std_error
	print("Résultat du test statistique : ", test_stat)
	print("Difference : ", test_stat - best)

	# Computation of pvalue
	pvalue = 2*dist.norm.cdf(-np.abs(test_stat))
	print("Pvalue = ", pvalue)

	if pvalue <= fold : 
		print("Pvalue < 0.05 --> pas de différence significative")
	else : 
		print("Pvalue > 0.05 --> existe une différence significative")


def Mean_test(df1, df2, col) :
	
	print("Mean Test (Student t-test)")
	print("----------------", col, "----------------")

	#Control of the normaliy of the samples --> pvalue of shapiro test must be > fold
	# (null hypothesis of a normal distribution)
	x1, pval1 = shapiro(df1[col])
	x2, pval2 = shapiro(df2[col])


	#Control of the variance equality
	# Fisher-Snedecor F-test
	f(df1[col], df2[col])

	#Anova unidirectional
	stats.f_oneway(df1[col], df2[col])

	#We compare the samples in order to see if there is a significant difference
	y = stats.ttest_ind(df1[col], df2[col])
	print("Computed pvalue : ", y[1])
	pvalue = y[1]

	if pvalue < 0.05 :
		if pvalue < 0.001 : 
			print("pvalue significant --> there is a difference, ***")
		elif pvalue < 0.01 : 
			print("pvalue significant --> there is a difference, **")
		else :
			print("pvalue significant --> there is a difference, *")
	else : 
		print("There is no significant difference")



def Sep_long_proteom(path, pattern1, pattern2, fold) :
	
	proteom = glob.glob(path+pattern1)
	print(len(proteom))

	os.chdir(path+pattern2)

	for p in proteom :
		dico = {}
		print("-", basename(p))
		with open(p, 'r') as filin :
			for line in filin :
				if line.startswith('>') :
					idt = line.strip()
					dico[idt] = ""
				else :
					dico[idt] += line.strip()
		print(len(dico.keys()))
		if len(dico.keys()) > fold :
			print("File too big ----->", basename(p), len(dico.keys()))
			new_dic = {}
			i = 0
			k = 0
			id_list = []
			for idt, seq in dico.items() :
				id_list.append(idt)
				if len(seq) < 120 :
					print("!!!!LEN(SEQ) < 120!!!!")
					print(idt, len(seq))
			for idt in dico.keys() :
				if i <= fold :
					new_dic[idt] = dico[idt]
					i += 1
					if idt == id_list[-1] :
						k += 1
						with open('sep_proteom_'+basename(p)+'_'+str(k)+'.txt', 'w') as filout :
							for idt, seq in new_dic.items() :
								filout.write(idt+"\n"+seq+"\n")
						print(len(new_dic.keys()))
				elif i > fold :
					i = 0
					k += 1
					new_dic[idt] = dico[idt]
					with open('sep_proteom_'+basename(p)+'_'+str(k)+'.txt', 'w') as filout :
						for idt, seq in new_dic.items() :
							filout.write(idt+"\n"+seq+"\n")
					print(len(new_dic.keys()))
					new_dic = {}
					#i += 1


def concat_ard2(path, pattern1, pattern2, pattern3) :
	
	res = glob.glob(path+pattern1)
	res.sort()
	print(res, len(res), type(res))

	os.chdir(path+pattern3)
	new_file = []


	for f in res :
		with open(f, 'r') as filin :
			for line in filin :
				new_file.append(line)
	print(len(new_file), type(new_file))
	#print(new_file)

	a = 0
	for i in range(len(new_file)) :
		if new_file[i].startswith('>') :
			a += 1
			new_file[i] = '>'+str(a)+'\t R\n'
	#print(new_file)

	with open('STDOUT_'+pattern2+'.txt', 'w') as filout :
		for line in new_file :
			filout.write(line)



def rows_acc(path, file) :

	files = glob.glob(path+file)
	files.sort()
	print(files)

	idt = []

	for f in files :
		with open(f, 'r') as filin :
			for line in filin :
				#pass
				line = line.strip()
				line = line.split('"')[1]
				line = '>'+line
				if line != '>1' :
					idt.append(line)
	
	return idt


def add_df(idt, pattern1, path, file) :
	
	print("---------df add---------")
	#df['acc'] = 0

	df = pd.read_csv(path+file, sep = '\t')

	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']


	list_idt = list(df.index)

	print(df)

	to_del = []
	for ident in idt :
		if ident not in list_idt :
			to_del.append(ident)
	#print(to_del, len(to_del))

	with open(path_output+'to_del_acc.txt', 'w') as filout :
		for elem in to_del :
			filout.write(elem+'\n')

	to_del_ind = []
	z = 0
	for ident in idt :
		if ident in to_del :
			to_del_ind.append(int(z))
		z += 1

	print(df)
	print(to_del_ind, len(to_del_ind))
	
	#print(df[df['type'] == 0])
	#print(list_idt, len(list_idt))


	print("---------------ACC---------------")

	for i in range(36) :
		df['acc'+str(i)] = 0


	
	fich = glob.glob(path_output+pattern1)
	print(fich)
	fich.sort()
	print(fich)

	k = 0
	element = []
	for f in fich :
		print(basename(f))
		with open(f, 'r') as filin :
			for line in filin :
				if k not in to_del_ind :
					element.append(line)
				k += 1
	print(len(element))

	#print(element)
	with open(path_output+'ACC/acc_file_all.txt', 'w') as filout :
		for elem in element :
			filout.write(elem)

	k = 0
	with open(path_output+'ACC/acc_file_all.txt', 'r') as filin :
		for line in filin :
			line = line.split('\t')
			line[-1] = line[-1].strip()
			#print(line, len(line))
			for i in range(len(line)) :
				df['acc'+str(i)].iloc[k] = float(line[i])
			k += 1

	#for col in df :
	#	print(df[col])
	#print(df)

	print("---------------FREQ---------------")

	for aa in list_of_aa :
		df[aa] = 0

	dico_all = Proteom_all(path_output)

	#print(dico_all)
	#print(len(dico_all.keys()), len(df))


	dico = {}
	for idt, seq in dico_all.items() :
		dico[idt] = ""
		dico[idt] = freq_aa(seq)
			
	#print(dico)
	#print(len(dico.keys()), len(df))


	for idt, dic in dico.items() :
		for index in list_idt :
			if idt == index :
				for col in df :
					for aa, prop in dic.items() :
						if col == aa :
							df.loc[index, aa] = prop
	

	for col in df :
		print(df[col])

	print(df)
	df.to_csv('dataframe_all.csv', sep = '\t', header = True, index = True)


	return df


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



if __name__ == '__main__' :

	path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/proteomes_diatom/outputs/TMHMM/Pour_Celine/"
	path_output = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/proteomes_diatom/outputs/"
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']
	
	os.chdir(path_output)


	# Phaedodactylum

	# ARD2
	path_ard2 = path_output+"ARD2/Phaedodactylum/"

	# RADAR
	path_radar = path_output+"RADAR/Phaedodactylum/"

	results_ard2, results_radar = Data_Create("STDOUT_Phaedodactylum", 'TMHMM/Pour_Celine/*.txt', 'Phaedodactylum_RADAR.txt', 'TMHMM/files/')
	final_results = dataframe_maker(results_ard2, results_radar, path_output+'TMHMM/files/')
	df = Modification(final_results, path_output+'TMHMM/files/')
	#idt = rows_acc(path_output, 'ACC/rownames_*')
	#df_f = add_df(idt, 'ACC/Acc_output_*', path_output, 'dataframe_interm.csv')
	#writing(df_f)

	'''
	# All Diatoms 

	# ARD2
	path_ard2 = path_output+"ARD2/*/"
	#path_tmhmm = path_output+"TMHMM/"
	#Long_prot_sep = Sep_long_proteom(path_output, 'TMHMM/New_prot/*.txt', 'TMHMM/New_prot/Separated/', int(25000))
	#concat(path_output, 'ARD2/Arabi/old_1_2/*/_STDOUT_*', 'Arabi', 'ARD2/Arabi/concat/')
	
	# WOLFPSORT
	#path_wpsort = path_output+"WOLFPSORT/*/"

	# TARGETP2
	#path_trgp2 = path_output+"TARGETP2/*/"

	# DEEPLOC
	#path_dploc = path_output+"DEEPLOC/*/"

	# LOCALIZER
	#path_loca = path_output+"LOCALIZER/*/"

	# RADAR
	path_radar = path_output+"RADAR/*/"

	results_ard2, results_radar = Data_Create("STDOUT_"+"*", 'TMHMM/Pour_Celine/*.txt', '*'+'RADAR', 'TMHMM/files/')
	final_results = dataframe_maker(results_ard2, results_radar, path_output+'TMHMM/files/')
	df = Modification(final_results, path_output+'TMHMM/files/')
	#idt = rows_acc(path_output, 'ACC/rownames_*')
	#df_f = add_df(idt, 'ACC/Acc_output_*', path_output, 'dataframe_interm.csv')
	#writing(df_f)
	#tsne = Tsne(df_f)
	'''














