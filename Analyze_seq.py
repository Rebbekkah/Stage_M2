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
Stage M2 - supervisor : Ingrid Lafontaine
------------------------------------------------------------------

"""

# Modules
import pandas as pd
import numpy as np
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

	print("Total prot : ", idt_all, len(idt_all))
	print("PROT TO KEEP : ", id_to_keep, len(id_to_keep))

	id_to_delete = []	
	for idt in idt_all :
		if idt in dico_h.keys() :
			id_to_delete.append(idt)
	print("TO DELETE : ", id_to_delete)
	
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
			print(idt)
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

	Returns
	-------
	None

	Writes
	------
	New_Proteom.txt : txt file
		File of the new proteoms (negative and positive sets)

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
		for idt in idt_ok[i] :
			idt = ">"+idt
			if idt in dico_all.keys() :
				dico2[idt] = ""
				dico2[idt] += dico_all[idt]
		print("LEN DICO2", len(dico2))
		with open("New_Proteom_"+basename(p)+".txt", "w") as filout :
			for idt, seq in dico2.items() :
				if len(seq) > 120 :
					filout.write(idt+"\n"+seq+"\n")
		i += 1


def ard2(file) :
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


	proteoms = glob.glob(path_tmhmm+'*.fasta_line.txt')

	dico_else = {}
	dico_pos = {}
	dico_neg = {}
	for p in proteoms :
		with open(p, "r") as filin :
			
			if basename(p) == 'New_Proteom_1196_tem_neg.fasta_line.txt' :
				for line in filin :
					if line.startswith('>') :
						idt = line.strip()
						dico_neg[idt] = ""
					else :
						dico_neg[idt] += line.strip()
			elif basename(p) == 'New_Proteom_1081_tem_pos.fasta_line.txt' :
				for line in filin :
					if line.startswith('>') :
						idt = line.strip()
						dico_pos[idt] = ""
					else :
						dico_pos[idt] += line.strip()
			'''
			else : 
			
			for line in filin :
				if line.startswith('>') :
					idt = line.strip()
					dico_else[idt] = ""
				else :
					dico_else[idt] += line.strip()
			'''
	dico_f = {}
	link = []
	
	if basename(file) == 'STDOUT_neg' :
		for idt in dico_neg.keys() :
			dico_f[idt] = {}

		for linker in dico.values() :
			link.append(linker)

		for index, key in enumerate(dico_f) :
			dico_f[key] = link[index]

		
	elif basename(file) == 'STDOUT_pos' :
		for idt in dico_pos.keys() :
			dico_f[idt] = {}

		for linker in dico.values() :
			link.append(linker)

		for index, key in enumerate(dico_f) :
			dico_f[key] = link[index]
	
	'''
	else : 
	
	if basename(file) == 'STDOUT' :
		for idt in dico_else.keys() :
			dico_f[idt] = {}

		for linker in dico.values() :
			link.append(linker)
	
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


def radar(file) :
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
		

	dico_all = Proteom_all(path_output)

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
	with open(path+"tmhmm_filtred/New_Proteom_All.txt", "r") as filin :
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


def Data_Create() :
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

	file_ard2 = glob.glob(path_ard2+"STDOUT_"+"*")
	file_wlf = glob.glob(path_wpsort+"*.wolfpsort")
	file_trgp2 = glob.glob(path_trgp2+"short_output_"+"*")
	file_dploc = glob.glob(path_dploc+"*"+"deeploc"+"*")
	file_loca = glob.glob(path_loca+'*'+'localizer')
	file_radar = glob.glob(path_radar+'*'+'radar')

	dico_trgp2 = {}
	for file in file_trgp2 :
		if basename(file) == 'short_output_neg' :
			dico_trgp2['neg'] = {}
			dico_trgp2['neg'] = targetp2(file)
		elif basename(file) == 'short_output_pos' :
			dico_trgp2['pos'] = {}
			dico_trgp2['pos'] = targetp2(file)
		else : 
			dico_trp2[basename(file)] = targetp2(file)


	dico_wlf = {}
	for file in file_wlf :
		if basename(file) == 'output_neg.wolfpsort' :
			dico_wlf['neg'] = {}
			dico_wlf['neg'] = wolfpsort(file)
		elif basename(file) == 'output_pos.wolfpsort' :
			dico_wlf['pos'] = {}
			dico_wlf['pos'] = wolfpsort(file)
		else : 
			dico_wlf[basename(file)] = wolfpsort(file)


	dico_ard2 = {}
	for file in file_ard2 :
		if basename(file) == 'STDOUT_neg' :
			dico_ard2['neg'] = {}
			dico_ard2['neg'] = ard2(file)
		elif basename(file) == 'STDOUT_pos' :
			dico_ard2['pos'] = {}
			dico_ard2['pos'] = ard2(file)
		else : 
			dico_ard2[basename(file)] = ard2(file)	


	dico_loca = {}
	for file in file_loca : 
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_localizer' :
			dico_loca['neg'] = {}
			dico_loca['neg'] = localizer(file)
		elif basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_localizer' :
			dico_loca['pos'] = {}
			dico_loca['pos'] = localizer(file)
		else : 
			dico_loca[basename(file)] = localizer(file)
	

	dico_dploc = {}
	for file in file_dploc :
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_deeploc.txt' :
			dico_dploc['neg'] = {}
			dico_dploc['neg'] = deeploc(file)
		elif basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_deeploc.txt' :
			dico_dploc['pos'] = {}
			dico_dploc['pos'] = deeploc(file)	
		else : 
			dico_dploc[basename(file)] = deeploc(file)


	dico_radar = {}
	for file in file_radar :
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_radar' :
			dico_radar['neg'] = {}
			dico_radar['neg'] = radar(file)
		elif basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_radar' :
			dico_radar['pos'] = {}
			dico_radar['pos'] = radar(file)
		else : 
			dico_radar[basename(file)] = radar(file)


	return dico_trgp2, dico_wlf, dico_ard2, dico_loca, dico_dploc, dico_radar




def dataframe_maker(dico_trgp2, dico_wlf, dico_ard2, dico_loca, dico_dploc, \
	dico_radar) :
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

	dico_all = Proteom_all(path_output)
	idt_all = []

	for idt in dico_all.keys() :
		idt_all.append(idt)

	idt_all = np.array(idt_all)

	col = ['type', 'trp2', 'wolfpsort', 'ard2', 'localizer', 'deeploc', \
	'radar']
	df = pd.DataFrame(0, index = idt_all, columns = col)



	for prote, dic in dico_ard2.items() :
		idt_l = []
		#i = 0
		for idt, val in dic.items() :
			idt_l.append(idt)
		if prote == 'neg' :
			for idt in df.index :
				idt = str(idt)
				if idt in idt_l :
					df.loc[idt, 'type'] = 1
		#else : 
		#	for idt in df.index :
		#		df.loc[idt, 'type'] = i
		#		i += 1
				

	for prote, dic in dico_trgp2.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'trp2'] = res

	for prote, dic in dico_wlf.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'wolfpsort'] = res

	for prote, dic in dico_ard2.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'ard2'] = [res]	

	for prote, dic in dico_loca.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'localizer'] = res

	for prote, dic in dico_dploc.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'deeploc'] = [res]

	for prote, dic in dico_radar.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'radar'] = [res]


	print("---------df1")
	print(df)
	#print("------------------------DICO WOLFPSORT------------------------")
	#print(dico_wlf)
	#print("------------------------END DICO WOLFPSORT------------------------")
	return df


def Modification(dataframe) :
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

	dico_all = Proteom_all(path_output)
	
	for col in dataframe :
		print(dataframe[col], type(dataframe[col]))	

	for elem in dataframe['type'] :
		dataframe.replace(elem, float(elem), inplace = True)

	possible = []
	for res in dataframe['trp2'] :
		if res not in possible :
			possible.append(res)


	for index, poss in enumerate(possible) :
		for elem in dataframe['trp2'] :
			if elem == poss :
				dataframe.replace(elem, float(index), inplace = True)
	'''
	for elem in dataframe['wolfpsort'] :
		if elem == 'chlo' or elem == 'mito' :
			dataframe.replace(elem, float(1), inplace = True)
		else : 
			dataframe.replace(elem, float(0), inplace = True)
	'''

	for elem in dataframe['localizer'] :
		if elem == 'Y' :
			dataframe.replace(elem, float(1), inplace = True)
		else :
			dataframe.replace(elem, float(0), inplace = True)
	

	for index, elem in enumerate(dataframe['deeploc']) :
		probability = []
		for key, prob in elem[0].items() :
			probability.append(float(prob))
		p = max(probability)
		dataframe['deeploc'].iloc[index] = p
		#if p == 0 :
		#	dataframe['deeploc'].iloc[index] = p
		#else : 
		#	dataframe['deeploc'].iloc[index] = p


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

	dataframe.to_csv('dataframe_all.csv', sep = '\t', header = True, index = True)

	return dataframe



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


def Tsne_all(pos, neg, to_plot) :

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

		for index, elem in enumerate(df['ard2']) :
			df['ard2'].iloc[index] = elem[to_plot] 

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

	if col == 'ard2' :
		for index, elem in enumerate(df1[col]) :
				df1['ard2'].iloc[index] = elem[to_plot]
		for index, elem in enumerate(df2[col]) :
				df2['ard2'].iloc[index] = elem[to_plot]

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


def Sep_long_proteom(pattern1, pattern2, fold) :
	
	proteom = glob.glob(path_output+pattern1)
	print(len(proteom))

	os.chdir(path_output+pattern2)

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


def add_df(df) :
	
	print("---------df add---------")
	#df['acc'] = 0
	list_idt = list(df.index)

	print(df)
	#print(list_idt, len(list_idt))

	for i in range(36) :
		k = i
		df['acc'+str(k)] = 0


	k = 0
	fich = glob.glob(path_output+'acc/v2/Acc_output_*')

	print(fich)
	if fich[0] == path_output+'acc/v2/Acc_output_New_Proteom_1196_tem_neg.fasta_line.txt.txt' :
		fich = fich[::-1]
	print(fich)

	for f in fich :
		#k = 0
		print(basename(f))
		with open(f, 'r') as filin :
			#df_acc = pd.read_csv(f, sep = '\t')
			for line in filin :
				
				#lt = []
				#df['acc'+str(k)] = 0
				line = line.split('\t')
				#print(line)
				line[-1] = line[-1].strip()
				#print(line, len(line))
				#lt.append(line)
				#print(lt, len(lt))
				#k+=1

				for i in range(len(line)) :
					df['acc'+str(i)].iloc[k] = line[i]


				'''
				for i in range(len(line)) :
					if basename(f) == 'Acc_output_New_Proteom_1196_tem_neg.fasta_line.txt.txt' :
						if df['type'].iloc[k] == 1 :
							df['acc'+str(i)].iloc[k] = line[i]
					elif basename(f) == 'Acc_output_New_Proteom_1081_tem_pos.fasta_line.txt.txt' :
						if df['type'].iloc[k] == 0 :
							df['acc'+str(i)].iloc[k] = line[i]
				'''
				'''
				if basename(f) == 'Acc_output_New_Proteom_1196_tem_neg.fasta_line.txt.txt' :
					for i in range(len(line)) :
						if df['type'].iloc[k] == 1 :
							#print("ok1")
							df['acc'+str(i)].iloc[k] = line[i]
							#print("----------------------")
					#print(df.iloc[k])
				elif basename(f) == 'Acc_output_New_Proteom_1081_tem_pos.fasta_line.txt.txt' :
					#print("oui1")
					#print(df['type'])
					for i in range(len(line)) :
						#print(df['type'].iloc[k] == 0)
						if df['type'].iloc[k] == 0 :
							#print(df.iloc[k])
							print("oui2")
							#print("ok2")
							df['acc'+str(i)].iloc[k] = line[i]
					#print(df.iloc[k])
				'''
				k += 1
				

			#print(df_acc)
		#if basename(f) == 'Acc_output_New_Proteom_1196_tem_neg.fasta_line.txt.txt' :
		#	pass

				#for idt in list_idt :
				#	for ind in range(len(list_idt)) :


	#for col in df :
	#	print(df[col])

	print(df)


def Plotting_by_col(df, to_plot) :

	#for index, elem in enumerate(df['ard2']) :
	#	df['ard2'].iloc[index] = elem[to_plot]

	#print(df['ard2'])

	sns.set(style = "darkgrid")

	#for elem in df['ard2'] :
		#print(elem[to_plot])
	#	sns.histplot(data = df, x = elem[to_plot], y = df.index)
	#plt.show()


	for elem in df['radar'] :
		if type(elem) != type(float(1)) :
			print(elem, type(elem))


	for column in df.columns :
		#values = [0:1]
		plt.figure()
		print("--------------", column, "--------------")
		if column != 'ard2' :
			#sns.histplot(df[column])
			sns.countplot(df[column])
			#if column == 'trp2' : 
				#plt.legend(['NoTP', 'cTP', 'mTP', 'LuTP'])
				#plt.legend()
			if column == 'wolfpsort' :
				#values = list(range(2))
				#plt.xticks(ticks = values)
				#plt.xlim([0, 1])
				#plt.xlim(np.arange(0, 1, step = 0.001))
				plt.xticks(np.arange(0, 1, step = 0.002))
				plt.ylim(0, 100)
				#plt.xticks(rotation = 'vertical')
			if column == 'radar' :
				#plt.xlim([0, 1])
				#plt.bar(x = column, col = 'darkred', height = )
				#sns.countplot(df[column], color = 'darkred')
				plt.xticks(np.arange(0, 1, step = 0.002))
				#plt.xticks(rotation = 'vertical')
				#plt.yticks(np.arange(0, 1))
				#plt.xlim(0, 1, step = 0.002)
				plt.ylim(0, 10)

		elif column == 'ard2' : 
			for index, elem in enumerate(df['ard2']) :
				df['ard2'].iloc[index] = elem[to_plot]
			#sns.histplot(df[column])
			sns.countplot(df[column])
			#plt.xlim([0, 1])
			plt.xticks(np.arange(0, 1, step = 0.001))
			#plt.xlim(0, 2)
			#plt.xticks(rotation = 'vertical')

		plt.legend(column, prop = {'size' : 5.7})

		if df is df_pos :
			plt.title("Histogram of the distribution of postitive samples", fontsize = 15)
		elif df is df_neg : 
			plt.title("Histogram of the distribution of negative samples", fontsize = 15)
		else :
			plt.title("Histogram of the distribution of samples", fontsize = 15)

		plt.show()



	#sns.set(style = "darkgrid")

	#sns.histplot(data = df, x = col, kde = True)
	#plt.show()

	'''
	plt.hist(df['ard2'], hist = True, kde = False, \
		hist_kws = {'edgecolor' : 'black'})
	plt.show()
	'''



def Plotting_pos_neg(df, df_pos, def_neg, col, to_plot) :

	#plt.figure()
	#plt.subplots(2, 2)

	if col == 'ard2' :
		nb_pos = []
		med_pos = []
		for index, elem in enumerate(df_pos['ard2']) :
				df_pos['ard2'].iloc[index] = elem[to_plot]
				nb_pos.append(elem[0])
				med_pos.append(elem[2])
		sns.scatterplot(data = df_pos, x = range(len(df_pos)), y = df_pos[col], color = 'darkred', label = 'pos')
		
		nb_neg = []
		med_neg = []

		for index, elem in enumerate(df_neg['ard2']) :
				df_neg['ard2'].iloc[index] = elem[to_plot]
				nb_neg.append(elem[0])
				med_neg.append(elem[2])
		sns.scatterplot(data = df_neg, x = range(len(df_neg)), y = df_neg[col], label = 'neg')
		if to_plot == 0 :
			plt.suptitle('with number of linker')
		elif to_plot == 2 :
			plt.suptitle('with median of linker')
	else :
		sns.scatterplot(data = df_pos, x = range(len(df_pos)), y = df_pos[col], color = 'darkred', label = 'pos')
		sns.scatterplot(data = df_neg, x = range(len(df_neg)), y = df_neg[col], label = 'neg')
		plt.legend()
	plt.title("Scatterplot of the distribution between positive and negative samples")
	plt.xticks([])

	plt.show()

	for index, elem in enumerate(df['ard2']) :
		df['ard2'].iloc[index] = elem[to_plot]
	sns.boxplot(x = df['type'], y = df[col], hue = df['type'])
	plt.legend()
	plt.title('Boxplot on all sequences')

	if col == 'ard2' :
		if to_plot == 0 :
			plt.suptitle('with number of linker')
		elif to_plot == 2 :
			plt.suptitle('with median of linker')
	plt.show()

	if nb_pos or nb_neg :
		sns.scatterplot(x = med_pos, y = nb_pos, color = 'darkred', label = 'pos')
		sns.scatterplot(x = med_neg, y = nb_neg, label = 'neg')
		plt.xlabel("median")
		plt.ylabel("number")
		plt.title("Scatterplot for linkers")
		plt.legend()
		plt.show()

if __name__ == '__main__' :

	path_output = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/"
	to_script = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script"
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']

	path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/"



	os.chdir(path_output)

	# TMHMM
	#proteins = listing(path_proteom, '*.tmhmm')
	#new_proteom = proteome_maker(proteins, path_proteom)
	#separateur = sep(new_proteom, path_proteom)
	#separateur = sep(path_proteom)

	# ARD2
	path_ard2 = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/ard2_outputs/"
	path_tmhmm = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/tmhmm_filtred/"
	
	# WOLFPSORT
	path_wpsort = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/wolfpsort_output/"
	path_wpsort_essai = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/wolfpsort_output/output_pos.wolfpsort"
	#wolfpsort(path_wpsort_essai)

	# TARGETP2
	path_trgp2 = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/targetp2_outputs/"

	# DEEPLOC
	path_dploc = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/outputs_deeploc/"

	# LOCALIZER
	path_loca = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/output_localizer/"

	# RADAR
	path_radar = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/output_radar/"
	path_pb = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/output_radar/idt_neg.txt"

	results_trgp2, results_wlf, results_ard2, results_loca, results_dploc, results_radar = Data_Create()
	final_results = dataframe_maker(results_trgp2, results_wlf, results_ard2, results_loca, results_dploc, results_radar)
	df = Modification(final_results)
	df_f = add_df(df)
	#df_pos, df_neg = splitting(df_f)
	#Plotting_by_col(df_pos, 2)
	#Plotting_pos_neg(df_f, df_pos, df_neg, 'ard2', 2)
	#test_of_proportion = Prop_Test(df_pos, df_neg, 0.05, 'radar', 2)
	#tsne = Tsne(df)
	#Tsne_all(df_pos, df_neg, 2)

	'''
	path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/proteomes_diatom/"
	path_output = path_proteom+"outputs/"
	print(path_output)
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']


	os.chdir(path_output+'TMHMM/files/')

	# TMHMM
	#proteins = listing(path_output, 'TMHMM/old/*.tmhmm')
	#new_proteom = proteome_maker(proteins, path_proteom, '*/*.f'+'*a')
	#separateur = sep(path_proteom, '*/*.f'+'*a', 'outputs/TMHMM/old/*.tmhmm', 'TMHMM/New_proteom_all/')
	#Long_prot_sep = Sep_long_proteom('TMHMM/New_proteom_all/*.txt', 'TMHMM/sep_prot/', int(32000))


	# ARD2
	path_ard2 = path_output+"ARD2/"
	path_tmhmm = path_output+"TMHMM/files/"
	
	# WOLFPSORT
	path_wpsort = path_output+"WOLFPSORT/"

	# TARGETP2
	path_trgp2 = path_output+"TRP2/"

	# DEEPLOC
	path_dploc = path_output+"DEEPLOC/"

	# LOCALIZER
	path_loca = path_output+"LOCALIZER/"

	# RADAR
	path_radar = path_output+"RADAR/"

	results_trgp2, results_wlf, results_ard2, results_loca, results_dploc, results_radar = Data_Create()
	final_results = dataframe_maker(results_trgp2, results_wlf, results_ard2, results_loca, results_dploc, results_radar)
	#df_f = Modification(final_results)
	#tsne = Tsne(df_f)
	'''

	#df_pos, df_neg = splitting(df_f)
	#tsne = Tsne(df_pos)
	#tsne = Tsne(df_neg)
	#test_of_proportion = Prop_Test(df_pos, df_neg, 0.05, 'ard2')

