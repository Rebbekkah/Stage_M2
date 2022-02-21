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
from os.path import basename
import glob
from operator import itemgetter


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
	
	os.chdir(path_output)
	with open("output_Analyzeseq_ToKeep_"+basename(fichier)+".txt", "w") as filout :
		for idt in id_to_keep :
			filout.write(idt+"\n")

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



def proteome_maker(ToKeep, path_proteom) :
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
	proteom = glob.glob(path_proteom+"*.fasta_line")

	os.chdir(path_output)
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
			filout.write(idt+"\n"+seq+"\n")
	
	return dico2


def sep(path_proteom) :
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
	proteom = glob.glob(path_proteom+"*.fasta_line")
	proteom = proteom[::-1]
	file = glob.glob(path_proteom+"*.tmhmm")
	print(proteom, "\n", file)

	idt_ok = []
	for f in file :
		idt_ok.append(reading(f))


	os.chdir(path_output)
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
				#print("-------------")
				#print(window, len(window))
				l = []

				for elem in window :
					if elem[1] > 0.1 :
						l.append(elem[1])

				if not not l :
					p = max(l)
					#print(window, p)

					for elem in window :
						if p == elem[1] :
							if elem[0] in dico_linker.keys() :
								dico_linker[elem[0]].extend([elem[1:]])
							else :
								dico_linker[elem[0]] = [elem[1:]]


				window = []
			dico[idt] = dico_linker

	#print(dico)
	#print("-----------------------------------------")


	proteoms = glob.glob(path_tmhmm+'*.fasta_line.txt')

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

	
	dico_f = {}
	if basename(file) == 'STDOUT_neg' :
		for idt in dico_neg.keys() :
			#print(idt)
			dico_f[idt] = {}
			print(dico_f)
			for key, linker in dico.items() :
				#print("--------------")
				#print(basename(file))
				#print(linker)
				#dico_f[idt] = {}
				#print(dico_f)
				dico_f[idt] = dico[key]
	elif basename(file) == 'STDOUT_pos' :
		print("------------------------------")
		for idt in dico_pos.keys() :
			for key, linker in dico.items() :
				#dico_f[idt] = dico_pos[idt]
				#dico_f[idt] = linker
				dico_f[idt] = dico[key]

	print(dico_f)
	print(basename(file))

	'''
	for i in range(len(window)-1) :
		p = window[i][1]
		#print(p)
		if p > 0.1 :
			k = window[i]
			if window[i+1][1] > p :
				k = window[i+1]
			#print(k[1:])
			to_keep.append(k)
#print(to_keep, len(to_keep))
for i in range(len(to_keep)) :
	if to_keep[i][0] in dico_linker.keys() :
		dico_linker[to_keep[i][0]].extend([k[1:]])
	else :
		dico_linker[to_keep[i][0]] = [k[1:]]
dico[idt] = dico_linker
	'''
				#l = []
				#for l in window :
				#for i in range(len(window)-1) :
				#	l.append(window[i][1])
				#print(l)
				#p = max(l)
				#print(p)
				#for i in range(len(window)) :
				#	if window[i][1] == p and p > 0.1 :
				#		to_keep.append(window[i])
				#	else :
				#		continue
			#print(to_keep, len(to_keep))
			#dico[idt] = to_keep

	'''
	#if window[i][1] > 0.1 :
	#	p = window[i][1] 
	#if window[i][1] > 0.1 or window[i+1][1] > window[i][1] : 
	if window[i][1] > 0.1 :
		#print(window[i], window[i+1])
		#print(window)
		#to_keep.append(window[i])
		to_keep = window[i]
		#keep.append(to_keep)
	else :
		to_keep = 0
	if window[i+1][1] > window[i][1] :
		to_keep = window[i+1]
		#to_keep.append(window[i+1])
	if to_keep != 0 :
		keep.append(to_keep)
	'''

	'''
	if l[1] > 0.1 :
		p = l[1]
		if l[1] >= p :
			keep.append(l)
	'''

	#print(keep, len(keep))

	'''
	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if line.startswith('>') :
				pos = 0
				idt = line.split()[0]
				dico[idt] = ""
				dico_linker = {} 
			else :
				pos += 1
				elem = line.split("\t")[0]
				elem = elem.split()[0]
				aa = elem[0]
				proba = elem[1:5]
				if float(proba) > 0.10 :
					if aa in dico_linker.keys() :
						dico_linker[aa].extend([proba, pos])
					else :
						dico_linker[aa] = [proba, pos]
			dico[idt] = dico_linker


	proteoms = glob.glob(path_tmhmm+'*.fasta_line.txt')

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


	dico_f = {}
	if basename(file) == 'STDOUT_neg' :
		for idt in dico_neg.keys() :
			for linker in dico.values() :
				dico_f[idt] = dico_neg[idt]
				dico_f[idt] = linker
	elif basename(file) == 'STDOUT_pos' :
		for idt in dico_pos.keys() :
			for linker in dico.values() :
				dico_f[idt] = dico_pos[idt]
				dico_f[idt] = linker
	'''

	'''
	for tem, dic in dico_f.items() :
		for idt, liste in dic.items() :
			l = []
			for i in range(len(liste)-2) :
				if type(liste[i]) == type(0) :
					if liste[i+2] == liste[i]+1 :
						print(liste)
						if int(liste[i-1]) > int(liste[i+1]) : 
							pass
							l.append(liste[i-1])
							l.append(liste[i])
							liste.remove(liste[i+1], liste[i+2])
						elif int(liste[i+1]) > int(liste[i-1]) :
							pass
							l.append(liste[i+1])
							l.append(liste[i+2])
							liste.remove(liste[i-1], liste[i])					

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
			if not line.startswith('#') :
				elem = line.split(" ")
				idt = elem[0]
				adressage = elem[1]
				if adressage == 'chlo' or adressage == 'mito' :
					dico[idt] = ""
					dico[idt] += str(adressage)

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
				dico[idt] = {'l_seq', 'l_rep', 'pos', 'aa_prop', 'seq'}
			if '- ' in line :
				ligne = line.strip()
				ligne = ligne.split()
				debut = ligne[0].split('-')
				pos = [debut[0], ligne[1]]
				repet.append(pos)


			dico[idt] = {'pos' : repet, 'l_rep' : l, 'aa_prop' : aa_prop, 'seq' : ""}

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
	with open(path+"/tmhmm_filtred/New_Proteom_All.txt", "r") as filin :
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
		if basename(file) == 'short_output_pos' :
			dico_trgp2['pos'] = {}
			dico_trgp2['pos'] = targetp2(file)


	dico_wlf = {}
	for file in file_wlf :
		if basename(file) == 'output_neg.wolfpsort' :
			dico_wlf['neg'] = {}
			dico_wlf['neg'] = wolfpsort(file)
		if basename(file) == 'output_pos.wolfpsort' :
			dico_wlf['pos'] = {}
			dico_wlf['pos'] = wolfpsort(file)


	dico_ard2 = {}
	for file in file_ard2 :
		if basename(file) == 'STDOUT_neg' :
			dico_ard2['neg'] = {}
			dico_ard2['neg'] = ard2(file)
		if basename(file) == 'STDOUT_pos' :
			dico_ard2['pos'] = {}
			dico_ard2['pos'] = ard2(file)
	

	dico_loca = {}
	for file in file_loca : 
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_localizer' :
			dico_loca['neg'] = {}
			dico_loca['neg'] = localizer(file)
		if basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_localizer' :
			dico_loca['pos'] = {}
			dico_loca['pos'] = localizer(file)

	
	dico_dploc = {}
	for file in file_dploc :
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_deeploc.txt' :
			dico_dploc['neg'] = {}
			dico_dploc['neg'] = deeploc(file)
		if basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_deeploc.txt' :
			dico_dploc['pos'] = {}
			dico_dploc['pos'] = deeploc(file)	


	dico_radar = {}
	for file in file_radar :
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_radar' :
			dico_radar['neg'] = {}
			dico_radar['neg'] = radar(file)
		if basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_radar' :
			dico_radar['pos'] = {}
			dico_radar['pos'] = radar(file)


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


	for tem, dic in dico_ard2.items() :
		idt_l = []
		for idt, val in dic.items() :
			idt_l.append(idt)
		if tem == 'neg' :
			for idt in df.index :
				idt = str(idt)
				if idt in idt_l :
					df.loc[idt, 'type'] = 1

	for tem, dic in dico_trgp2.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'trp2'] = res

	for tem, dic in dico_wlf.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'wolfpsort'] = res

	for tem, dic in dico_ard2.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'ard2'] = [res]	

	for tem, dic in dico_loca.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'localizer'] = res

	for tem, dic in dico_dploc.items() :
		for idt, res in dic.items() :
			idt = ">"+idt
			if idt in idt_all : 
				df.loc[idt, 'deeploc'] = [res]

	for tem, dic in dico_radar.items() :
		for idt, res in dic.items() :
			if idt in idt_all : 
				df.loc[idt, 'radar'] = [res]


	return df


def Tsne(dataframe) :
	pass

	#for col in dataframe.columns :
	#	print(col)

	#print(dataframe['ard2'])





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
	tsne = Tsne(final_results)
