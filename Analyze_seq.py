import pandas as pd
import os
from os.path import basename
import glob


def reading(fichier) :
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
		#print(start)
		#print(end)

		if start > 68 or end > 68 :
			dico_h[cle] = {}
			dico_h[cle] = {'start' : start, 'end' : end}
		elif end < 68 :
			#print(cle)
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
	#os.chdir(to_script)

	return id_to_keep


def listing(path, pattern) :

	fich = glob.glob(path+pattern)
	print(fich)
	prot = []
	
	for f in fich :
		print("----------", basename(f), "----------")
		prot.append(reading(f))

	return prot



def proteome_maker(ToKeep, path_proteom) :
	print("-----------------proteom maker-----------------")
	#print(ToKeep)
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
			#print(ident)
			if ident in dico.keys() :
				dico2[ident] = ""
				dico2[ident] = dico[ident]
				

	#print(len(dico2.keys()))
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
	
	#print(len(ToKeep[0])+len(ToKeep[1]))
	#print(idt, len(idt))
	#print(seq, len(seq))
	#print(dico)
	return dico2

def sep(path_proteom) :
	print("-----------------separateur-----------------")
	proteom = glob.glob(path_proteom+"*.fasta_line")
	proteom = proteom[::-1]
	file = glob.glob(path_proteom+"*.tmhmm")
	print(proteom, "\n", file)

	idt_ok = []
	for f in file :
		idt_ok.append(reading(f))
	#for i in idt_ok :
	#	print("LEN IDT OK ", len(i))

	os.chdir(path_output)
	idt_all = []
	#print("LEN DICO", len(dico))	
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
			#print("LEN IDT", len(idt_ok[i]))
			#print(idt)
			#break
			idt = ">"+idt
			if idt in dico_all.keys() :
				#print("oui")
				dico2[idt] = ""
				dico2[idt] += dico_all[idt]
		print("LEN DICO2", len(dico2))
		with open("New_Proteom_"+basename(p)+".txt", "w") as filout :
			for idt, seq in dico2.items() :
				filout.write(idt+"\n"+seq+"\n")
		i += 1



def ard2(file) :
	dico = {}
	#dico_linker = {} 
	with open(file, "r") as filin :
		for line in filin :
			if line.startswith('>') :
				dico[line] = ""
				dico_linker = {} 
			else :
				#dico_linker = {} 
				#dico[line] += line.split("\t")[0]
				#elem = line.split("\t")[0]
				#aa = elem.split()[0]
				aa = line.split()[0]
				proba = line.split()[1:4]
				#proba = elem.split()[1]
			if proba > 0.10 :
				dico_linker[aa] = ""
				dico_linker[aa] = proba 
			dico[line] += dico_linker


	####changer les cles du dico en identifiant protéique

	return dico


def wolfpsort(file) :
	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if not line.startswith('#') :
				elem = line.split(" ")
				idt = elem[0]
				adressage = elem[1]
			if adressage == 'chlo' or adressage == 'mito' :
				dico[idt] = ""
				dico[idt] += adressage

	return dico




if __name__ == '__main__' :

	#path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests/proteome_diatom.tmhmm"
	path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests_small/"
	path_small_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests_small/"
	path_output = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/"
	to_script = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script"
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']

	path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/"

	#read = reading(path)

	# TMHMM
	proteins = listing(path_proteom, '*.tmhmm')
	new_proteom = proteome_maker(proteins, path_proteom)
	#separateur = sep(new_proteom, path_proteom)
	separateur = sep(path_proteom)

	# ARD2


	# WOLFPSORT








