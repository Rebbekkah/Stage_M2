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
	print(basename(file))
	dico = {}
	#idt = []
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
					#dico_linker[aa] = []
					#dico_linker[aa].extend([proba, pos])
					#dico_linker[aa] = [proba, pos]
					#for cle, val in dico_linker.items() :
					if aa in dico_linker.keys() :
						#dico_linker[aa].append([proba, pos])
						dico_linker[aa].extend([proba, pos])
					else :
						dico_linker[aa] = [proba, pos]
			dico[idt] = dico_linker


	#### changer les cles du dico en identifiant protéique en comparant les seq
	#### + récupérer les positions avec les lignes et faire attention aux idt 
	#### qui prennent des lignes

	print(dico)
	return dico


def wolfpsort(file) :

	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if not line.startswith('#') :
				elem = line.split(" ")
				#print(elem)
				idt = elem[0]
				#print(idt)
				adressage = elem[1]
				#print(adressage)
				#print(adressage, type(adressage))
				if adressage == 'chlo' or adressage == 'mito' :
					dico[idt] = ""
					dico[idt] += str(adressage)

	#print(dico)
	return dico


def targetp2(file) :

	dico = {}
	with open(file, "r") as filin :
		for line in filin :
			if not line.startswith('#') :
				idt = line.split()[0]
				dico[idt] = ""
				dico[idt] += line.split()[1]

	return dico


def dataFrame() :

	df_pos = pd.DataFrame()
	df_neg = pd.DataFrame()

	file_ard2 = glob.glob(path_ard2+"STDOUT_"+"*")
	print(file_ard2)
	file_wlf = glob.glob(path_wpsort+"*.wolfpsort")
	#print(file_wlf)
	file_trgp2 = glob.glob(path_trgp2+"short_output_"+"*")
	#print(file_trgp2)

	dico_trgp2 = {}
	for file in file_trgp2 :
		if basename(file) == 'short_output_neg' :
			dico_trgp2['neg'] = {}
		dico_trgp2['neg'] = targetp2(file)
		if basename(file) == 'short_output_pos' :
			dico_trgp2['pos'] = {}
		dico_trgp2['pos'] = targetp2(file)
	#print("------------------")
	#print(dico_trgp2, len(dico_trgp2))
	
	#wolfpsort(file_wlf[0])

	dico_wlf = {}
	for file in file_wlf :
		if basename(file) == 'output_neg.wolfpsort' :
			dico_wlf['neg'] = {}
		dico_wlf['neg'] = wolfpsort(file)
		if basename(file) == 'output_pos.wolfpsort' :
			dico_wlf['pos'] = {}
		dico_wlf['pos'] = wolfpsort(file)

	#print(dico_wlf)

	ard2(file_ard2[0])


	#return dico_trgp2


if __name__ == '__main__' :

	#path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests/proteome_diatom.tmhmm"
	path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests_small/"
	path_small_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests_small/"
	path_output = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/"
	to_script = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script"
	list_of_aa = ['M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H']

	path_proteom = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/"

	#read = reading(path)



	os.chdir(path_output)

	# TMHMM
	#proteins = listing(path_proteom, '*.tmhmm')
	#new_proteom = proteome_maker(proteins, path_proteom)
	#separateur = sep(new_proteom, path_proteom)
	#separateur = sep(path_proteom)

	# ARD2
	path_ard2 = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/ard2_outputs/"

	# WOLFPSORT
	path_wpsort = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/wolfpsort_output/"

	# TARGETP2
	path_trgp2 = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/targetp2_outputs/"


	data_final = dataFrame()


