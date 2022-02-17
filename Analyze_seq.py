import pandas as pd
import sys
import os
from os.path import basename
import glob
#import codecs
#import unicodedata
#import unidecode
from operator import itemgetter


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

	fich = glob.glob(path+pattern)
	print(fich)
	prot = []
	
	for f in fich :
		print("----------", basename(f), "----------")
		prot.append(reading(f))

	return prot



def proteome_maker(ToKeep, path_proteom) :
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
	#print(basename(file))
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
	#print(proteoms)

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


	return dico_f


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
					dico[idt] += str(adressage)

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

	#print(basename(file))
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

	idt_list = []
	#repet = []
	#l = []
	dico = {}
	with open(file, "r", encoding = 'ascii', errors = 'ignore') as filin :
		for line in filin :
			if line.startswith('>') or '>' in line :
				repet = []
				l = []
				aa_prop = []
				idt_list.append(line.strip())
				idt = line.strip()
				dico[idt] = {'l_seq', 'l_rep', 'pos', 'aa_prop'}
				#print("----------------", idt)
			if '- ' in line :
				ligne = line.strip()
				ligne = ligne.split()
				debut = ligne[0].split('-')
				pos = [debut[0], ligne[1]]
				#longr = int(ligne[1]) - int(debut[0]) ####### positions plus bonnes après avoir modif les hevauchement \
													  ####### à revoir 
				#l.append(longr)
				repet.append(pos)


			dico[idt] = {'pos' : repet, 'l_rep' : l, 'aa_prop' : aa_prop}

		k = 0
		for idt, dic in dico.items() :
			print(idt)
		#for dic in dico.values() :
			to_remove = []
			rep = []
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
					print(rep2, len(rep2))
					print(rep2[i], rep2[i+1])
					to_remove.append(rep2[i+1])
			
				if rep2[i][0] > rep2[i+1][0] and rep2[i][1] < rep2[i+1][0] :
					nb += 1
					print(rep2)
					print(rep2[i], rep2[i+1])
					to_remove.append(rep2[i])

				if rep2[i][0] < rep2[i+1][1] and rep2[i][1] > rep2[i+1][1] \
				and rep2[i+1][0] < rep2[i][0] :
					nb += 1
					print(rep2)
					print(rep2[i], rep2[i+1])

				if rep2[i][1] < rep2[i+1][1] and rep2[i+1][0] > rep2[i+1][0] \
				and rep2[i][1] > rep2[i+1][0] :
					nb += 1
					print(rep2)
					print(rep2[i], rep2[i+1])

			if to_remove :
				print("REPT À SUPP ----> ", to_remove, len(to_remove))
			if nb != 0 :
				print("nb de chevauchement dans la séquence :", nb)
				k += 1
				for elem in to_remove :
					print("REMOVING")
					print(elem)
					rep2.remove(elem)
					print(rep2, len(rep2))
				print("\n============================\n")

			dic['pos'] = rep2

		print('nb de seq avec chevauchement :', k)
		

		#print(dico)


	print(basename(file))
	

def verif() : 
	
	fich = path_output+"output_Analyzeseq_ToKeep_tem_neg.tmhmm.txt"
	idt_radar = radar(path_pb)
	#print(idt_radar, len(idt_radar))
	#print(fich, len(fich))

	idt_keep = []

	with open(fich, 'r') as filin :
		for line in filin :
			idt_keep.append(line.strip())

	#print(idt_keep, len(idt_keep))

	idt_not = []
	for idt in idt_keep :
		idt = '>'+idt
		if idt not in idt_radar :
			idt_not.append(idt)
	print("-----------", idt_not, len(idt_not))


def Data_Create() :

	#df_pos = pd.DataFrame()
	#df_neg = pd.DataFrame()

	file_ard2 = glob.glob(path_ard2+"STDOUT_"+"*")
	#print(file_ard2)
	file_wlf = glob.glob(path_wpsort+"*.wolfpsort")
	#print(file_wlf)
	file_trgp2 = glob.glob(path_trgp2+"short_output_"+"*")
	#print(file_trgp2)
	file_dploc = glob.glob(path_dploc+"*"+"deeploc"+"*")
	#print(file_dploc)
	file_loca = glob.glob(path_loca+'*'+'localizer')
	#print(file_loca)
	file_radar = glob.glob(path_radar+'*'+'radar')
	#print(file_radar)

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

	#ard2(file_ard2[0])

	dico_ard2 = {}
	for file in file_ard2 :
		if basename(file) == 'STDOUT_neg' :
			dico_ard2['neg'] = {}
			dico_ard2['neg'] = ard2(file)
		if basename(file) == 'STDOUT_pos' :
			dico_ard2['pos'] = {}
			dico_ard2['pos'] = ard2(file)
	
	#print(dico_ard2)

	#localizer(file_loca[0])
	dico_loca = {}
	for file in file_loca : 
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_localizer' :
			dico_loca['neg'] = {}
			dico_loca['neg'] = localizer(file)
		if basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_localizer' :
			dico_loca['pos'] = {}
			dico_loca['pos'] = localizer(file)
	#print(dico_loca)
	#return dico_trgp2

	#deeploc(file_dploc[0])
	
	dico_dploc = {}
	for file in file_dploc :
		if basename(file) == 'New_Proteom_1196_tem_neg.fasta_line.txt_deeploc.txt' :
			dico_dploc['neg'] = {}
			dico_dploc['neg'] = deeploc(file)
		if basename(file) == 'New_Proteom_1081_tem_pos.fasta_line.txt_deeploc.txt' :
			dico_dploc['pos'] = {}
			dico_dploc['pos'] = deeploc(file)	
	#print(dico_dploc)


	radar(file_radar[0])
	#radar(path_pb)


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

	data_final = Data_Create()
	#verif()


