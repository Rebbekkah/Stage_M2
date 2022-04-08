import pandas as pd
import numpy as np
import math as m
import os
import glob
import seaborn as sns
from os.path import basename
import matplotlib.pyplot as plt



def which_proteom() :

	os.chdir(path_Chlamy_arabi+'Predictions/')

	file = path_Chlamy_arabi+'Predictions/prot_alpha.txt'
	print(file)

	proteoms = glob.glob(path_Chlamy_arabi+'TMHMM/prote/*.txt')
	proteoms.sort()
	print(proteoms)

	alpha = []
	with open(file, 'r') as filin :
		for line in filin :
			line = line.strip()
			alpha.append(line)
	print(len(alpha))


	dico = {}
	for p in proteoms :
		print(p, type(p))
		dico[basename(p)] = {}
		dico[basename(p)] = read_proteom(p)
	print(dico.keys(), len(dico.keys()))


	alpha_Chlamy = []
	alpha_Arabi = []
	for a in alpha :
		for org, dic in dico.items() :
			if a in dic.keys() :
				if org == 'New_Proteom_proteome_Arabidopsis_thaliana.faa.txt' :
					alpha_Arabi.append(a)
				elif org == 'New_Proteom_proteome_Chlamydomonas.fa.txt' :
					alpha_Chlamy.append(a)
	print(alpha_Arabi[:100])
	print(len(alpha_Arabi))
	print("----------------------")
	print(alpha_Chlamy[:100])
	print(len(alpha_Chlamy))

	with open('alpha_Arabi.txt', 'w') as filout :
		for a in alpha_Arabi :
			filout.write(a+'\n')

	with open('alpha_Chlamy.txt', 'w') as filout :
		for a in alpha_Chlamy :
			filout.write(a+'\n')

	#return alpha


def read_proteom(file) :

	dico = {}
	with open(file, 'r') as filin :
		for line in filin :
			if line.startswith('>') :
				idt = line.split(' ')[0]
				dico[idt] = ""
			else :
				dico[idt] += line.strip()

	return dico


def is_pos_neg() :

	dico = {}
	proteoms = glob.glob(path_pos_neg+'*.fasta_line')
	proteoms.sort()
	print(proteoms)

	for p in proteoms :
		dico[basename(p)] = {}
		dico[basename(p)] = read_proteom(p)
	print(dico.keys(), len(dico.keys()))

	is_pos = []
	is_neg = []

	#AT != NP !!!!!!!!! demander à céline les bons protéomes


def select_imp(file) :

	os.chdir(path_Chlamy_arabi+'Predictions/')

	df = pd.read_csv(path_Chlamy_arabi+file, sep = '\t')
	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']

	print(df)


	dico = {}
	for index, elem in enumerate(df['importance']) :
		if elem > 0.02 :
			dico[df.loc[index, 'variable']] = 0
			dico[df.loc[index, 'variable']] += float(elem)
	print(dico.keys(), len(dico.keys()))

	with open('Most_imp_desc.txt', 'w') as filout :
		for cle, value in dico.items() :
			filout.write(cle+'\t'+str(value)+'\n')


	return dico



def get_idt(file) :

	idt = []
	with open(file, 'r') as filin :
		for line in filin :
			idt.append(line.strip())
	print("idt ", len(idt))

	return idt


def comp_res_Celine(path) :

	print("----------- MÉTHODE 1&2 -----------")
	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/')


	file = path_Chlamy_arabi+'Predictions/prot_alpha.txt'
	print(file, len(file), type(file))

	res = glob.glob(path)
	res.sort()
	print(res, len(res))

	dico = {}
	for r in res :
		dico[basename(r)] = {}
		dico[basename(r)] = read_proteom(r)


	#alpha = []
	#with open(file, 'r') as filin :
	#	for line in filin :
	#		alpha.append(line.strip())
	#print("ALPHA ", len(alpha))

	alpha = get_idt(file)
	print("ALPHA ", len(alpha))

	'''
	yes = []
	no = []
	new_pred = []
	dickeys = []
	keys = []
	for org, dic in dico.items() :
		dickeys.append(list(dic.keys()))
		for idt in dic.keys() :
			if idt in alpha :
				yes.append(idt)
			if idt not in alpha :
				no.append(idt)
	for elem in dickeys :
		keys = keys+elem
	for a in alpha :
		if a not in keys :
			if a not in new_pred :
				new_pred.append(a)
	'''


	yes, no, new_pred = is_in(dico, alpha)

	print("PRED & PRED ", len(yes))
	print("CEL PRED MAIS PAS MOI ", len(no))
	print("NEW PRED PAR MOI ", len(new_pred))

	with open('pred_pred.txt', 'w') as filout :
		for ident in yes :
			filout.write(ident+'\n')
	with open('Cel_pred_moi_non.txt', 'w') as filout :
		for ident in no :
			filout.write(ident+'\n')
	with open('new_pred.txt', 'w') as filout :
		for ident in new_pred :
			filout.write(ident+'\n')

	with open('new_pred_Chlamy.txt', 'w') as filout_1 :
		with open('new_pred_Arabidopsis.txt', 'w') as filout_2 :
			for ident in new_pred :
				if ident.startswith('>Cre') :
					filout_1.write(ident+'\n')
				else :
					filout_2.write(ident+'\n')


def comp_methode_2(path) :

	print("----------- MÉTHODE 2 -----------")

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/comp_M2/')


	file = path_Chlamy_arabi+'Predictions/prot_alpha.txt'
	print(file, len(file), type(file))

	res = glob.glob(path)
	res.sort()
	print(res, len(res))

	dico = {}
	for r in res :
		dico[basename(r)] = {}
		dico[basename(r)] = read_proteom(r)

	print(dico.keys())
	for ident in dico.values() :
		print(len(ident))
	#print(dico.values())

	alpha = get_idt(file)
	print("len aplha : ", len(alpha))

	pred_pred, Cel_pred, new_pred_m2 = is_in(dico, alpha) 

	print("PRED & PRED ", len(pred_pred))
	print("CEL PRED MAIS PAS MOI ", len(Cel_pred))
	print("NEW PRED PAR MOI ", len(new_pred_m2))


	with open('pred_pred_Chlamy_m2.txt', 'w') as filout_1 :
		with open('pred_pred_Arabi_m2.txt', 'w') as filout_2 :
			for ident in pred_pred :
				if ident.startswith('>Cre') :
					filout_1.write(ident+'\n')
				else :
					filout_2.write(ident+'\n')

	with open('Cel_pred_moi_non_Chlamy_m2.txt', 'w') as filout_1 :
		with open('Cel_pred_moi_non_Arabi_m2.txt', 'w') as filout_2 :
			for ident in Cel_pred :
				if ident.startswith('>Cre') :
					filout_1.write(ident+'\n')
				else :
					filout_2.write(ident+'\n')

	with open('new_pred_vs_m2.txt', 'w') as filout :
		for ident in new_pred_m2 :
			filout.write(ident+'\n')

	with open('new_pred_Chlamy_vs_m2.txt', 'w') as filout_1 :
		with open('new_pred_Arabidopsis_vs_m2.txt', 'w') as filout_2 :
			for ident in new_pred_m2 :
				if ident.startswith('>Cre') :
					filout_1.write(ident+'\n')
				else :
					filout_2.write(ident+'\n')


def is_in(dico, list_idt) :

	yes = []
	no = []
	new_pred = []
	dickeys = []
	keys = []
	for org, dic in dico.items() :
		dickeys.append(list(dic.keys()))
		for idt in dic.keys() :
			if idt in list_idt :
				yes.append(idt)
			if idt not in list_idt :
				no.append(idt)
	for elem in dickeys :
		keys = keys+elem
	for a in list_idt :
		if a not in keys :
			if a not in new_pred :
				new_pred.append(a)

	return yes, no, new_pred


def sep_alpha() :

	print("---------- SEP ALPHA ----------")
	os.chdir(path_Chlamy_arabi+'Predictions/')

	file = path_Chlamy_arabi+'Predictions/prot_alpha.txt'

	alpha = get_idt(file)

	with open('alpha_Arabi.txt', 'w') as filout_1 :
		with open('alpha_Chlamy.txt', 'w') as filout_2 :
			for a in alpha :
				if a.startswith('>Cre') :
					filout_2.write(a+'\n')
				else :
					filout_1.write(a+'\n')



def proteom_alpha() :

	os.chdir(path_Chlamy_arabi+'Predictions/')
	
	all_proteom = Proteom_all(path_Chlamy_arabi)

	alpha_Chlamy = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/alpha_Chlamy.txt')
	alpha_Arabi = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/alpha_Arabi.txt')

	list_alpha = [alpha_Arabi, alpha_Chlamy]


	proteom = {}
	for l in list_alpha :
		dic = {}
		if l is alpha_Arabi :
			print("Arabi")
			for prot in l :
				dic[prot] = all_proteom[prot]
			proteom['Arabidopsis'] = dic
		if l is alpha_Chlamy :
			print("Chlamy")
			for prot in l :
				dic[prot] = all_proteom[prot]
			proteom['Chlamydomonas'] = dic


	with open('Proteom_alpha_Arabidopsis.txt', 'w') as filout_1 :
		with open('Proteom_alpha_Chlamydomonas.txt', 'w') as filout_2 :
				for org, dic in proteom.items() :
					if org == 'Arabidopsis' :
						for prot, seq in dic.items() :
							filout_1.write(prot+'\n'+seq+'\n')
					if org == 'Chlamydomonas' :
						for prot, seq in dic.items() :
							filout_2.write(prot+'\n'+seq+'\n')




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
	with open(path+'New_Proteom_All.txt', "r") as filin :
		for line in filin : 
			if line.startswith('>') :
				idt = line.strip()
				dico[idt] = ""
			else :
				dico[idt] += line.strip()

	return dico


def minus_log_evalue(pattern) :


	os.chdir(path_Chlamy_arabi+'Predictions/Pour_Celine_comp/db_Arabi/')

	files = glob.glob(path_Chlamy_arabi+pattern)
	files.sort()
	print(files)

	evalue = []
	l = []
	with open(files[0], 'r') as filin :
		for line in filin :
			l.append(line)
			res = line.split('\t')[2]
			evalue.append(float(res))
	print(len(evalue))

	print(evalue[:10], len(evalue))
	mlog = []
	for ev in evalue :
		if ev == 0.0 :
			ev = 10**-300
		res = -m.log(ev)
		mlog.append(res)
	print(mlog[:10], len(mlog))


	print(l[:10], len(l))
	i = 0
	new = []
	for prot in l :
		#print(prot.split('\t')[2])
		#prot.split('\t')[2] = mlog[i]
		#print(prot)
		new_line = prot.split('\t')
		#print(new_line)
		new_line[2] = str(mlog[i])
		#print(new_line)
		new_line = '\t'.join(new_line)
		#print(new_line)
		new.append(new_line)
		i += 1
	print(new[:10], len(new))

	with open('for_cytoscape_Arabi.csv', 'w') as filout :
		for line in new :
			filout.write(line)


def correspondance_acc(file) :


	os.chdir(path_Chlamy_arabi+'Predictions/')

	df = pd.read_csv(path_Chlamy_arabi+file, sep = '\t')

	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']

	acc_list = []
	for col in df :
		if col.startswith('acc') :
			acc_list.append(col)

	real_acc = []
	with open(path_to_script+'acc/colnames_acc.txt', 'r') as filin :
		for line in filin :
			real_acc.append(line.strip())

	acc_df = pd.DataFrame(acc_list, columns = ['acc_df'])
	df_real = pd.DataFrame(real_acc, columns = ['real'])
	df_acc = pd.concat([acc_df, df_real], axis = 1)
	
	df_acc.to_csv('acc.csv', sep = '\t', header = True, index = True)



if __name__ == '__main__' :

	#path = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/neg_pos/'
	#path_prote = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/other/'
	path_Chlamy_arabi = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/RF/Chlamy_Arabi/results/"
	path_pos_neg = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/"
	path_method_Cel = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/methode_1_2_Celine/'
	path_to_script = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/"

	os.chdir(path_Chlamy_arabi)

	#which_proteom()
	#is_pos_neg()
	#dico_imp = select_imp('Importance_desc.csv')
	#comp_res_Celine(path_method_Cel+'*/*')
	#comp_methode_2(path_method_Cel+'M2_*/*')
	#sep_alpha()
	#proteom_alpha()
	#minus_log_evalue('Predictions/Pour_celine_comp/db_*/*.out')
	correspondance_acc('Predictions/dataframe_all.csv')


