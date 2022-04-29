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


def read_df(path_df) :

	df = pd.read_csv(path_df+'dataframe_all.csv', sep = '\t')
	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']

	return df


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
	print("Number of idt --> ", len(idt))

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


	#os.chdir(path_Chlamy_arabi+'Predictions/Pour_Celine_comp/db_Chlamy')

	files = glob.glob(path_Chlamy_arabi+pattern)
	files.sort()
	print(files)
	lp = ['db_Arabi', 'db_Chlamy']

	i = 0
	for f in files :
		p = lp[i]
		os.chdir(path_Chlamy_arabi+'Predictions/Pour_Celine_comp/'+p+'/')
		i += 1
		evalue = []
		l = []
		with open(f, 'r') as filin :
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
			res = -m.log10(ev)
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


		with open('for_cytoscape_'+p+'.csv', 'w') as filout :
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



def adressage_alpha(file1, file2) :

	df = read_df(path_Chlamy_arabi+'Predictions/')

	#df_2 = df.iloc[:150, :]

	print(df)
	adress = ['wolfpsort', 'deeploc', 'trp2', 'localizer']

	dico = {'Arabi' : get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/'+file1), \
	'Chlamy' : get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/'+file2)}

	#print(dico)

	df_adr_Chl = pd.DataFrame(index = dico['Chlamy'], columns = [adress])
	#print(df_adr_Chl)
	df_adr_Arabi = pd.DataFrame(index = dico['Arabi'], columns = [adress])
	df_adr = pd.DataFrame(index = dico['Chlamy']+dico['Arabi'], columns = [adress])
	#print("-------------------\n", df_adr)
	

	for org, lidt in dico.items() :
		for idt in lidt :
			for software in adress :
				#df_adr.loc[idt, software] = df_adr.loc[idt, software]
				df_adr.loc[idt, software] = df.loc[idt, software]
				if org == 'Chlamy' :
					df_adr_Chl.loc[idt, software] = df.loc[idt, software]
				elif org == 'Arabi' :
					df_adr_Arabi.loc[idt, software] = df.loc[idt, software]

	print("DF ADR\n", df_adr)
	#print(df_adr_Chl)
	#print(df_adr_Arabi)

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	df_adr.to_csv('df_adr.csv', sep = '\t', header = True, index = True)
	df_adr_Chl.to_csv('df_adr_Chl.csv', sep = '\t', header = True, index = True)
	df_adr_Arabi.to_csv('df_adr_Arabi.csv', sep = '\t', header = True, index = True)

	ldf = [df_adr, df_adr_Chl, df_adr_Arabi]

	idt_adr = []
	idt_adr_Chl = []
	idt_adr_Arabi = []
	truc = []
	#print("--------------")
	for d in ldf :
		#print("--------------")
		#print(d)
		d_mask = d['deeploc'] >= 0.25
		filtered_df = d['deeploc'][d_mask].dropna().index
		new_df = d.loc[filtered_df]
		d_mask = new_df['wolfpsort'] >= 0.5
		new_df = new_df.loc[filtered_df]
		#print("NEW DF ---------------\n", new_df)
		#print(len(d))
		truc.append(new_df)

	print(truc)


	index = []
	for t in truc :
		index.append(list(t.index))

	print("----------------------")
	l = ['all', 'Chl', 'Arabi']

	i = 0
	for ind in index :
		with open('prot_adress_'+l[i]+'.txt', 'w') as filout :
			for idt in ind :
				filout.write(idt+"\n")
		i += 1


def is_ppr_opr(file) :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_Celine_comp/')

	df = pd.read_csv(file, sep = '\t')
	df = df.set_index(df['protein_id'], inplace = False)
	del df['protein_id']

	#for col in df :
	#	print(df[col])

	lidt_opr = []
	k = 0
	for index, elem in enumerate(df['source']) :
		k += 1
		print(elem)
		if elem.startswith('Cre') :
			lidt_opr.append(elem.split()[0])
		else :
			lelem = []
			print("-----------")
			print(elem)
			print(elem.split())
			lelem = elem.split()
			for el in lelem :
				if el.startswith('Cre') :
					lidt_opr.append(el)
	#print(lidt_opr, len(lidt_opr))

	for i in range(len(lidt_opr)) :
		if '|' in lidt_opr[i] :
			lidt_opr[i] = lidt_opr[i].split('|')[0]
	print(lidt_opr, len(lidt_opr))

	#print(df.loc['Chlre_OPR102', 'source'])
	lidt_alpha = []
	with open(path_Chlamy_arabi+'Predictions/Pour_celine_comp/alpha_Chlamy.txt') as filin :
		for line in filin :
			lidt_alpha.append(line.strip())

	print("nb of opr in file :", k)
	print(len(lidt_alpha))

	alpha_opr = []
	for idt in lidt_opr :
		idt = '>'+idt
		if idt in lidt_alpha :
			#print(idt)
			alpha_opr.append(idt)

	with open('opr_in_alpha_pred_Chlamy.txt', 'w') as filout :
		for idt in alpha_opr :
			filout.write(idt+'\n')
	

	
	fich = path_to_script+'Celine/methode_1_2_Celine/M1_OPR_Chlre/cat_loops+blastp_motifs.fasta_PROTEINS_2rep_chlre_M1OPR_CHL'
	Cel_opr = read_proteom(fich)

	Cel_opr_in = []
	print(len(Cel_opr.keys()))
	for idt in lidt_opr :
		idt = '>'+idt
		if idt in Cel_opr.keys() :
			#print(idt)
			Cel_opr_in.append(idt)


	print("MOI :", alpha_opr)
	print("CÉLINE :", Cel_opr_in)

	#for col in df :
	#	print(df[col])



def comp_Hedi(file) :

	os.chdir(path_Chlamy_arabi+'Predictions/Comp_Hedi/')

	alpha_Hedi = []
	non_alpha_Hedi = []
	with open(path_Chlamy_arabi+file, 'r') as filin :
		for line in filin :
			line = line.split(' ')
			line[-1] = line[-1].strip()
			if line[-1] == '1' :
				alpha_Hedi.append('>'+line[1])
			else : 
				non_alpha_Hedi.append('>'+line[1])
	print("prot alpha_Hedi", len(alpha_Hedi))
	print("prot non alpha_Hedi", len(non_alpha_Hedi))

	alpha = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/prot_alpha.txt')

	same = []
	not_same = []
	his_pred = []
	for a in alpha :
		if a in alpha_Hedi :
			same.append(a)
		else :
			not_same.append(a)

	for a in alpha_Hedi :
		if a not in alpha :
			his_pred.append(a)

	with open('same_pred_1746.txt', 'w') as filout1 :
		with open('my_pred_1025.txt', 'w') as filout2 :
			with open('his_pred_6853.txt', 'w') as filout3 :
				for prot in same :
					filout1.write(prot+'\n')
				for prot in not_same :
					filout2.write(prot+'\n')
				for prot in his_pred :
					filout3.write(prot+'\n')


def right_proteom_opr(file) :

	os.chdir(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/')

	dico = read_blast(file)

	print(len(dico.keys()))

	proteom_opr = read_proteom(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/chlamydomonas_opr.fasta')
	proteom_Chlamy = read_proteom(path_Chlamy_arabi+'Predictions/Proteom_alpha_Chlamydomonas.txt')
	#print(len(proteom))

	new_prot = {}
	for idt, seq in proteom_opr.items() :
		for old, new in dico.items() :
			old = '>'+old
			new = '>'+new
			if old == idt :
				new_prot[new] = seq

	print(len(new_prot))


	with open('New_proteom_OPR_alpha_pred_Chlamy.fasta', 'w') as filout1 :
		with open('prot_alpha_opr_Chlamy.txt', 'w') as filout2 :
			for idt, seq in new_prot.items() :
				filout1.write(idt+'\n'+seq+'\n')
			for idt in new_prot.keys() :
				filout2.write(idt+'\n')


	
	'''
	#print(proteom_Chlamy, len(proteom_Chlamy))
	proteom_C = {}
	key = list(proteom_Chlamy.keys())
	i = 0
	for k in key :
		k = k.strip()
		key[i] = k
		i += 1
	print(len(key))
	for k in key :
		for idt, seq in proteom_Chlamy.items() :
			proteom_C[k] = seq
	#print(proteom_C, len(proteom_C))

	opr_not_pred = []
	z = 0
	print(new_prot.keys(), len(new_prot))
	for idt in new_prot.keys() :
		if idt not in proteom_C.keys() :
			opr_not_pred.append(idt)
		else :
			z += 1
	print(opr_not_pred, len(opr_not_pred))
	print(z)
	'''

def read_blast(file) :

	dico = {}
	with open(file, 'r') as filin :
		for line in filin :
			line = line.split('\t')
			if line[2] == '0.0' :
				if line[0] not in dico.keys() :
					dico[line[0]] = line[1]

	return dico


def comp_pos_neg(Arabi_vs_neg, Arabi_vs_pos, Chlamy_vs_neg, Chlamy_vs_pos) :


	os.chdir(path_Chlamy_arabi+'Predictions/Res_blast/')

	dico_Arabi = {'positive' : read_blast(Arabi_vs_pos), 'negative' : read_blast(Arabi_vs_neg)}
	dico_Chlamy = {'positive' : read_blast(Chlamy_vs_pos), 'negative' : read_blast(Chlamy_vs_neg)}

	print(dico_Chlamy, len(dico_Chlamy))
	#print(dico_Arabi, len(dico_Arabi))

	ldico = [dico_Arabi, dico_Chlamy]
	p = 0
	n = 0
	for dico in ldico :
		for tpe, dic_idt in dico.items() :
			if dico is dico_Arabi :
				print("Arabi -------")
			else :
				print("Chlamy -------")
			print(tpe, '\n', len(dic_idt))
			#print(len(dic_idt))
			print("----------")
			if dico is dico_Arabi :
				with open('Arabi_alpha_pred_found_in_'+str(tpe)+'.txt', 'w') as filout :
					for old, new in dic_idt.items() :
						filout.write(new+'\n')
			else : 
				with open('Chlamy_alpha_pred_found_in'+str(tpe)+'.txt', 'w') as filout :
					for old, new in dic_idt.items() :
						if old.startswith('Cre') :
							if tpe == 'positive' :
								p += 1
							else : 
								n += 1
							filout.write(new+'\n')
	print(p, n)




def for_cytoscape() :

	os.chdir(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/')

	opr = []
	with open(path_Chlamy_arabi+'Predictions/Res_blast/Chlamy_alpha_pred_found_in_postitive.txt', 'r') as filin :
		for line in filin :
			opr.append('>'+line.strip())
	#print(opr, len(opr))

	alpha_chl = []
	with open(path_Chlamy_arabi+'Predictions/Pour_celine_comp/alpha_Chlamy.txt', 'r') as filin :
		for line in filin :
			alpha_chl.append(line.strip())
	#print(opr)
	#print(alpha_chl)

	with open('for_cytoscape_col_OPR_3.txt', 'w') as filout :
		for idt in alpha_chl :
			if idt in opr :
				filout.write(idt+'\t'+'Yes'+'\n')
			else :
				filout.write(idt+'\t'+'No'+'\n')


def for_cytoscape_2() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	new = []
	with open(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/new_new_arabi.txt', 'r') as filin :
		for line in filin :
			new.append(line.strip())
	#print(new, len(new))

	alpha_chl = []
	with open(path_Chlamy_arabi+'Predictions/Pour_celine_comp/alpha_Arabi.txt', 'r') as filin :
		for line in filin :
			alpha_chl.append(line.strip())
	#print(new)
	#print(alpha_chl)

	with open('for_cytoscape_col_new_arabi.txt', 'w') as filout :
		for idt in alpha_chl :
			if idt in new :
				filout.write(idt+'\t'+'New'+'\n')
			else :
				filout.write(idt+'\t'+'Else'+'\n')





def which_opr(file) :


	os.chdir(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/')
	yes = []
	with open(file, 'r') as filin :
		for line in filin :
			idt = '>'+line.split()[0]
			if line.split()[1].startswith('Yes') :
				yes.append(idt)

	with open('Chlamy_OPR_idt.txt', 'w') as filout :
		for idt in yes :
			filout.write(idt+'\n')



def adressage_alpha_deeploc(path1, path_deeploc, path_wolfpsort) :

	'''
	files_adr = glob.glob(path1+'*.csv')
	files_adr.sort()
	print(files_adr, len(files_adr))

	ldf_adr = []
	idt = []
	for f in files_adr :
		df = pd.read_csv(f, sep = '\t')
		df = df.set_index(df['Unnamed: 0'], inplace = False)
		del df['Unnamed: 0']
		idt.append(list(df.index))
		ldf_adr.append(df)
	#print(ldf_adr)
	#print(idt)
	'''

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	df = pd.read_csv(path1+'df_adr.csv', sep = '\t')
	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']
	#print(df)
	idt = list(df.index)
	#print(idt)
	#print(len(idt))

	files_dploc = glob.glob(path_deeploc+'*DEEPLOC.txt')
	files_dploc.sort()
	print(files_dploc, len(files_dploc))

	files_wlf = glob.glob(path_wolfpsort+'*.wolfpsort')
	files_wlf.sort()
	print(files_wlf, len(files_wlf))

	#idt_dploc = []
	#idt_wlf = []

	#k = 0
	print("----------")

	#for f in files_dploc : 
	k = 0
	z = 0
	#print("-----------------", basename(f), "-----------------")
	idt_dploc = []
	idt_wlf = []
	df_ = pd.read_csv(files_dploc[0], sep = '\t')
	df_ = df_.set_index('>'+df_['ID'], inplace = False)
	del df_['ID']
	#print(list(df_.index))
	print(df_)
	print(df)
	#print(files_dploc)
	#print("--------pppppp-------")
	#print(df_.index)
	#if f is files_dploc[1] :
	#	print(df_.loc['>Cre02.g095550.t1.1', :])
	print("----------tttttt-------")
	for ident in idt :
		#if f is files_dploc[1] :
		#	print("ok1")
		#print(ident)
		if ident in list(df_.index) :
			k += 1
			#print(ident)
			#print("ok")
			#print(basename(f))
			if df_.loc[ident, 'Location'] == 'Mitochondrion' or df_.loc[ident, 'Location'] == 'Plastid' :
				idt_dploc.append(ident)
				#print(df_.loc[ident, 'Location'])
		else :
			z += 1
			#print(ident)
	print("k = ", k)
	print("z = ", z)
	print(len(idt_dploc))
	'''
	if 'Arabidopsis' in f :
		with open('Filtrage_deeploc_alpha_Arabi.txt', 'w') as filout :
			for index in idt_dploc :
				filout.write(index+'\n')
	else :
	'''
	with open('Filtrage_deeploc_alpha_Arabi.txt', 'w') as filout :
		for index in idt_dploc :
			filout.write(index+'\n')

	new_ = df_.loc[idt_dploc, :].dropna().index
	new_df = df_.loc[new_]
	print("NEW DF ---->\n", new_df)

	new_df.to_csv('new_df_Arabi.csv', sep = '\t', header = True, index = True)


'''
col = new_df.columns[1:]
#print(col)
for ident in idt :
	print("-------------ICI")
	machin = list(new_df.loc[ident, col])
	#print(machin)
	#print(max(machin))
	if df.loc[ident, 'wolfpsort'] == max(machin) :
		idt_wlf.append(ident)
	break
	#print("--------")
	#print(max(new_df.loc[ident, 1], new_df.loc[ident, 1]))
	#pass
		#if df.loc[ident, 'wolfpsort'] is max()
'''


def adressage_alpha_wolfpsort(path_file) :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	#files = glob.glob(path_df+'new_df_*.csv')
	#files = glob.glob(path_file+'Filtrage_deeploc_*.txt')
	files = glob.glob(path_file+'*wolfpsort')
	files.sort()
	print(files, len(files))

	alpha = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/prot_alpha.txt')
	#print(alpha, len(alpha))


	ident = []
	for f in files :
		with open(f, 'r') as filin :
			for line in filin :
				idt = '>'+line.split()[0]
				res = line.split()[1]
				if res == 'mito' or res == 'chlo' :
					ident.append(idt)
	#print(ident, len(ident))

	new_alpha = []
	chl = []
	arabi = []
	for a in alpha :
		if a in ident :
			new_alpha.append(a)
			if a.startswith('>Cre') :
				chl.append(a)
			else :
				arabi.append(a)

	print(len(new_alpha))
	print(len(chl))
	print(len(arabi))

	with open('Filtrage_wolpsort_Chlamy.txt', 'w') as filout :
		for new in chl :
			filout.write(new+'\n')
	with open('Filtrage_wolpsort_Arabi.txt', 'w') as filout :
		for new in arabi :
			filout.write(new+'\n')


def adressage_alpha_localizer(path_file) :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')
	alpha = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/prot_alpha.txt')

	files = glob.glob(path_file+'*LOCALIZER')
	files.sort()
	print(files, len(files))

	ident = []
	chl = []
	arabi = []
	new = []
	non = ['Over', '#', 'Identifier', '-', '\n']
	for f in files :
		with open(f, 'r') as filin :
			for line in filin :
				first = line.split(" ")[0]
				if first not in non :
					one = line.split("\t")[0]
					if one[0] not in non :
						idt = line.split("\t")[0]
						idt = idt.split()[0]
						elem = line.split("\t")[1]
						if elem.split()[0] == 'Y' :
							ident.append('>'+idt)
					else :
						break
	print(ident, len(ident))
	for a in alpha :
		if a in ident :
			new.append(a)
			if a.startswith('>Cre') :
				chl.append(a)
			else :
				arabi.append(a)
	#print(arabi, len(arabi))

	with open('Filtrage_localizer_Chlamy.txt', 'w') as filout :
		for p in chl :
			filout.write(p+'\n')
	with open('Filtrage_localizer_Arabi.txt', 'w') as filout :
		for p in arabi :
			filout.write(p+'\n')



def intersection(path_files) :
	
	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	files = glob.glob(path_files+'Filtrage_*.txt')
	files.sort()
	print(files, len(files))

	chl = []
	arabi = []
	for f in files :
		if 'Chlamy.txt' in f :
			chl.append(get_idt(f))
		if 'Arabi.txt' in f :
			arabi.append(get_idt(f))

	inter_chl = []
	inter_arabi = []
	dp_wlf_chl = []
	dp_wlf_arabi = []
	dp_loca_chl = []
	dp_loca_arabi = []
	wlf_loca_chl = []
	wlf_loca_arabi = []


	
	for idt in chl[0] :
		if idt in chl[1] and idt in chl[2] :
			inter_chl.append(idt)
		if idt in chl[1] :
			dp_loca_chl.append(idt)
		if idt in chl[2] :
			dp_wlf_chl.append(idt)
	for idt in chl[1] :
		if idt in chl[2] :
			wlf_loca_chl.append(idt)
	for idt in arabi[0] :
		if idt in arabi[1] and idt in arabi[2] :
			inter_arabi.append(idt)
		if idt in arabi[1] :
			dp_loca_arabi.append(idt)
		if idt in arabi[2] :
			dp_wlf_arabi.append(idt)
	for idt in arabi[1] :
		if idt in arabi[2] :
			wlf_loca_arabi.append(idt)

	liste = [inter_chl, inter_arabi, dp_wlf_chl, dp_wlf_arabi, dp_loca_chl, dp_loca_arabi, wlf_loca_chl, wlf_loca_arabi]
	lstr = ['inter_chl', 'inter_arabi', 'dp_wlf_chl', 'dp_wlf_arabi', 'dp_loca_chl', 'dp_loca_arabi', 'wlf_loca_chl', 'wlf_loca_arabi']

	i = 0
	for l in liste :
		with open('Intersection_'+lstr[i]+'.txt', 'w') as filout :
			for idt in l :
				filout.write(idt+'\n')
		i += 1


def what_is_in_filtrage_deeploc_chl(path_file, file_pos) :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	files = glob.glob(path_file+'Filtrage_deeploc_*')
	files.sort()
	print(files, len(files))

	for f in files :
		if 'Chlamy.txt' in f :
			idt_chl = get_idt(f)
		elif 'Arabi.txt' in f :
			idt_arabi = get_idt(f)
	#print(idt_chl, len(idt_chl))

	dico = read_proteom(file_pos)
	#print(dico, len(dico))
	#print(dico.keys())

	in_pos_chl = []
	not_in_pos_chl = []
	for idt in idt_chl :
		if idt in dico.keys() :
			in_pos_chl.append(idt)
		else :
			not_in_pos_chl.append(idt)
	#print(not_in_pos_chl, len(not_in_pos_chl))
	#print(in_pos_chl, len(in_pos_chl))

	with open('in_pos_chl.txt', 'w') as filout :
		for idt in in_pos_chl :
			filout.write(idt+'\n')
	with open('not_in_pos_chl.txt', 'w') as filout :
		for idt in not_in_pos_chl :
			filout.write(idt+'\n')

	################################## arabi --> idt différents --> blast ?


def comp_new_Cel() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	files = glob.glob(path_method_Cel+'*/*')
	print(files, len(files))
	new = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/not_in_pos_arabi.txt')
	idt = []
	for f in files :
		dico = read_proteom(f)
		idt.append(list(dico.keys()))
	#print(idt, len(idt))

	already = []
	new_new = []
	her = []
	for liste in idt :
		for ident in liste :
			if ident.startswith('>NP') :
				her.append(ident)
	print(her, len(her))

	for ident in new :
		if ident in her :
			already.append(ident)
		else :
			new_new.append(ident)

	print(new_new, len(new_new))

	with open('new_new_arabi.txt', 'w') as filout_1 :
		with open('already_arabi.txt', 'w') as filout_2 :
			for prot in new_new :
				filout_1.write(prot+'\n')
			for prot in already :
				filout_2.write(prot+'\n')



def for_eggNOG() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')

	dico_all = Proteom_all(path_Chlamy_arabi)
	new = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/new_new_arabi.txt')

	dico_new = {}
	for idt in new : 
		if idt in dico_all.keys() :
			dico_new[idt] = ""
			dico_new[idt] = dico_all[idt]

	with open('new_Arabi.fasta', 'w') as filout :
		for idt, seq in dico_new.items() :
			filout.write(idt+'\n'+seq+'\n')


def diff_adr() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/Model_without_adr/')

	with_adr = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/prot_alpha.txt')
	not_adr = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/Model_without_adr/prot_alpha.txt')

	new_pred = []
	common = []
	adr_pred_only = []

	for idt in not_adr :
		if idt in with_adr :
			common.append(idt)
		else :
			new_pred.append(idt)
	for idt in with_adr :
		if idt not in not_adr :
			adr_pred_only.append(idt)

	with open('Diffrence_new_pred.txt', 'w') as filout_1 :
		with open('Diffrence_common.txt', 'w') as filout_2 :
			with open('Diffrence_adr_pred_only.txt', 'w') as filout_3 :
				for idt in new_pred :
					filout_1.write(idt+'\n')
				for idt in common :
					filout_2.write(idt+'\n')
				for idt in adr_pred_only :
					filout_3.write(idt+'\n')


def Proteom_arabi_Filtrage_dploc() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/blast/')

	idt = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/Filtrage_deeploc_alpha_Arabi.txt')
	dico_all = Proteom_all(path_Chlamy_arabi)

	dico = {}
	for ident in idt :
		if ident in dico_all.keys() :
			dico[ident] = ""
			dico[ident] = dico_all[ident]

	with open('Arabi_filtrage_dploc.fasta', 'w') as filout :
		for ident, seq in dico.items() :
			filout.write(ident+'\n'+seq+'\n')


def what_is_in_filtrage_deeploc_arabi() :

	os.chdir(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')
	dico = read_blast(path_pos_neg+'db_pos/Tem_pos_VS_alpha_Arabi_Filtrage_Dploc.out')
	#print(dico, len(dico))

	in_pos = list(dico.keys())
	#print(in_pos, len(in_pos))
	for i in range(len(in_pos)) :
		in_pos[i] = '>'+in_pos[i]
	#print(in_pos, len(in_pos))

	filtrage_dploc = get_idt(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/Filtrage_deeploc_alpha_Arabi.txt')
	#print(filtrage_dploc, len(filtrage_dploc))
	not_in_pos = []
	for idt in filtrage_dploc :
		if idt not in in_pos :
			not_in_pos.append(idt)
	#print(not_in_pos, len(not_in_pos))

	with open('in_pos_arabi.txt', 'w') as filout_1 :
		with open('not_in_pos_arabi.txt', 'w') as filout_2 :
			for idt in in_pos :
				filout_1.write(idt+'\n')
			for idt in not_in_pos :
				filout_2.write(idt+'\n')



def annotation_eggnog(path_files, path_new) :

	os.chdir(path_files+'Parsing/')

	files = glob.glob(path_files+'*/*.tsv')
	files.sort()
	#print(files, len(files))
	files_new = glob.glob(path_new+'*.txt')
	files_new.sort()
	#print(files_new, len(files_new))

	for f in files :
		name = f.split('/')[-2]
		print("---------------------", name, "---------------------")
		#df = pd.read_excel(f, header = None)
		#df = pd.read_excel(f)
		#df = pd.read_csv(f, sep = '\t')
		#print(df)
		idt_egg = []
		dico = {}
		with open(f, 'r') as filin :
			for line in filin :
				if line.startswith('#query') :
					col = list(line.split('\t'))
					col[-1] = col[-1].strip()
					print(col, len(col))
				if not line.startswith('#') :
					idt = '>'+line.split()[0]
					idt_egg.append(idt)
					#print(line.split()[0])
		#print(idt_egg, len(idt_egg))

		with open('idt_egg_'+name+'.txt', 'w') as filout :
			for prot in idt_egg :
				filout.write(prot+'\n')
		for n in files_new :
			if name+'.txt' in n :
				new = get_idt(n)
		#print(new, len(new))

		non_annot = []
		annot = []
		for prot in new :
			if prot not in idt_egg :
				non_annot.append(prot)
			else :
				annot.append(prot)
		#print(non_annot, len(non_annot))
		#print(annot, len(annot))

		with open('annoted_'+name+'.txt', 'w') as filout_1 :
			with open('non_annoted_'+name+'.txt', 'w') as filout_2 :
				for prot in annot :
					filout_1.write(prot+'\n')
				for prot in non_annot :
					filout_2.write(prot+'\n')

		#print([col])
		df = pd.DataFrame(index = [idt_egg], columns = [col[1:]])
		#df = pd.DataFrame(zip(idt_egg), columns = [col])
		print(df)



		

def parsing_eggnog(file) :

	pass




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
	#minus_log_evalue('Predictions/Pour_celine_comp/db_*/*_VS_*.out')
	#correspondance_acc('Predictions/dataframe_all.csv')
	#adressage_alpha('new_pred_Arabidopsis.txt', 'new_pred_Chlamy.txt') # utiliser plutot alpha_pred total
	#adressage_alpha('alpha_Arabi.txt', 'alpha_Chlamy.txt')
	#adressage_alpha_deeploc(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/', \
	#	path_Chlamy_arabi+'DEEPLOC/', path_Chlamy_arabi+'WOLFPSORT/')
	#adressage_alpha_wolfpsort(path_Chlamy_arabi+'WOLFPSORT/')
	#adressage_alpha_localizer(path_Chlamy_arabi+'LOCALIZER/')
	#intersection(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/')
	#what_is_in_filtrage_deeploc(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/', path_pos_neg+'1081_tem_pos.fasta_line')
	#is_ppr_opr(path_Chlamy_arabi+'Predictions/Pour_celine_comp/Chlamydomonas_opr_table961897.txt')
	#comp_Hedi('Predictions/comp_Hedi/arabi_chlamy_2022_02_24_predicition_ingrid.txt')
	#right_proteom_opr(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/alpha_Chlamy_VS_OPR_Chlamy.out')
	#for_cytoscape()
	#comp_pos_neg(path_Chlamy_arabi+'Predictions/Res_blast/alpha_Arabi_VS_neg.out', path_Chlamy_arabi+'Predictions/Res_blast/alpha_Arabi_VS_pos.out', \
	#path_Chlamy_arabi+'Predictions/Res_blast/alpha_Chlamy_VS_neg.out', path_Chlamy_arabi+'Predictions/Res_blast/alpha_Chlamy_VS_pos.out')
	#for_cytoscape()
	#which_opr(path_Chlamy_arabi+'Predictions/Chlamy_opr_blast/for_cytoscape_col_OPR_3.txt')
	#for_cytoscape_2()
	#comp_new_Cel()
	#for_eggNOG()
	#diff_adr()
	#Proteom_arabi_Filtrage_dploc()
	#what_is_in_filtrage_deeploc_arabi()
	#comp_new_Cel()
	#for_eggNOG()
	#for_cytoscape_2()
	annotation_eggnog(path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/eggNOG_res/', path_Chlamy_arabi+'Predictions/Pour_celine_comp/df_adressage/REEL_NEW/')







