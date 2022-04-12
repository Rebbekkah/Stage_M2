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
	print("--------------")
	for d in ldf :
		#truc.append(d[(d['wolfpsort'] > 0.5) & (d['deeploc'] > 0.5) & (d['trp2'] == 1.0) |\
		#	(d['deeploc'] == 2.0) & (d['localizer'] == 1.0)])
		truc.append(d[(d.wolfpsort > 0.5) & (d.deeploc > 0.5) & ((d.trp2 == 1.0) | \
							(d.trp2 == 2.0)) & (d.localizer == 1.0)])
	print(truc)

	print("-----LALAALLALALALALA-----")


	for t in truc :
		pass
		#for index, elem in enumerate(t) :
			#t.iloc[index, :] = df.loc[t.index[index], adress]
		#print(list(t.index), len(list(t.index)))
		#for index in list(t.index) :
			#print(df.loc[index, :])
			#print(t[t[index, :]])
			#t.loc[index, :] = df.loc[index, adress]
		#	t.loc[index, adress].replace(['NaN'], df.loc[index, adress], inplace = True)
			#t[t[index]]
	print(truc)



'''
		#for software in adress :
			#for elem in df[software] :
		#for k in range(len(d)) :
		#	print(d['wolfpsort'].iloc[k] > 0.5)
		#	if (d['wolfpsort'].iloc[k] > 0.5) == True :
		#		if (d['deeploc'].iloc[k] > 0.5) == True :
		#			if d['trp2'].iloc[k] == 1.0 :
		#				if d['trp2'].iloc[k] == 1.0 :
		#					if d['localizer'.iloc[k]] == 1.0 :
		#						print('ok')
					



			#print(d.iloc[k, 'wolfpsort'])
			#if d['wolfpsort'].iloc[k] > 0.5 | d['deeploc'].iloc[k] > 0.5 and \
			#	d.iloc[k, 'trp2'] == 1.0 | d.iloc[k, 'trp2'] == 2.0 and d.iloc[k, 'localizer'] == 1.0 :
			#	if d.index[k] not in idt_adr :
			#		idt_adr.append(d.index[k])
			#	if d is df_adr_Chl :
			#		if d.index[k] not in idt_adr_Chl :
			#			idt_adr_Chl.append(d.index[k])
			#	elif d is df_adr_Arabi :
			#		if d.index[k] not in idt_adr_Arabi :
			#			idt_adr_Arabi.append(d.index[k])

			
			if d.iloc[k, 'wolfpsort'] > 0.5 or d.iloc[k, 'deeploc'] > 0.5 and \
				d.iloc[k, 'trp2'] == 1.0 or d.iloc[k, 'trp2'] == 2.0 and d.iloc[k, 'localizer'] == 1.0 :
				if d.index[k] not in idt_adr :
					idt_adr.append(d.index[k])
				if d is df_adr_Chl :
					if d.index[k] not in idt_adr_Chl :
						idt_adr_Chl.append(d.index[k])
				elif d is df_adr_Arabi :
					if d.index[k] not in idt_adr_Arabi :
						idt_adr_Arabi.append(d.index[k])
			


				
				if software == 'wolfpsort' or software == 'deeploc' :
					if elem > 0.5 :
						if d.index[k] not in idt_adr :
							idt_adr.append(d.index[k])
						if d is df_adr_Chl :
							if d.index[k] not in idt_adr_Chl :
								idt_adr_Chl.append(d.index[k])
						elif d is df_adr_Arabi :
							if d.index[k] not in idt_adr_Arabi :
								idt_adr_Arabi.append(d.index[k])
				if software == 'trp2' :
					if elem == 1.0 or elem == 2.0 :
						if d.index[k] not in idt_adr :
							idt_adr.append(d.index[k])
						if d is df_adr_Chl :
							if d.index[k] not in idt_adr_Chl :
								idt_adr_Chl.append(d.index[k])
						elif d is df_adr_Arabi :
							if d.index[k] not in idt_adr_Arabi :
								idt_adr_Arabi.append(d.index[k])
				if software == 'localizer' :
					if elem == 1.0 :
						if d.index[k] not in idt_adr :
							idt_adr.append(d.index[k])
						if d is df_adr_Chl :
							if d.index[k] not in idt_adr_Chl :
								idt_adr_Chl.append(d.index[k])
						elif d is df_adr_Arabi :
							if d.index[k] not in idt_adr_Arabi :
								idt_adr_Arabi.append(d.index[k])
				
			#k += 1

	
	lidt = [idt_adr, idt_adr_Chl, idt_adr_Arabi]

	print("----------------------")
	print(idt_adr, len(idt_adr))
	print(idt_adr_Chl, len(idt_adr_Chl))
	print(idt_adr_Arabi, len(idt_adr_Arabi))

	l = ['all', 'Chl', 'Arabi']

	i = 0
	for newf in lidt :
		with open('prot_adress_'+l[i]+'.txt', 'w') as filout :
			for idt in newf :
				filout.write(idt+"\n")
		i += 1


	print(df.loc['>NP_001030886.1', adress])
	print(df.loc['>NP_001030703.1', adress])
	print(df.loc['>NP_001030656.1', adress])
	print(df.loc['>NP_001030887.1', adress])
	


	
	i = 0
	for l in lidt :
		new = []
		for elem in l :
			if elem not in new :
				new.append(elem)
		#print(len(new))
		#print("---")
		lidt[i] = new
	
	
	i = 0
	for idt in lidt :
		break
		new = []
		for ident in idt :
			if ident not in new :
				new.append(ident)
		print(len(new))
		print("lalalal")
		lidt[i] = new
		i += 1
	
	#print(lidt)
	#print(len(lidt[0]))
	#print("----------\n", idt_adr_Arabi)

	
	for idt in lidt :
		print(len(idt))
		idt = list(set(idt))
		print(len(idt))
		break
		lidt[i] = idt
		i += 1

	
'''


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
	#adressage_alpha('new_pred_Arabidopsis.txt', 'new_pred_Chlamy.txt')
	is_ppr_opr(path_Chlamy_arabi+'Predictions/Pour_celine_comp/Chlamydomonas_opr_table961897.txt')

