""" This allow perform a Random Forest model in order to analyze and classify
	proteins as alpha-solenoïds ROGEs (Regulator Of Genome Expression) in
	plants organelles

------------------------------------------------------------------
Rebecca GOULANCOURT
M2 BIOLOGIE - INFORMATIQUE
Université de Paris 2021 - 2022
Stage M2 - supervisor : Ingrid Lafontaine & Céline Cattelin
------------------------------------------------------------------

"""

import pandas as pd
import numpy as np
import os
import glob
import seaborn as sns
from os.path import basename
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix


def data_reading(file) :
	''' Read a csv file that contains a dataframe of the results for the
	"Analyze_seq.py" code 

	
	Parameter
	---------
	file : str
		file to be read

	Returns
	-------
	df : Dataframe
		dataframe that has been read

	'''

	df = pd.read_csv(file, sep = '\t')

	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']
	print(df)

	return df

def app_test_val(df, len_app_test, len_val, len_app, len_test) :
	''' Use a dataframe to split and create learning, test and validation dataframes
	--> Dire l'organisation des données


	Parameter
	---------
	df : Dataframe
		dataframe of results

	Returns
	-------
	df_shuffled : Dataframe
		the same dataframe in input but shuffled

	val_data : Dataframe
		dataframe that contains validation dataset

	df_app : Dataframe
		dataframe that contains train dataset

	df_test : Dataframe
		dataframe that contains test dataset

	'''

	df_shuffled = df.sample(frac = 1)
	print(df_shuffled)

	len1 = len_app_test*len(df_shuffled)
	len2 = len_val*len(df_shuffled)

	val_data = df_shuffled.iloc[:int(len2), :]
	df_shuffled = df_shuffled.iloc[int(len2):, :]
	print("--------------VAL DF--------------")
	print(val_data)
	print("--------------NEW DF SHUFFLED--------------")
	print(df_shuffled)


	#app_test = df_shuffled.iloc[int(len2):, :]
	app_test = df_shuffled
	print("--------------APP_TEST DF--------------")
	print(app_test)

	len3 = len_app*len(app_test)
	len4 = len_test*len(app_test)

	df_app = app_test.iloc[:int(len3), :]
	print("--------------APP DF--------------")
	print(df_app)

	df_test = app_test.iloc[int(len3):, :]
	print("--------------TEST DF--------------")
	print(df_test)


	return df_shuffled, val_data, df_app, df_test



def Optimal_parameters(train) :

	rf = RandomForestClassifier(max_features = 'auto', oob_score = True, random_state = 1, n_jobs = -1)

	param_grid = { "criterion" : ["gini", "entropy"], "min_samples_leaf" : [1, 5, 10],\
	 "min_samples_split" : [2, 4, 10, 12, 16], "n_estimators": [50, 100, 400, 700, 1000]}

	gs = GridSearchCV(estimator = rf, param_grid = param_grid, scoring = 'accuracy', cv = 3, n_jobs = -1)

	gs = gs.fit(train.iloc[:, 1:], train.iloc[:, 0])

	print(gs.best_score_)
	print(gs.best_params_)
	print(gs.cv_results_)



def model() :
	''' Set the model and its parameters

	Parameters
	----------
	None

	Returns
	-------
	rf : sklearn 'RandomForestClassifier'
		the model we will use to predict

	'''
	
	rf = RandomForestClassifier(criterion = 'gini', 
							 n_estimators = 500,
							 min_samples_split = 30,
							 min_samples_leaf = 1,
							 max_features = 'auto',
							 oob_score = True,
							 random_state = 1,
							 n_jobs = -1,
							 verbose = 1)

	'''
	rf = RandomForestClassifier(criterion = 'gini', 
							 n_estimators = 700,
							 min_samples_split = 10,
							 min_samples_leaf = 1,
							 max_features = 'auto',
							 oob_score = True,
							 random_state = 1,
							 n_jobs = -1)
	'''
	return rf


def Model_(rf, train, test, val) :
	''' Train the model on the train dataset, test it on the test dataset 
	and validation dataset

	Parameters
	----------
	rf : sklearn 'RandomForestClassifier'
		the model used to predict

	train : Dataframe
		train dataset

	test : Dataframe
		test dataset

	val : Dataframe
		validation dataset

	Returns
	-------
	res : sklearn 'RandomForestClassifier'
		results of the prediction on the train dataset

	score : float
		accuracy calculated by the rf itself

	imp : str
		importance of each descriptors

	df_pred : Dataframe
		dataframe of prediciton results on test dataset

	df_pred_val : Dataframe
		dataframe of prediciton results on validation dataset

	Writes
	------
	Predictions_res_test.csv : .csv file
		df_pred

	Predictions_res_val.csv : .csv file
		df_pred_val

	'''

	res = rf.fit(train.iloc[:, 1:], train.iloc[:, 0])
	score = rf.oob_score_
	imp = rf.feature_importances_
	print("%.4f" % score)
	print("%.4f" % (1-score))

	print(test.iloc[:, 0], type(test.iloc[:, 0]))
	print(train)
	pred = rf.predict(test.iloc[:, 1:])
	df_pred = pd.DataFrame(pred, columns = ['pred'])
	df_pred.index = test.index
	df_pred = pd.concat([test.iloc[:, 0], df_pred], axis = 1)
	print("--------PRED\n", df_pred)



	pred_val = rf.predict(val.iloc[:, 1:])
	df_pred_val = pd.DataFrame(pred_val, columns = ['pred'])
	df_pred_val.index = val.index
	df_pred_val = pd.concat([val.iloc[:, 0], df_pred_val], axis = 1)

	df_pred.to_csv('Predictions_res_test.csv', sep = '\t', header = True, index = True)
	df_pred_val.to_csv('Predictions_res_val.csv', sep = '\t', header = True, index = True)

	return res, score, imp, df_pred, df_pred_val



def Importance(rf, train, important) :
	''' Find and writes as all the descriptors that the model model considers
		the most important for the prediction (helps him)
		Rank the descriptors by ascending

	Parameters
	----------
	rf : sklearn 'RandomForestClassifier'
		the model we will use to predict

	train : Dataframe
		train dataset
	
	important : str
		importance of each descriptors


	Writes
	------
	Importance_desc.csv : .csv file
		csv file that contains the dataframe of col1 : number of column, col2 : descriptor name,
		col3 : importance


	Returns
	-------
	df_desc : Dataframe
		dataframe of the ascended rank descriptors
	
	'''


	df_desc = pd.concat((pd.DataFrame(train.iloc[:, 1:].columns, columns = ['variable']),
		pd.DataFrame(important, columns = ['importance'])), 
		axis = 1).sort_values(by = 'importance', ascending = False)

	print(df_desc)
	df_desc.to_csv('Importance_desc.csv', sep = '\t', header = True, index = True)

	#sns.barplot(x = df_desc[0], y = df_desc[1])
	sns.barplot(x = 'variable', y = 'importance', data = df_desc)
	plt.xticks(rotation = 'vertical')
	plt.tick_params(axis = 'x', labelsize = 8)
	plt.xlabel('Feature importance Score')
	plt.ylabel('Features')
	plt.title('Visualizing important features')
	#plt.show()


	return df_desc


def Perf_calculator(pred_test, pred_val) :
	''' Function that computes the performance of the model

	Parameters
	----------
	pred_test : Dataframe
		dataframe of the results of the predictions on the test dataset

	pred_val : Dataframe
		dataframe of the results of the predictions on the validation dataset

	Plots
	-----
	Heatmap of performance

	Writes
	------
	Good_pred.txt : .txt file
		protein identifier that were well predicted

	Bad_pred.txt : .txt file
		protein identifier that were wrong predicted		

	Returns
	-------
	index_wp : list
		list of protein identifier well predicted

	indew_bp : list
		list of protein identifier wrong predicted

	'''

	predictions = [pred_test, pred_val]
	index_wp = []
	index_bp = []

	for p in predictions :
		TP = 0
		FP = 0
		TN = 0
		FN = 0
		BP = 0
		GP = 0
		if p is pred_test :
			print("-----------TEST-----------")
		elif p is pred_val :
			print("-----------VAL-----------")
		print(len(p))
		for i in range(len(p)) :
			if p.iloc[i, 0] == 1 and p.iloc[i, 1] == 1 :
				GP += 1
				TP += 1
				index_wp.append(p.index[i])
			elif p.iloc[i, 0] == 0 and p.iloc[i, 1] == 0 :
				GP += 1
				TN += 1
				index_wp.append(p.index[i])
			elif p.iloc[i, 0] == 0 and p.iloc[i, 1] == 1 :
				BP += 1
				FP += 1
				index_bp.append(p.index[i])
			elif p.iloc[i, 0] == 1 and p.iloc[i, 1] == 0 : 
				BP += 1
				FN += 1
				index_bp.append(p.index[i])


		total = len(p)
		print(TP, TN, FP, FN)
		print(TP/total, TN/total, FP/total, FN/total)
		print(GP, BP)

		Accuracy = (TP+TN)/total
		Sensibility = TP/(TP+FN)
		Specificity = TN/(TN+FP)
		print("SENSIBILITY :", Sensibility)
		print("SPECIFICITY :", Specificity)
		print("ACCURACY :", Accuracy)


		y_true = p['type']
		y_pred = p['pred']
		data = confusion_matrix(y_true, y_pred)
		df_cm = pd.DataFrame(data, columns = np.unique(y_true), index = np.unique(y_true))
		df_cm.index.name = 'Actual'
		df_cm.columns.name = 'Predicted'
		plt.figure()
		sns.heatmap(df_cm, cmap = "Blues", annot = True)
		if p is pred_test :
			plt.title("Heatmap of Performances on the test dataset")
		elif p is pred_val :
			plt.title("Heatmap of Performances on the validation dataset")
		#plt.show()

		print("Good pred : ", index_wp, "\n", "Bad pred : ", index_bp)

	with open('Good_pred.txt', 'w') as filout :
		for idt in index_wp :
			filout.write(idt+"\n")
	with open('Bad_pred.txt', 'w') as filout :
		for idt in index_bp :
			filout.write(idt+"\n")

	return index_wp, index_bp


def which_clade(df_test, df_val) :
	''' Find the organisme for each identifier
		then computes the perfomances of the model for Chlamydomonas and Arabidopsis

	Parameters
	----------
	df_test : Dataframe
		test dataset


	df_val : Dataframe
		validation dataset

	Returns
	-------
	None

	'''


	dico = idt_(path_prote)
	#print(dico.keys(), len(dico.keys()))
	#print(df_test)
	#print(df_val)

	ldf = [df_test, df_val]
	for df in ldf :
		df['org'] = 0

	for df in ldf :
		for ind in list(df.index) :
			#print(ind)
			for org, ident in dico.items() :
				#print(org)
				if ind in ident :
				#	print("oui")
					df.loc[ind, 'org'] = org
				#else :
				#	print("non")
			if df.loc[ind, 'org'] == 0 :
				df.loc[ind, 'org'] = 'other'


		for ind in list(df.index) :
			#print(ind.split(' ')[-1])
			#print(ind, type(ind))
			if ind.startswith('>Cre') :
				df.loc[ind, 'org'] = 'proteome_Chlamydomonas.fa'
			elif ind.startswith('>AT') :
				df.loc[ind, 'org'] = 'proteome_Arabidopsis_thaliana.faa'



	print(ldf)

	lorg = ['proteome_Chlamydomonas.fa', 'proteome_Arabidopsis_thaliana.faa', 'other']

	for df in ldf :
		TP_other = 0
		TP_chlamy = 0
		TP_arabi = 0
		FP_other = 0
		FP_chlamy = 0
		FP_arabi = 0
		TN_other = 0
		TN_chlamy = 0
		TN_arabi = 0
		FN_other = 0
		FN_chlamy = 0
		FN_arabi = 0
		if df is df_test :
			print("-----------TEST-----------")
		elif df is df_val :
			print("-----------VAL-----------")
		for i in range(len(df)) :
			if df.iloc[i, 0] == 1 and df.iloc[i, 1] == 1 :
				if df.iloc[i, 2] == lorg[0] :
					TP_chlamy += 1
				elif df.iloc[i, 2] == lorg[1] :
					TP_arabi += 1
				else :
					TP_other += 1
			elif df.iloc[i, 0] == 0 and df.iloc[i, 1] == 0 :
				if df.iloc[i, 2] == lorg[0] :
					TN_chlamy += 1
				elif df.iloc[i, 2] == lorg[1] :
					TN_arabi += 1
				else :
					TN_other += 1
			elif df.iloc[i, 0] == 0 and df.iloc[i, 1] == 1 :
				if df.iloc[i, 2] == lorg[0] :
					FP_chlamy += 1
				elif df.iloc[i, 2] == lorg[1] :
					FP_arabi += 1
				else :
					FP_other += 1
			elif df.iloc[i, 0] == 1 and df.iloc[i, 1] == 0 : 
				if df.iloc[i, 2] == lorg[0] :
					FN_chlamy += 1
				elif df.iloc[i, 2] == lorg[1] :
					FN_arabi += 1
				else :
					FN_other += 1


		total_chlamy = len(df[df['org'] == lorg[0]])
		total_arabi = len(df[df['org'] == lorg[1]])
		total_other = len(df[df['org'] == lorg[2]])

		Accuracy_chlamy = (TP_chlamy+TN_chlamy)/total_chlamy
		Accuracy_arabi = (TP_arabi+TN_arabi)/total_arabi
		Accuracy_other = (TP_other+TN_other)/total_other
		print("ACCURACY CHLAMY:", Accuracy_chlamy)
		print("ACCURACY ARABIDOPSIS:", Accuracy_arabi)
		print("ACCURACY OTHER:", Accuracy_other)
		print("---------------------")

		if (TP_chlamy+FN_chlamy) != 0 :
			Sensibility_chlamy = TP_chlamy/(TP_chlamy+FN_chlamy)
			print("SENSIBILITY CHLAMY:", Sensibility_chlamy)
		else :
			print("SENSIBILITY CHLAMY: non calculable --> (TP_chlamy+FN_chlamy) = 0")
			print("TP = ", TP_chlamy, "FN = ", FN_chlamy)
		if (TP_arabi+FN_arabi) != 0 :
			Sensibility_arabi = TP_arabi/(TP_arabi+FN_arabi)
			print("SENSIBILITY ARABIDOPSIS:", Sensibility_arabi)
		else :
			print("SENSIBILITY ARABIDOPSIS: non calculable --> (TP_arabi+FN_arabi) = 0")
			print("TP = ", TP_arabi, "FN = ", FN_arabi)
		if (TP_other+FN_other) != 0 :
			Sensibility_other = TP_other/(TP_other+FN_other)
			print("SENSIBILITY OTHER:", Sensibility_other)
		else :
			print("SENSIBILITY OTHER: non calculable --> (TP_other+FN_other) = 0")
			print("TP = ", TP_other, "FN = ", FN_other)
		
		print("---------------------")


		if (TN_chlamy+FP_chlamy) != 0 :
			Specificity_chlamy = TN_chlamy/(TN_chlamy+FP_chlamy)
			print("SPECIFICITY CHLAMY:", Specificity_chlamy)
		else : 
			print("SPECIFICITY CHLAMY: non calculable --> (TN_chlamy+FP_chlamy) = 0")
			print("TN =", TN_chlamy, "FP = ", FP_chlamy)
		if (TN_arabi+FP_arabi) != 0 :
			Specificity_arabi = TN_arabi/(TN_arabi+FP_arabi)
			print("SPECIFICITY ARABIDOPSIS:", Specificity_arabi)
		else : 
			print("SPECIFICITY ARABIDOPSIS: non calculable --> (TN_arabi+FP_arabi) = 0")
			print("TN =", TN_arabi, "FP = ", FP_arabi)
		if (TN_other+FP_other) != 0 :
			Specificity_other = TN_other/(TN_other+FP_other)
			print("SPECIFICITY OTHER:", Specificity_other)
		else : 
			print("SPECIFICITY OTHER: non calculable --> (TN_other+FP_other) = 0")
			print("TN =", TN_other, "FP = ", FP_other)

def idt_(path) :
	''' Function that lists for each organism in a repertory all its protein identifiers

	Parameters
	----------
	path : str
		path to the proteoms

	Returns
	-------
	dico : dict
		dictionnary (key = organism and values = identifiers (as a list))

	'''

	file = glob.glob(path+'*.f*'+'a')

	dico = {}
	for f in file :
		print(basename(f))
		dico[basename(f)] = ""
		l = []
		with open(f, 'r') as filin :
			for line in filin :
				if line.startswith('>') :
					l.append(line.split(' ')[0])
			dico[basename(f)] = l

	return dico



def To_Predict(path, rf, file, name) :
	''' Use the model to predict its alpha_solenoid carateristics on a dataset

	Parameters
	----------
	path : str
		path to the dataframe

	rf : sklearn 'RandomForestClassifier'
		the model we will use to predict

	file : str
		name of the file that contains the dataframe

	name : str
		name of the output file

	Returns
	-------
	None

	Writes
	------
	df_alpha : dataframe
		dataframe of the proteins predicted as alpha-solenoïds

	df_other : dataframe
		dataframe of the proteins predicted as non alpha-solenoïds

	'''

	print('----------NEW PREDICTIONS----------')
	os.chdir(path_Chlamy_arabi+'Predictions/')

	df = data_reading(path+file)

	pred = rf.predict(df.iloc[:, 1:])
	print('PRED\n', pred, len(pred))

	df_pred = pd.DataFrame(pred, columns = ['pred'])

	df_pred.index = df.index

	print('DF_PRED\n', df_pred)
	
	df_pred.to_csv('Predictions_'+name+'.csv', sep = '\t', header = True, index = True)

	#df_alpha = pd.DataFrame()
	#df_other = pd.DataFrame()
	alpha = []
	other = []
	#i = 0
	#j = 0
	for index, elem in enumerate(df_pred['pred']) :
		if elem == 0 :
			alpha.append(df_pred.index[index])
			#df_alpha.iloc[i] = df_pred.iloc[index]
			#i =+ 1
		elif elem == 1 :
			other.append(df_pred.index[index])
			#df_other.iloc[j] = df_pred.iloc[index]
			#j += 1

	print('nb alpha :', len(alpha))
	print('nb non alpha :', len(other))
	#df_alpha.to_csv('Predictions_alpha_'+name+'.csv', sep = '\t', header = True, index = True)
	#df_other.to_csv('Predictions_other_'+name+'.csv', sep = '\t', header = True, index = True)

	with open('prot_alpha.txt', 'w') as filout :
		for a in alpha :
			filout.write(a+'\n')

	with open('prot_non_alpha.txt', 'w') as filout :
		for o in other :
			filout.write(o+'\n')


	return alpha, other, df_pred


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



if __name__ == '__main__' :

	#path = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/neg_pos/'
	#path_prote = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/other/'
	path_Chlamy_arabi = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/RF/Chlamy_Arabi/results/"
	path_pos_neg = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/"
	path_method_Cel = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/methode_1_2_Celine/'

	os.chdir(path_Chlamy_arabi)

	'''
	df = data_reading(path+'dataframe_all.csv')
	sh_df, val, app, test = app_test_val(df, 0.90, 0.10, 0.80, 0.20)

	random_forest = model()
	#Optimal_parameters(app)

	model_res_, score_app, importance, predictions, val_pred = Model_(random_forest, app, test, val)
	df_imp = Importance(random_forest, app, importance)
	Perf_calculator(predictions, val_pred)


	#clade = which_clade(predictions, val_pred)

	alphasol, nonalphasol, df_pred = To_Predict(path_Chlamy_arabi, random_forest, 'dataframe_all.csv', 'Chlamy_Arabi')
	'''
	#which_proteom()
	#is_pos_neg()
	#dico_imp = select_imp('Importance_desc.csv')
	#comp_res_Celine(path_method_Cel+'*/*')
	#comp_methode_2(path_method_Cel+'M2_*/*')
	#sep_alpha()
	#proteom_alpha()
