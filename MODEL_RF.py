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

	df = pd.read_csv(path+file, sep = '\t')

	df = df.set_index(df['Unnamed: 0'], inplace = False)
	del df['Unnamed: 0']
	print(df)

	return df

def app_test_val(df, len_app_test, len_val, len_app, len_test) :

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

	df_desc = pd.concat((pd.DataFrame(train.iloc[:, 1:].columns, columns = ['variable']),
		pd.DataFrame(important, columns = ['importance'])), 
		axis = 1).sort_values(by = 'importance', ascending = False)

	print(df_desc)
	df_desc.to_csv('Importance_desc.csv', sep = '\t', header = True, index = True)

	return df_desc


def Perf_calculator(pred_test, pred_val) :

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
			#elif ind.split(' ')[-1]



	print(ldf)

	lorg = ['proteome_Chlamydomonas.fa', 'proteome_Arabidopsis_thaliana.faa', 'other']
	print(lorg[0], type(lorg[0]), len(org[0]))

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
		#total_arabi = len(df[df['org'] == lorg[1]])
		total_other = len(df[df['org'] == lorg[2]])
		#print(total_chlamy)
		print(TP_chlamy, TN_chlamy, FP_chlamy, FN_chlamy)
		print(TP_chlamy/total_chlamy, TN_chlamy/total_chlamy, FP_chlamy/total_chlamy, FN_chlamy/total_chlamy)

		Accuracy_chlamy = (TP_chlamy+TN_chlamy)/total_chlamy
		#Accuracy_arabi = (TP_arabi+TN_arabi)/total_arabi
		Accuracy_other = (TP_other+TN_other)/total_other
		print("ACCURACY CHLAMY:", Accuracy_chlamy)
		#print("ACCURACY ARABIDOPSIS:", Accuracy_arabi)
		print("ACCURACY OTHER:", Accuracy_other)
		print("---------------------")

		print(TP_arabi, TN_arabi, FP_arabi, FN_arabi)
		Sensibility_chlamy = TP_chlamy/(TP_chlamy+FN_chlamy)
		#Sensibility_arabi = TP_arabi/(TP_arabi+FN_arabi)
		Sensibility_other = TP_other/(TP_other+FN_other)
		print("SENSIBILITY CHLAMY:", Sensibility_chlamy)
		#print("SENSIBILITY ARABIDOPSIS:", Sensibility_arabi)
		print("SENSIBILITY OTHER:", Sensibility_other)
		print("---------------------")



		Specificity_chlamy = TN_chlamy/(TN_chlamy+FP_chlamy)
		#Specificity_arabi = TN_arabi/(TN_arabi+FP_arabi)
		Specificity_other = TN_other/(TN_other+FP_other)
		print("SPECIFICITY CHLAMY:", Specificity_chlamy)
		#print("SPECIFICITY ARABIDOPSIS:", Specificity_arabi)
		print("SPECIFICITY OTHER:", Specificity_other)
		#print("SENSIBILITY :", Sensibility)
		#print("SPECIFICITY :", Specificity)
		#print("ACCURACY :", Accuracy)


def idt_(path) :

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



if __name__ == '__main__' :

	path = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/'
	path_prote = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/other/'


	df = data_reading('dataframe_all.csv')
	sh_df, val, app, test = app_test_val(df, 0.90, 0.10, 0.80, 0.20)

	random_forest = model()
	#Optimal_parameters(app)

	model_res_app, score_app, importance, predictions, val_pred = Model_(random_forest, app, test, val)
	df_imp = Importance(random_forest, app, importance)
	Perf_calculator(predictions, val_pred)


	clade = which_clade(predictions, val_pred)







