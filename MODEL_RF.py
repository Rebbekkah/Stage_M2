import pandas as pd
import numpy as np
import os
import seaborn as sns
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
	df_pred.to_csv('Predictions_res_val.csv', sep = '\t', header = True, index = True)

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
			print("-----------TRAIN-----------")
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








if __name__ == '__main__' :

	path = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/'

	df = data_reading('dataframe_all.csv')
	sh_df, val, app, test = app_test_val(df, 0.90, 0.10, 0.80, 0.20)

	random_forest = model()
	Optimal_parameters(app)

	#model_res_app, score_app, importance, predictions, val_pred = Model_(random_forest, app, test, val)
	#df_imp = Importance(random_forest, app, importance)
	#Perf_calculator(predictions, val_pred)










