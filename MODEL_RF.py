import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV


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
	print("--------------VAL DF--------------")
	print(val_data)

	app_test = df_shuffled.iloc[int(len2):, :]
	print("--------------APP_TEST DF--------------")
	print(app_test)

	len3 = len_app*len(app_test)
	len4 = len_test*len(app_test)

	df_test = app_test.iloc[:int(len3), :]
	print("--------------TEST DF--------------")
	print(df_test)

	df_app = app_test.iloc[int(len3):, :]
	print("--------------APP DF--------------")
	print(df_app)


	return df_shuffled, val_data, df_app, df_test



def Optimal_parameters(train) :
	rf = RandomForestClassifier(max_features='auto', oob_score=True, random_state=1, n_jobs=-1)

	param_grid = { "criterion" : ["gini", "entropy"], "min_samples_leaf" : [1, 5, 10],\
	 "min_samples_split" : [2, 4, 10, 12, 16], "n_estimators": [50, 100, 400, 700, 1000]}

	gs = GridSearchCV(estimator=rf, param_grid=param_grid, scoring='accuracy', cv=3, n_jobs=-1)

	gs = gs.fit(train.iloc[:, 1:], train.iloc[:, 0])

	print(gs.bestscore)
	print(gs.bestparams)
	print(gs.cvresults)



def model() :
	pass













if __name__ == '__main__' :

	path = '/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/TEMOINS_POS_NEG/outputs/'
	
	df = data_reading('dataframe_all.csv')
	sh_df, val, app, test = app_test_val(df, 0.90, 0.10, 0.80, 0.20)
	Optimal_parameters(app)










