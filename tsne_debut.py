import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.manifold import TSNE

'''
from scipy.spatial.distance import pdist
from sklearn.manifold.t_sne import _joint_probabilities
from scipy import linalg
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
'''
'''
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.manifold.t_sne import _joint_probabilities
from scipy import linalg
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
import seaborn as sns
'''

sns.set(rc={'figure.figsize':(11.7,8.27)})
palette = sns.color_palette("bright", 10)

sns.set(rc={'figure.figsize':(11.7,8.27)})
palette = sns.color_palette("bright", 10)




def read_seq(*fasta):

	files = list(fasta)
	#print(files)
	#print(type(files))

	dico = {}

	for file in files :
		for seq_record in SeqIO.parse(file, "fasta") :
			record = SeqIO.read(file, "fasta")
			print(record.id)
			#print(record.seq)

			seq_name = record.id
			#dico[seq_name] = ""

			my_seq = ProteinAnalysis(str(seq_record.seq))
			freq = my_seq.get_amino_acids_percent()
			dico[seq_name] = freq
	#print(dico.values())

	print("----------------------")

	return dico, df
		
def tsne(dico) :

	X = []
	#print(len(dico.keys()))

	for aa_freq in dico.values() :
		#print(aa_freq)
		#print(aa_freq.values())
		#print(list(aa_freq.values()))
		X.append(list(aa_freq.values()))
		print("-------------------------------")
		#for freq in aa_freq.values() :
			#print(freq)
	print(type(X))
	X = tuple(X)
	#df = pd.DataFrame(dico, columns = ['AA', 'FREQ'])
	#dico_2 = dico.values()
	#print(dico_2)
	#print(type(dico_2))
	df = pd.DataFrame(list(dico.items()), columns = ['ID', 'FREQ'])
	#df = pd.DataFrame([dico_2])
	print(df.iloc[0, 1])

	#tsne_em = TSNE(n_components = 2, perplexity = 30.0, verbose = 1).fit_transform(df.iloc[:, 1])

	#print(X)
	#print(type(X))
	#print(X[1])		
	
	#n_samples = len(X)
	#distances = pairwise_distances(X, metric = 'euclidean', squared = True)
	#print(distances)

	'''
	for freq in X :
		print(freq)
		
		tsne_em = TSNE(n_components = 2, perplexity = 30.0, verbose = 1).fit_transform(freq)
	'''


	'''
	l = []
	for aa_freq in dico.values() :
		#print(aa_freq)
		for freq in aa_freq.values() :
			#print(freq)
			l.append(freq)
	print(l)
	'''





if __name__ == '__main__' :

	dico_freq = read_seq("Q9SNB7.fasta", "F4HXU3.fasta")
	TSNE = tsne(dico_freq)




