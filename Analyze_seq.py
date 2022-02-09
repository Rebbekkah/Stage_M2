import pandas as pd
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


		print("----dico helix :", dico)

	for cle, val in dico.items() :
		print(type(val))
		start = int(val.split(' ')[0])
		end = int(val.split(' ')[-1])
		print(start)
		print(end)

		if start > 70 :
			dico_h[cle] = {}
			dico_h[cle] = {'start' : start, 'end' : end}
		elif end < 70 :
			if cle not in id_list_ok :
				id_list_ok.append(cle)
			dico_ok = {}
			dico_ok[cle] = {'start' : start, 'end' : end}


	for idt in idt_all :
		#print(idt)
		if idt in dico_h.keys() :
			print("À SUPP : ", idt)


	#print(dico_h)
	#print(dico_ok)
	#print(id_list_ok, len(id_list_ok))
	print(idt_all, len(idt_all))



	'''
	for dic in dico.values() :
		if len(dic.keys()) == 0 :
			print("vide")
		else : 
			print("TBH")
	print(id_list)
	print(dico)
	print(len(id_list), len(dico.keys()))
	'''
	return id_list_ok

def listing(path) :

	fich = glob.glob(path+'*.tmhmm')
	prot = []
	
	for f in fich :
		prot.append(reading(path))

	return prot



if __name__ == '__main__' :

	path = "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/Celine/tests/proteome_diatom.tmhmm"
	read = reading(path)
	#proteins = listing(path)