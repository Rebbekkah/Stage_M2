# __Annotation des protéines candidates à la régulation post-transcriptionelle des génomes des organites chez les organismes photosynthétiques__

## Organismes étudiés

* Les [diatomées](https://www.researchgate.net/publication/338043776_Diatom_Molecular_Research_Comes_of_Age_Model_Species_for_Studying_Phytoplankton_Biology_and_Diversity)
* [*Chlamydomonas reinhardtii*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6713297/)
* [*Arabidopsis thaliana*](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.13687)


## But de l'étude

Le but de ce stage a été de développer un modèle de machine learning afin de __détecter les protéines en α-solénoïdes candidates à la régulation post-transcriptionelle dans le chloroplaste et la mitochondrie des organismes photosynthétiques__ en utilisant les propriétés physico-chimiques des protéines.

-----

### Modèle de classification

Le modèle utilisé est un modèle d'arbre décisionnel de type [Random Forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html).   
Pour la prédiction nous utilisons les propriétés structurelles des protéines en α-solénoïdes adressées aux organites (appelées [ROGEs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558696/)) : ces protéines sont constituées de répétitions de séquences correspondant à des répétitions d'hélices α et reliées par des coudes, formant leur structure. Aussi, elles ne possèdent pas de domaines transmembranaires.   
    
Le modèle se base sur les résultats de logiciels de prédiction utilisés en standalone :
* [Radar](https://www.ebi.ac.uk/Tools/pfa/radar/) (Rapid Automatic Detection adn Alignment of Repeats) pour la détection de répétitions de séquences au sein des protéines.
* [Ard2](https://bio.tools/ard2) (Alpha-rod Repeat Detector) permet la détection des coudes entre hélices α.
* Des logiciels de prédiction intracellulaire :   
[`Deeploc`](https://services.healthtech.dtu.dk/service.php?DeepLoc-1.0)     
[`Wolfpsort`](https://wolfpsort.hgc.jp/)     
[`Localizer`](https://localizer.csiro.au/)     
[`Targetp2`](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)     

Chaque outil prends en entrée un protéome filtré au préalable par [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) pour ne garder que les protéines sans domaines transmembranaire. Toute protéine avec un ou plusieurs domaines transmembranaires après le 68e acide aminé est supprimé du jeu de données. Avant cet acide aminé, nous considérons qu'une telle prédiction peut correspondre à un peptide d'adressage qui sera clivé et perdu à l'entrée de l'organite et nous n'en tenons pas compte.

À cela nous ajoutons le calcul des fréquences de chaque acide aminé pour chaque séquence et les valeurs d'[ACC sur Z-scales](https://pubmed.ncbi.nlm.nih.gov/32731621/) avec un lag de 4.

-----

### Scripts

Les codes répondant aux problématiques ont été rédigés principalement en python (v3.9) et en R.
Sur ce github tous les codes utilisés sont disponibles dont :       
* acc.R -> code permettant le calcul des ACC sur chacune des séquences d'un protéome. Pour utiliser le script il faut indiquer le chemin vers les protéomes étudiés qu'il lira. En sortie nous obtenons deux fichiers utilisables pour le script `Analyze_sequence.py` nommés `ACC_output_nomdufichier.txt` et `rownames_nomdufichier.txt`.    Le premier contient les résultats du calcul d'ACC sur Z-scales pour chaque séquence et le second les noms des séquences correspondantes.    
    
* Analyze_sequence.py -> Ce script prends les différents chemins jusqu'aux résultats des outils et des protéomes sur lesquels ils ont été lancés pour les lire, traiter et les analyser. En sortie ce script fournit les "nouveaux" protéomes (après analyse de TMHMM) qui seront réutilisés par la suite.        
Il produit aussi deux dataframe : `dataframe_interm.csv`(avec seulement les résultats des outils de prédiction) et `dataframe_all.csv` (outils de prédiciton, ACC sur Z-scales et fréquences d'acides aminés) qui sont les dataframes contenant les résultats des traitements de données. Seule `dataframe_all.csv` sera donnée au modèle. 
    
* MODEL_RF.py -> Script contenant le code du modèle Random Forest. Il faut lui donner le chemin vers la matrice `dataframe_all.csv`. En sortie il donne :
-> Pour les témoins postifis et négatifs : les résultats des prédictions (0 = positif et 1 = négatif) ainsi que des fichiers contenant les protéines correctement prédites et celles non correctement prédites. Il fournit aussi des fichiers contenant les protéines prédites comme ROGEs et les prédites prédites comme n'étant pas des ROGEs. Avec cela s'additionne un histogramme de l'importance des descripteurs. Plus un descripteur est important plus il permet au modèle de prédire correctement la nature des protéines.
-> Pour les protéomes à analyser : Fournit deux fichiers contenant les protéines prédites comme ROGEs ou non.
* Analyze_results.py -> Script qui a servit à analyser les résultats du modèle. Pour les diatomées il faut supprimer les colonnes de la dataframe contenant les résultats des logiciels de prédiction d'adressage à l'aide de la fonction "`for_Diatoms()`".
    
* TSNE_plot_freq_ACC.py -> script qui permet d'avoir une visualisation graphique sous forme de TSNE (T-distribution Neighbor joining Embbedding) qui permet la réduction de dimensionnalité. Le graphe représente la TSNE des fréquences d'acides aminés et les ACC sur Z-scales pour les protéomes d'entrée. Il faut lui fournir le chemin vers les protéomes et des résultats d'ACC.


-----

### Résultats

Nous avons appris notre modèle sur deux jeux de données (un positif contenant des ROGEs et un négatif ne contenant aucune ROGE). Puis nous l'avons testé et calculé la moyenne des performances sur 1000 run (ci-dessous) :    
    
<p align="center"> 
    
| Précision   |      Sensibilité      |  Spécificité |
|----------|:-------------:|------:|
| 0.94 |  0.94 | 0.92 |     
    
</p>
    
À la suite de l'apprentissage il est possible d'utiliser le modèle sur les protéomes cibles et de déterminer les potentielles candidates ROGEs.
    
-----

## Autres fichiers importants

Dans ce github vous trouverez également les protéomes utilisés ainsi que des exemples de sorties pour chacun des logiciels.   
Le rapport du stage est aussi disponible au format PDF.

-----

## Remerciements

Je remercie Angela Falciatore ainsi que toute l'équipe Photosynthèse de l'IBPC pour m'avoir acceuillie dans leur laboratoire.  
Je remercie particulièrement __*Ingrid Lafontaine*__ et __*Céline Cattelin*__ pour leur dévouement et leur acceuil, ainsi que pour le suivi qu'elles m'ont apporté tout au long de mon stage.   
Je remercie également __*Clotilde Garrido*__ pour les directives qu'elle a pu m'apporter.

-----

## Bibliographie

La bibliographie est disponible [ici](https://github.com/Rebbekkah/Stage_M2/tree/main/Bibliographie).
