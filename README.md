# __Annotation des protéines candidates à la régulation post-transcriptionelle des génomes des organites chez les diatomées et autres organismes photosynthétiques__

## Organismes étudiés

* Les [diatomées](https://www.researchgate.net/publication/338043776_Diatom_Molecular_Research_Comes_of_Age_Model_Species_for_Studying_Phytoplankton_Biology_and_Diversity)
* [*Chlamydomonas reinhardtii*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6713297/)
* [*Arabidopsis thaliana*](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.13687)


## But de l'étude

Le but de ce stage a été de développer un modèle de machine learning afin de pouvoir détecter les protéines en α-solénoïdes candidates à la régulation post-transcriptionelle chez les chloroplastes des diatomées en se basant sur les propriétés physico-chimiques des protéines.

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

Chaque outil prends en entrée un protéome qui a été filtré au préalable par [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) pour ne garder que les protéines sans domaines transmembranaire. Toute protéine avec un ou plusieurs domaines transmembranaires après le 68e acide aminé (avant cela correspond à un peptide d'adressage qui sera clivé et perdu à l'entrée de l'organite) sera supprimé du jeu de données.

À cela nous ajoutons le calcul des fréquences de chaque acide aminé pour chaque séquence et les valeurs d'[ACC sur Z-scales](https://pubmed.ncbi.nlm.nih.gov/32731621/) avec un lag de 4.

-----

### Scripts

Les codes répondant aux problématiques ont été rédigés principalement en python (v3.9) et en R.
Sur ce github tous les codes utilisés sont disponibles dont :    
* acc.R -> code permettant le calcul des ACC sur chacune des séquences d'un protéome. Pour utiliser le script il faut indiquer le chemin vers les protéomes étudiés qu'il lira. En sortie nous obtenons un fichier nommé

#### Calcul des ACC sur Z-scales

### Résultats

Nous avons appris notre modèle sur deux jeux de données (un positif contenant des ROGEs et un négatif ne contenant aucune ROGE). Puis nous l'avons testé et calculé la moyenne des performances visibles ci-dessous sur 1000 run :    
    
<p align="center"> 
    
| Précision   |      Sensibilité      |  Spécificité |
|----------|:-------------:|------:|
| 0.94 |  0.94 | 0.92 |     
    
</p>
    
    
    
## Remerciements

## Bibliographie


