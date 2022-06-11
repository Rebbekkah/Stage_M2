# __Annotation des protéines candidates à la régulation post-transcriptionelle des génomes des organites chez les diatomées et autres organismes photosynthétiques__

## Organismes étudiés

* Les [diatomées](https://www.researchgate.net/publication/338043776_Diatom_Molecular_Research_Comes_of_Age_Model_Species_for_Studying_Phytoplankton_Biology_and_Diversity)
* [*Chlamydomonas reinhardtii*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6713297/)
* [*Arabidopsis thaliana*](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.13687)


## But de l'étude

Le but de ce stage a été de développer un modèle de machine learning afin de pouvoir détecter les protéines en α-solénoïdes candidates à la régulation post-transcriptionelle chez les chloroplastes des diatomées en se basant sur les propriétés physico-chimiques des protéines.

-----

### Modèle de classification

Le modèle utilisé est un modèle de type [Random Forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html).   
Pour la prédiction nous utilisons les propriétés structurelles des protéines en α-solénoïdes adressées aux organites (appelées [ROGEs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558696/)) : ces protéines sont constituées de répétitions de séquences correspondant à des répétitions d'hélices α et reliées par des coudes, formant leur structure.   
    
Le modèle se base sur les résultats de logiciels de prédiction utilisés en standalone :
* [Radar](https://www.ebi.ac.uk/Tools/pfa/radar/) (Rapid Automatic Detection adn Alignment of Repeats) pour la détection de répétitions de séquences au sein des protéines.
* [Ard2](https://bio.tools/ard2) (Alpha-rod Repeat Detector) permet la détection des coudes entre hélices α.
* Des logiciels de prédiction intracellulaire :   
`* Deeploc`     
`* Wolfpsort`     
`* Localizer`     
`* Targetp2`     



### Scripts 

### Résultats
