1. Utiliser le script compile_latex.sh pour compiler le fichier latex
2. Le code permettant de calculer les éuqations équivalentes est dans le fichier /code/ et les dépendances python nécessaires sont précisées dans le requirements.txt 
3. Le latex est articulé dans le main.tex et chaque paragraphe est écrit dans un .tex du dossier /parts/

.
├── code                            # Code pour calculer les équations équivalentes
│   ├── libs                            # Fonctions propres utilisées dans le code
│   ├── notebook_equ_equiv.ipynb        # Code pour calculer les équations équivalentes
│   └── requirements.txt                # Dépendances necessaires
├── compile_latex.sh                # Script pour compiler le Latex
├── main.pdf                        # PDF de l'étude
├── main.tex                        # fichier latex principal 
├── parts                           # fichiers latex annexe 
│   ├── 0_annexes
│   ├── 0_Introduciton.tex
│   ├── 1_problemeCible
│   ├── 2_equations_equivalentes
│   ├── 3_PerteOrdre
│   └── 5_Conclusion.tex
├── README.MD
└── references.bib                  # références bibliographiques