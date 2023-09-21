###########################################################
#      #### ####   ###   #### ##### #####		  #
#      #  # #  #  #   #    #  #       #			  #
#      #### ####  #   #    #  ###     #      ##           #
#      #    #  #  #   # #  #  #       #     #  #  #  #    #
#      #    #   #  ###   ##   #####   #  #  ###    ##     #
########################################### # #### # ######
                                            #     #

Auteur : DELATTRE MATHYS 
Numéro étudiant : 28606807

######################
#  I. INTRODUCTION   #
######################

Le scripte Projet.py a été conçu pour répondre à un problème biologique : 
celui de comparer des génomes entre plusieurs espèces bactériennes afin de
venir étudier leur synthénie.

##########################
#  II. FONCTIONNALITES   #
##########################

1. Téléchargement de fichiers fasta contenant le protéome des espèces voulues
2. Réalisation d'un Blast sur le serveur obiwan
3. Implémentation des données dans la base de donnée BaseBgMuscle
4. Tester l'homologie entre chaque protéine des deux génomes
5. Réalisation d'un Dotplot entre les deux génomes (un point équivaut à une 
   homologie entre deux protéines)

#####################################
#  III. PREREQUIS ET CONFIGURATION  #
#####################################

    o Afin d'assurer le bon fonctionnement du scripte, il faut s'assurer de l'installation de 
      plusieurs librairies python.

        - tkinter
        - numpy
        - psycopg2
        - matplotlib.pyplot
        - os

    o Assurez vous d'avoir extrait le fichier.gz sur le serveur obiwan 
      Tous les fichiers requis pour le fonctionnement du scripte sont contenus dedans.

        #
       # #      ATTENTION !! Un fichier DATA_Blast va être créé à l'emplacement d'extraction
      # ! #                  Son contenu va être régulièrement vidé donc il ne faut rien 
     #######                 y mettre d'important.


    o Créez la base de donnée grâce au fichier dump

            NB : Le nom donné à la base de donnée ainsi que l'emplacement du dossier 
                 DATA_Blast devront être spécifiés lors de chaque lancement de l'application.
                 Pour gagner du temps vous pouvez directement aller renommer les variables 
                 dans le scripte :
                 -maBase : ligne 17, -LocalDirectory : ligne 19, -DATAFolderName : ligne 21

    o Vérifiez votre connexion internet si vous voulez implémenter la base de donnée

    o Vous êtes prêts à utiliser l'application !


######################################
# IV. UTILISATION DE L'APPLICATION   #
######################################

    #Nous vous conseillons de garder un terminal ouvert à coté de l'application.
    #Cela permet de voir l'avancement de plusieurs étapes

    ___________________________
    1- PREMIER ONGLET 'DOTPLOT' 

        Cet onglet sert à visualiser les données présentes sur la base de donnée.
        Le résultat de la recherche se fait sous forme d'un dotplot.
        Les points jaunes correspondent à une homologie entre les protéines des
        espèces comparées, ces protéines sont classées par rapport à leur rang (= place
        dans le génome)

        a) Choisir les espèces à comparer (Espèce 1 & Espèce 2)
           (Le jeu d'espèce disponible est l'entièreté des espèces pour lesquelles le protéome est
           présent dans la base de donnée, il n'y a pas forcément de données de BLAST pour chaque combinaison)

        b) Cocher les critères que l'on veut appliquer à notre jeu de donnée. Et écrire 
           le seuil d'homologie désiré :

           E-value : homologie si donnée < seuil

           P-ident : homologie si donnée > seuil
                     valeur entre 0 et 100 correspondant au pourcentage d'identité entre les deux séquences protéiques
           
           Q-cover : homologie si donnée > seuil
                     valeur entre 0 et 1 représentant la proportion de la séquence query recouverte par l'alignement trouvé en BLAST 
                     (1= recouvrement total, 0= pas de recouvrement)

           Q-cover : homologie si donnée > seuil
                     valeur entre 0 et 1 représentant la proportion de la séquence subject recouverte par l'alignement trouvé en BLAST 
                     (1= recouvrement total, 0= pas de recouvrement)

        c) Sélectionner les paramètres de la fenêtre de lecture :

           Window Size : int, taille de la fenêtre de lecture

           Stringence : int, seuil pour considérer les séquences protéiques comme homologues

                ex: Window Size = 10, Stringence = 4 : deux protéines seront considérées comme homologues si sur ce couple de protéine et sur les 9 suivants au moins 4 couples sont homologues.

           
        d) Appuyez sur le bouton 'Go', un dotplot devrait s'afficher.
           Si il y a des erreurs dans la saisie ou si les données ne sont pas présentes, 
           un message d'erreur apparaîtra en rouge au dessus du bouton.
        
        e) Si il n'y a pas de données pour le couple sélectionné, remplissez la base de donnée
           (onglet 'Fill DB')

        
    ____________________________
    2- DEUXIEME ONGLET 'FILL DB' 

        Cet onglet est très simple et sert à remplire la table 'all_proteomes' de la base de donnée en 
        allant télécharger les protéomes des espèces (grâce au fichier fileCSV = 'prokaryotes_complete-genomes.csv').
        Ensuite un BLAST est réalisé sur le serveur obiwan permettant de remplir la table 'blast'.
        Les fichiers sont stockés dans le dossier DATAFolderName = 'DATA_Blast/'
        Le dossier est vidé avant chaque téléchargement.

        a) Sélectionner les espèces sur lesquels on veut effectuer un BLAST

                NB: Lorsqu'il y a plusieurs données pour une même espèce, on télécharge 
                    les données les plus récentes


        b) Cliquer sur le bouton 'Lancer Recherche'

        c) La base de donnée devrait se remplir et le message 'BLAST terminé' devrait s'afficher au dessus du bouton

