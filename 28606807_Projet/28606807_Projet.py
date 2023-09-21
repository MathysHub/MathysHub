
# LE CODE EST DONNE A BUT INFORMATIF UNIQUEMENT IL FAUT LE FAIRE TOURNER SUR UN SERVEUR SPECIAL POUR QU'IL FONCTIONNE
# ssh ...

import tkinter as tk 
from tkinter import ttk
import numpy as np
import psycopg2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import os

########################
# VARIABLES A RENOMMER #
########################

maBase ='BaseBgMuscle'              #Le nom de votre base de donnée

LocalDirectory = '/home/mathys/'    #L'adresse où le zip est extrait (ne pas oublier le '/' à la fin)

DATAFolderName = 'DATA_Blast/'      #Le nom du dossier où les données seront temporairement rangées 
                                    #(si vous avez renommé le dossier)
                                    #(ne pas oublier le '/' à la fin)


fileCSV = 'prokaryotes_complete-genomes.csv'    #Pas besoin de renommer


###############################################################################
# CODE
###############################################################################


def readfaa(nomfi) :
    ''' Cette fonction a pour objectif de lire un fichier qui contient le protéome d'une espèce
    en format fasta. 
    INPUT1  : nom du fichier
    OUTPUT1 : array de id Prot encodée
    OUTPUT2 : array de lg de seq
    OUTPUT3 : array de rang de la prot (=position)'''

    fichier = open(nomfi,'r')

    IdProts = []
    LenSeq = []
    rangProts = []

    
    lg = 0
    for line in fichier :
        line = line.rstrip()
        if line[0] == '>' :

            IdProts.append(line.split()[0][5:])
            if lg != 0 : 
                LenSeq.append(lg)
            lg = 0

        else : 
            lg += len(line)


    LenSeq.append(lg)
    IdProts = np.array(IdProts)
    LenSeq = np.array(LenSeq)
    RgProts = np.arange(1,len(IdProts)+1)

    fichier.close()
    
    return IdProts ,LenSeq, RgProts

def fill_prots_db(nomfi, name_specie, DB) :

    conn = psycopg2.connect(database = DB)
    cur = conn.cursor()
    sql_cmd = "INSERT INTO all_proteomes (protein_id, len_seq, rang_seq, name_specie) VALUES (%s, %s, %s, %s)"

    IdProts ,LenSeq, RgProts = readfaa(nomfi)

    for i in range(len(IdProts)) :
        cur.execute(sql_cmd , (IdProts[i], int(LenSeq[i]), int(RgProts[i]), name_specie))

    conn.commit()
    cur.close()
    conn.close()

def fill_blast_db(nomfi, DB) :
    ''' Cette fonction a pour objectif de lire un fichier qui contient le blast entre deux espèces. 
    et de remplir notre base de donnée (table blast)
    INPUT1  : nom du fichier
    OUTPUT1 : array de id de seq query
    OUTPUT2 : array de id de seq subject
    OUTPUT3 : pourcentage d'identite entre les deux seqs
    OUTPUT4 : e-value
    OUTPUT5 : query-cover proportion de query recouvert par le hit
    OUTPUT6 : subject-cover proportion de subject recouvert par le hit
    Mets ces 6 outputs dans la base de donnée.
    '''

    conn = psycopg2.connect(database = DB)
    cur = conn.cursor()
    sql_cmd = "INSERT INTO blast (query_id, subject_id, pidentite, e_value, query_cover, subject_cover, query_specie, subject_specie, query_rang, subject_rang) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    sql_search = "SELECT len_seq, rang_seq, name_specie FROM all_proteomes WHERE protein_id = %s"

    fichier = open(nomfi,'r')

    ct = 0
    for line in fichier :
        if line[0] != "#" :
            split = line.split()
            

            qstart = split[6]
            qend = split[7]
            sstart = split[8]
            send = split[9]

            cur.execute(sql_search, (split[0][4:],) )
            qlen,q_rang,q_specie =  cur.fetchall()[0]
            

            cur.execute(sql_search, (split[1][4:],) )
            slen,s_rang,s_specie =  cur.fetchall()[0]

            
            cur.execute(sql_cmd , (split[0][4:], split[1][4:], float(split[2]), float(split[10]) , (int(qend)-int(qstart)) / int(qlen) , (int(send)-int(sstart)) / int(slen) , q_specie, s_specie,int(q_rang) , int(s_rang)))
            
            ct+=1
            
            
        
        if ct%2000 == 0 :
            print('remplissage DB : ',ct)
    
    conn.commit()
    cur.close()
    conn.close()
    fichier.close()
    
class ConfigureWindow :
    def __init__(self,DB,dir,fold,CSV) :

        self.maBase = DB
        self.LocalDirectory = dir
        self.DATAFolderName = fold
        self.fileCSV = CSV
        
        self.root = tk.Tk()

        #Lab et Laberror

        self.LabHead = tk.Label(self.root, text = 'Veuillez vérifier les éléments suivants :')
        self.LabHead.grid(row = 0, column = 2)
        self.LabErrorInit = tk.Label(self.root, text = '')
        self.LabErrorInit.grid(row = 1, column = 2)

    

        self.txt_mabase = tk.StringVar()
        self.txt_mabase.set(self.maBase)
        self.bool_mabase = tk.BooleanVar()
        self.chk_mabase = tk.Checkbutton(self.root,  variable = self.bool_mabase, command = lambda : self.toggle_entry_state('mabase') ,text="   DataBase :")
        self.entry_mabase = tk.Entry(self.root, width = 20, textvariable = self.txt_mabase, state = 'disabled')
        
        self.bool_directory = tk.BooleanVar()
        self.chk_directory = tk.Checkbutton(self.root,  variable = self.bool_directory, command = lambda : self.toggle_entry_state('directory') , text="  /Directory/ :")
        self.txt_directory = tk.StringVar()
        self.txt_directory.set(self.LocalDirectory)
        self.entry_directory = tk.Entry(self.root, width = 20 , textvariable = self.txt_directory, state = 'disabled')
        
        self.bool_folder = tk.BooleanVar()
        self.chk_folder = tk.Checkbutton(self.root,  variable = self.bool_folder, command = lambda : self.toggle_entry_state('folder'), text="     Dossier/ :")
        self.txt_folder = tk.StringVar()
        self.txt_folder.set(self.DATAFolderName)
        self.entry_folder = tk.Entry(self.root, width = 20, textvariable = self.txt_folder,state = 'disabled')

        self.entry_mabase.grid(row = 2, column = 3)
        self.chk_mabase.grid(row = 2, column = 2)
        self.entry_directory.grid(row = 3, column = 3)
        self.chk_directory.grid(row = 3, column = 2)
        self.entry_folder.grid(row = 4, column = 3)
        self.chk_folder.grid(row = 4, column = 2)
        
        # Bouton pour lancer l'interface

        self.button_close = tk.Button(self.root, text="Valider", command=self.execute_appli)
        self.button_close.grid()


        self.root.mainloop()

    def toggle_entry_state(self, nom_entry) :
        ''' attention l'attribut doit être passé en string'''
        if getattr(self,'bool_' + nom_entry).get() :
            getattr(self, 'entry_' + nom_entry).configure(state = 'normal')
        else :
            getattr(self, 'entry_' + nom_entry).configure(state = 'disabled', text = '')

    def execute_appli(self):
        self.LabErrorInit.config(text = '')

        self.maBase = self.entry_mabase.get()
        self.LocalDirectory = self.entry_directory.get()
        self.DATAFolderName = self.entry_folder.get()

        try:
            conn = psycopg2.connect(database=self.maBase)

            print('Connexion à la base réussie !')
            conn.close()

        except : 
            self.LabErrorInit.config(text = 'Impossible de trouver la base ' + self.maBase,foreground = "red" )
            return
        
        if os.system("test -d " + self.LocalDirectory + self.DATAFolderName) == 0 :
            print('Directory fonctionnel !')

            self.root.destroy()
            GUI = application(self.maBase, self.LocalDirectory, self.DATAFolderName, self.fileCSV)
            GUI.interface()

        else :
            print('ERROR : Wrong Directory Path')
            self.LabErrorInit.config(text = "L'adresse " + self.LocalDirectory + self.DATAFolderName +" ne fonctionne pas",foreground = "red" )




        


class application :
    def __init__(self,DB,dir,fold,CSV)  :

        self.maBase = DB
        self.LocalDirectory = dir
        self.DATAFolderName = fold
        self.fileCSV = CSV

        conn = psycopg2.connect(database = self.maBase)
        cur = conn.cursor()
        sql_cmd = "SELECT DISTINCT name_specie FROM all_proteomes ;"
        cur.execute(sql_cmd)

        tuples = cur.fetchall()[:]

        self.all_esp = [truc[0] for truc in tuples]
        cur.close()

        conn.close()


        fichier = open(self.fileCSV,'r') 

        self.all_DL_esp = []

        for line in fichier :

            if line[0] != '#' :

                split = line.split(',')
    
                if split[-2] != '' and len(split) == 16 :
                    self.all_DL_esp.append(split[0][1:-1])

        self.all_DL_esp = np.array(self.all_DL_esp, dtype = str)
        self.all_DL_esp = np.unique(self.all_DL_esp)
        self.all_DL_esp = [i for i in self.all_DL_esp]

        pass
    
    #----------------------------------------------------
    #fonction pour la combobox 1

    def on_keyrelease1(self, event):
        # Récupérer le texte actuel dans la combobox
        text = event.widget.get()
        
        filtered_data = [item for item in self.all_esp if text.lower() in item.lower()]
        # Mettre à jour les valeurs dans la combobox
        self.combo1['values'] = filtered_data

    #----------------------------------------------------
    #fonction pour la combobox 2

    def on_keyrelease2(self, event):
        # Récupérer le texte actuel dans la combobox
        text = event.widget.get()
        
        filtered_data = [item for item in self.all_esp if text.lower() in item.lower()]
        # Mettre à jour les valeurs dans la combobox
        self.combo2['values'] = filtered_data

    #----------------------------------------------------
    #fonction pour la combobox 3

    def on_keyrelease3(self, event):
        # Récupérer le texte actuel dans la combobox
        text = event.widget.get()
        
        filtered_data = [item for item in self.all_DL_esp if text.lower() in item.lower()]
        # Mettre à jour les valeurs dans la combobox
        self.combo3['values'] = filtered_data

    #----------------------------------------------------
    #fonction pour la combobox 4

    def on_keyrelease4(self, event):
        # Récupérer le texte actuel dans la combobox
        text = event.widget.get()
        
        filtered_data = [item for item in self.all_DL_esp if text.lower() in item.lower()]
        # Mettre à jour les valeurs dans la combobox
        self.combo4['values'] = filtered_data

    #----------------------------------------------------
    #fonction pour les CheckButton
    
    def toggle_entry_state(self, nom_entry) :
        ''' attention l'attribut doit être passé en string'''
        if getattr(self,'bool_' + nom_entry).get() :
            getattr(self, 'entry_' + nom_entry).configure(state = 'normal')
        else :
            getattr(self, 'entry_' + nom_entry).configure(state = 'disabled', text = '')
	
    def interface(self) :

        #----------------------------------------------------
        #Creation root

        self.fenetre = tk.Tk()
        self.fenetre.configure(bg = '#2B2B2B')
        self.fenetre.title('Interface Projet')
        self.fenetre.geometry("620x400")

        #----------------------------------------------------
        #Creation Notebook

        # Initialize style
        s = ttk.Style()
        # Create style used by default for all Frames
        # s.configure('TFrame', background='grey')
        s.configure('TNotebook', background='#8F8F8F')
 

        self.notebook = ttk.Notebook(self.fenetre)
        self.notebook.grid(row = 0, column = 0, sticky = 'NSEW')

        self.frame1 = ttk.Frame(self.notebook, width=600, height=380)
        self.frame2 = ttk.Frame(self.notebook, width=600, height=380)
        #----------------------------------------------------
        #bloc combobox et titre au dessus
        
        N = 10
        for i in range(N):
            self.frame1.columnconfigure(i, weight = 1)
            self.frame1.rowconfigure(i, weight = 1)


        self.labE1 = tk.Label(self.frame1,  text ='Espèce 1')
        self.labE2 = tk.Label(self.frame1, text ='Espèce 2')
        self.labE1.grid(row=0,column = 0)
        self.labE2.grid(row=1,column = 0)
        
        

        self.combo1 = ttk.Combobox(self.frame1, values = self.all_esp, width = 40)
        self.combo2 = ttk.Combobox(self.frame1, values=self.all_esp, width = 40)

        self.combo1.bind('<KeyRelease>', self.on_keyrelease1)
        self.combo2.bind('<KeyRelease>', self.on_keyrelease2)

        self.combo1.grid(row = 0, column = 1, columnspan = 4)
        self.combo2.grid(row = 1, column = 1, columnspan = 4)     

        #----------------------------------------------------
        #création du bouton qui va lancer le dotplot

        self.plotButton = tk.Button(self.frame1, text = "DOTPLOT", command = self.dotplot)
        self.plotButton.grid(row=3,column = 2, columnspan = 1)

        #----------------------------------------------------
        #création d'une fenetre commentaires sous mes combobox

        self.LabError = tk.Label(self.frame1, text = '')
        self.LabError.grid(row=2,column = 0,columnspan = 5)

        #----------------------------------------------------
        #création des cases pour les critères d'homologies

        self.entry_evalue = tk.Entry(self.frame1, state = 'disabled', width = 10)
        self.bool_evalue = tk.BooleanVar()
        self.chk_evalue = tk.Checkbutton(self.frame1, text = 'E-value', variable = self.bool_evalue, command = lambda : self.toggle_entry_state('evalue') , )

        self.entry_pident = tk.Entry(self.frame1, state = 'disabled', width = 10)
        self.bool_pident = tk.BooleanVar()
        self.chk_pident = tk.Checkbutton(self.frame1, text = 'P-ident', variable = self.bool_pident, command = lambda : self.toggle_entry_state('pident') )

        self.entry_qcover = tk.Entry(self.frame1, state = 'disabled', width = 10)
        self.bool_qcover = tk.BooleanVar()
        self.chk_qcover = tk.Checkbutton(self.frame1, text = 'Q-cover', variable = self.bool_qcover, command = lambda : self.toggle_entry_state('qcover') )

        self.entry_scover = tk.Entry(self.frame1, state = 'disabled', width = 10)
        self.bool_scover = tk.BooleanVar()
        self.chk_scover = tk.Checkbutton(self.frame1, text = 'S-cover', variable = self.bool_scover, command = lambda : self.toggle_entry_state('scover') )

        self.chk_evalue.grid(row = 4, column = 0,columnspan = 1)
        self.entry_evalue.grid(row =4, column = 1)

        self.chk_pident.grid(row = 5, column = 0,columnspan = 1)
        self.entry_pident.grid(row =5, column = 1)

        self.chk_qcover.grid(row = 6, column = 0,columnspan = 1)
        self.entry_qcover.grid(row =6, column = 1)

        self.chk_scover.grid(row = 7, column = 0,columnspan = 1)
        self.entry_scover.grid(row =7, column = 1)

        #----------------------------------------------------
        #création d'une ligne séparatrice

        self.sep_line = tk.Canvas(self.frame1, width=2,height = 250)
        self.sep_line.grid(row = 4,column = 2, columnspan = 1, rowspan = 4)
        self.sep_line.create_line(1,0,1,250,fill = "black", width = 2)

        #----------------------------------------------------
        #création des entrées pour la fenetre et la stringence

        self.entry_winsize = tk.Entry(self.frame1, width = 10)

        self.entry_stringence = tk.Entry(self.frame1, width = 10)

        self.lab_winsize = tk.Label(self.frame1, text = "Window Size :")

        self.lab_stringence = tk.Label(self.frame1, text = "Stringence :")

        self.lab_winsize.grid(row = 4, column = 3,columnspan = 1)

        self.lab_stringence.grid(row = 5, column = 3,columnspan = 1)

        self.entry_winsize.grid(row = 4, column = 4,columnspan = 1)

        self.entry_stringence.grid(row = 5, column = 4,columnspan = 1)

        #----------------------------------------------------
        #Création de la 2eme frame pour remplir la database

        N = 10
        for i in range(N):
            self.frame2.rowconfigure(i, weight = 1)
            self.frame2.columnconfigure(i, weight = 1)

        #----------------------------------------------------
        #Labels et Combo
        
        self.headlab = tk.Label(self.frame2, text = 'Veuillez choisir les deux espèces sur lesquelles vous effectuer un blast')
        self.lab1 = tk.Label(self.frame2,  text ='Espèce 1')
        self.lab2 = tk.Label(self.frame2, text ='Espèce 2')
        self.headlab.grid(row=0,column = 0, columnspan = 5)
        self.lab1.grid(row=1,column = 0, columnspan = 2)
        self.lab2.grid(row=2,column = 0, columnspan = 2)

        self.combo3 = ttk.Combobox(self.frame2, values = self.all_DL_esp, width = 40) 
        self.combo4 = ttk.Combobox(self.frame2, values= self.all_DL_esp, width = 40)

        self.combo3.bind('<KeyRelease>', self.on_keyrelease3)
        self.combo4.bind('<KeyRelease>', self.on_keyrelease4)

        self.combo3.grid(row = 1, column = 2, columnspan = 3)
        self.combo4.grid(row = 2, column = 2, columnspan = 3)     

        #----------------------------------------------------
        #création du bouton qui va remplir la DB

        self.fillButton = tk.Button(self.frame2, text = "Lancer Recherche", command = self.Fill_DB)
        self.fillButton.grid(row=5,column = 0, columnspan = N-1)

        #----------------------------------------------------
        #Label error 

        self.LabError2 = tk.Label(self.frame2, text = '')
        self.LabError2.grid(row=4,column = 0,columnspan = N-1 )

        #----------------------------------------------------
        #Création de la 3eme Frame 'Options'

        #----------------------------------------------------
        #pack des frames du notebook

        self.frame1.grid(row = 0, column = 0)
        self.frame2.grid(row = 0, column = 0)

        self.notebook.add(self.frame1, text='DOTPLOT')
        self.notebook.add(self.frame2, text='FILL DB')

        self.fenetre.mainloop()

    def dotplot(self) :

        self.LabError.config(text='')

        self.espy = self.combo1.get()
        self.espx = self.combo2.get()

        # self.espy = 'Escherichia coli IAI1'
        # self.espx = 'Escherichia coli str. K-12 substr. MG1655'

        if (self.espy == '' or self.espx == '' ) :
            self.LabError.config(text="Veuillez sélectioner deux espèces", foreground = "red")
            return
        elif self.espy == self.espx :
            self.LabError.config(text="Veuillez sélectioner des espèces différentes", foreground = "red")
            return

        
        
        else :             
            

            #l'espece 1 correspond à l'axe 0 et l'esp2 à l'axe 1
            dicoSpecie = {}
            dicoSpecie[self.espy] = 0
            dicoSpecie[self.espx] = 1
        

            sql_cmd = "SELECT " 
            
            nb_crit = 0
            array_test = []

            if self.bool_evalue.get(): 
                sql_cmd += "e_value "
                array_test.append(float(self.entry_evalue.get()))
                nb_crit += 1
            if self.bool_pident.get() :
                sql_cmd += "pidentite "
                array_test.append(float(self.entry_pident.get()))
                nb_crit += 1
            if self.bool_qcover.get() :
                sql_cmd += "query_cover "
                array_test.append(float(self.entry_qcover.get()))
                nb_crit += 1
            if self.bool_scover.get() :
                sql_cmd += "subject_cover "
                array_test.append(float(self.entry_scover.get()))
                nb_crit += 1
            
            split = sql_cmd.split()

            if len(split) == 1 :
                self.LabError.config(text="Veuillez cocher au moins 1 critère", foreground = "red")
                return
            
            else :
                
                sql_cmd = "SELECT " + split[1]

                for criter in split[2:] :
                    sql_cmd += ", "+ criter 

                sql_cmd += ", query_specie, query_rang, subject_rang FROM blast WHERE (query_specie = %s AND subject_specie = %s) OR (query_specie = %s AND subject_specie = %s) ;"
                
                conn = psycopg2.connect(database = self.maBase)
                cur = conn.cursor()
                cur.execute(sql_cmd, (self.espy,self.espx,self.espx,self.espy))
                data = np.array(cur.fetchall())
                
                winsize = self.entry_winsize.get()
                stringence = self.entry_stringence.get()

                if len(data) == 0 :
                    self.LabError.config(text="ERROR : pas de données pour cette combinaison", foreground = "red")
                    return

                elif winsize == '' or stringence == '' :
                    self.LabError.config(text="Veuillez sélectioner une taille de fenêtre et une stringence", foreground = "red")
                    return
        

                else :
                    #data = 'eval, pid, qcov, scov, qspec, qrg,srg'
                    
                    winsize = int(self.entry_winsize.get())
                    stringence = int(self.entry_stringence.get())

                    array_test = np.array(array_test)

                    data_test = data[:,:nb_crit].astype(np.float64)
             

                    if self.bool_evalue.get(): 
                        signAjust = np.ones(nb_crit)
                        signAjust[0] = -1
                        rows = np.where(np.amin((data_test - array_test) * signAjust , axis = 1) > 0)[0]
                    else :   
                        rows = np.where(np.amin((data_test - array_test), axis = 1) > 0)[0] 

                   

                    cur.execute("SELECT COUNT(*) FROM all_proteomes WHERE name_specie = %s", (self.espy,))
                    lenE1 = cur.fetchall()[0][0]
                    cur.execute("SELECT COUNT(*) FROM all_proteomes WHERE name_specie = %s", (self.espx,))
                    lenE2 = cur.fetchall()[0][0]


                    self.dotplot = np.zeros((lenE1,lenE2))

                    #query_specie, query_rang, subject_rang

                    #Version avec dico qui considère qu'il peut y avoir esp 1 en query et subject
                    # for line in data[rows] :

                    #     if dicoSpecie[line[-3]] == 0 :
                            
                    #         self.dotplot[int(line[-2])-1 , int(line[-1])-1] = 1
                          
                    #     else :
                    #         self.dotplot[int(line[-1])-1  , int(line[-2])-1] = 1

                
                    # self.dotplotbis = 
                    if dicoSpecie[data[0,-3]] == 0 :
                        self.dotplot[(data[rows][:,-2].astype(int) -1),(data[rows][:,-1].astype(int) -1)] = 1
                    else :
                        
                        self.dotplot[(data[rows][:,-1].astype(int) -1),(data[rows][:,-2].astype(int) -1)] = 1

                    newdot = np.copy(self.dotplot) 

                    shape = newdot.shape

                    for i in range(1,winsize) :
                        newdot[:-i, :-i] += self.dotplot[i:, i:]

                   

                    self.x , self.y = np.where(newdot >= stringence) ## attention verifier si int ou float
                    print(len(self.x), ' Hits')

                    # newdot = np.where(newdot >= stringence, 1, 0)
                    # self.dotplot = newdot
             


                    self.plot_dotplot()

                cur.close()
                conn.close()  

    def plot_dotplot(self) :

        self.canva = tk.Tk()
        plt.rcParams['axes.facecolor'] = 'darkslategrey'

        fig = Figure(figsize=(10,10), dpi=100)
        ax = fig.add_subplot(111)        

        ax.set_xlabel('{}'.format(self.espx))
        ax.set_ylabel('{}'.format(self.espy))

        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        

        ax.invert_yaxis()

        ax.scatter(self.x,self.y, s = 1, c = 'yellow' , marker = ',')
        # ax.imshow(self.dotplot[:, :], interpolation = None, cmap = 'viridis')

        

        canvas = FigureCanvasTkAgg(fig, master=self.canva)
        # canvas.print_figure('ma_figure.png', dpi=100)

        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.canva.mainloop() 



        #idée : remplir dictionnaire dico[(espèce,rang)] : protid
        #puis parcourir le blast 

    def Fill_DB(self) :

        self.LabError2.config(text="")

        os.system('rm ' + self.LocalDirectory + self.DATAFolderName + '*')

      
        self.DL_esp1 = self.combo3.get()
        self.DL_esp2 = self.combo4.get()


        self.LabError2.config(text="", foreground = "red")

        if self.DL_esp1 == '' or  self.DL_esp2 == '':
            self.LabError2.config(text="Veuillez sélectioner deux espèces", foreground = "red")
            return

        elif self.DL_esp1 == self.DL_esp2 :
            self.LabError2.config(text="Veuillez sélectioner deux espèces différentes", foreground = "red")
            return
        

        conn = psycopg2.connect(database = self.maBase)
        cur = conn.cursor()

        sql_cmd = "SELECT * FROM all_proteomes WHERE name_specie = %s LIMIT 1;"

        cur.execute(sql_cmd, (self.DL_esp1 ,))
        isHere1 = len(cur.fetchall())
        cur.execute(sql_cmd, (self.DL_esp2 ,))
        isHere2 = len(cur.fetchall())

        #----------------------------------------------------
        #aller télécharger grace au fichier self.fileCSV
        #Si espèce non présente : remplir DB
        #Attention, lorsque plusieurs fois la meme espèce, télécharge données plus récente

            
        fichier = open(self.fileCSV, 'r')

        search = []

        for line in fichier :
            split = line.split(',')

            if split[0][1:-1] == self.DL_esp1 :
                date = split[13][1:-11].split('-')
                date = date[0] + date[1] + date[2]
                   
                search.append((int(date), split[-2][1:-1]))

        if len(search) > 1 :
               
            d = 0
            for i in search :
                if i[0] > d :
                    link = i[1]
            
        else :
            link = search[0][1]

        name_file =  link.split('/')[-1]+ '_translated_cds.faa.gz'
        link = link + '/' + name_file

        os.system('wget -O ' + self.LocalDirectory + self.DATAFolderName + 'cds_query.faa.gz ' + link)
        os.system('gunzip ' + self.LocalDirectory + self.DATAFolderName + 'cds_query.faa.gz')


        if isHere1 == 0 :
            fill_prots_db( self.LocalDirectory + self.DATAFolderName + 'cds_query.faa' , self.DL_esp1, self.maBase) 

        
            sql_cmd = "SELECT DISTINCT name_specie FROM all_proteomes ;"
            cur.execute(sql_cmd)
            tuples = cur.fetchall()[:]
            self.all_esp = [truc[0] for truc in tuples]
            self.combo1['values'] = self.all_esp
            self.combo2['values'] = self.all_esp
         

        fichier.close()
        #----------------------------------------------------
        #Fichier Subject

        fichier = open(self.fileCSV, 'r')
            
        search = []

        for line in fichier :
            split = line.split(',')

            if split[0][1:-1] == self.DL_esp2 :

                date = split[13][1:-11].split('-')
                date = date[0] + date[1] + date[2]
                   
                search.append((int(date), split[-2][1:-1]))

        fichier.close()

        if len(search) > 1 :
            d = 0
            for i in search :
                if i[0] > d :
                    link = i[1]
            
        else :
            print('search : ', search)
            link = search[0][1]

        name_file =  link.split('/')[-1]+ '_translated_cds.faa.gz'
        link = link + '/' + name_file
        os.system('wget -O ' + self.LocalDirectory + self.DATAFolderName + 'cds_subject.faa.gz ' + link)
        os.system('gunzip ' + self.LocalDirectory + self.DATAFolderName + 'cds_subject.faa.gz')

        if isHere2 == 0 :
            fill_prots_db(  self.LocalDirectory + self.DATAFolderName + 'cds_subject.faa' , self.DL_esp2, self.maBase) 

       
            sql_cmd = "SELECT DISTINCT name_specie FROM all_proteomes ;"
            cur.execute(sql_cmd)
            tuples = cur.fetchall()[:]
            self.all_esp = [truc[0] for truc in tuples]
           
            self.combo1['values'] = self.all_esp
            self.combo2['values'] = self.all_esp

                    

        #----------------------------------------------------
        #vérification que le blast n'existe pas déja

        sql_cmd = 'SELECT * FROM blast WHERE (query_specie = %s AND subject_specie = %s) OR (query_specie = %s AND subject_specie = %s) LIMIT 1'

        cur.execute(sql_cmd, (self.DL_esp1,self.DL_esp2,self.DL_esp2,self.DL_esp1))
        isBlastHere = len(cur.fetchall())

        if isBlastHere > 0 :
            self.LabError2.config(text="Un Blast existe déja pour ces espèces", foreground = "red")

        else : 
            os.system('blastp -query ' + self.LocalDirectory + self.DATAFolderName + 'cds_query.faa -subject ' + self.LocalDirectory + self.DATAFolderName + 'cds_subject.faa -evalue 1e-10 -out ' + self.LocalDirectory + self.DATAFolderName + 'BlastResult.out -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore"')
            fill_blast_db(self.LocalDirectory + self.DATAFolderName + 'BlastResult.out', self.maBase)
            
            print('BLAST terminé')
            self.LabError2.config(text="BLAST terminé", foreground = "black")
            


        cur.close()
        conn.close()



Configure = ConfigureWindow(maBase, LocalDirectory, DATAFolderName, fileCSV)


#Choses à améliorer :

#Penser retirer self.dotplot

#NB combo faire en sorte que seuleemnt les espèces avec blast ressortent

