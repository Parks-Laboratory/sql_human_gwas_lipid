from sqlalchemy import create_engine
import numpy as np 
import pyodbc
import pandas as pd
import os
import sys

try:
    conn =  pyodbc.connect("""
                            DRIVER={ODBC Driver 11 for SQL Server};
                            SERVER=PARKSLAB;
                            DATABASE=Human_GWAS;
                            Trusted_Connection=Yes;
                            """)
except:
    print("Unable to connect to Server :(")
    sys.exit()

def error_messeage():
    print("""
          ***********************************
          *****Please enter valid input!*****
          ***********************************
          """)

def getGene():
    while True:
            gene = input("\nHello! To get started enter gene name(s)\nOr press ENTER to input txt file.\n")
            if  gene.isdigit():
                error_messeage()
                continue
            elif gene == "":
                filepath = input("input your file path: ")
                file = open(filepath, 'r')
                gene_name_list = file.read()
                gene = gene_name_list.split(",")
                print ("input genes: ")
                print (gene)
                file.close()
                break
            else:
                gene = gene.split(",")
                break   
    return(gene)    

def getPValue():
    cutoff = input("""
Press just 'ENTER' to use the default value of 0.05.
Otherwise, specify the cutoff for p-value:\n""")
    if cutoff == '':
        return(0.05)
    while True:
        try:
            cutoff = float(cutoff)
            while (cutoff < 0) or (cutoff > 1):
                print('Enter a valid input.\n')
                cutoff = float(input("""
Press just 'ENTER' to use the default value of 0.05.\n
Otherwise, specify the cutoff for p-value:\n"""))

        except ValueError:
            print('Enter a valid input.\n')
            cutoff = input("""
Press just 'ENTER' to use the default value of 0.05.\n
Otherwise, specify the cutoff for p-value:\n""")

        return(cutoff)

def getMargin():
    while True:
        try:
            margin_input = input("""
Press just 'ENTER' to use default value of 200,000.
Otherwise, specify your margin:""")

            if margin_input == '':
                return(200000)
            else:
                margin = int(margin_input)
                break
        except ValueError:
            print('Enter a valid input.\n')
    return(margin)

def SQLViews():
    views = []
    
    bivariate = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Bivariate_Siewert_CGPM ORDER BY chr ASC"
    views.append(bivariate)
    chinese = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Chinese_exome_Lu_HMG ORDER BY chr ASC"
    views.append(chinese)
    ehr = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_EHR_Hoffmann_NG ORDER BY chr ASC"
    views.append(ehr)
    surakka = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Engage_Surakka_NG ORDER BY chr ASC"
    views.append(surakka)
    eastAsian = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Exome_Lu_East_Asian_NG ORDER BY chr ASC"
    views.append(eastAsian)
    euroAsian = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Exome_Lu_European_East_Asian_NG ORDER BY chr ASC"
    views.append(euroAsian)
    finnish = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Finnish_Locke ORDER BY chr ASC"
    views.append(finnish)
    glgcExome = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_GLGC_Exome_Liu_NG ORDER BY chr ASC"
    views.append(glgcExome)
    glgcWiller = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_GLGC_Willer_NG ORDER BY chr ASC"
    views.append(glgcWiller)
    japanese = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Japanese_lipid_trait_Kanai_NG ORDER BY chr ASC"
    views.append(japanese)
    korean = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Korean_Exome_Moon_SR ORDER BY chr ASC"
    views.append(korean)
    mvp = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_MVP_Klarin_NG ORDER BY chr ASC"
    views.append(mvp)
    spracklen = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Spracklen_Hum_Mol_Genetics ORDER BY chr ASC"
    views.append(spracklen)
    teslovich = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_Teslovich_Nature ORDER BY chr ASC"
    views.append(teslovich)
    uk10k = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_UK10K_prins_SR ORDER BY chr ASC"
    views.append(uk10k)
    ukbbHigh = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup ORDER BY chr ASC"
    views.append(ukbbHigh)
    ukbbLipid = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_UKBB_lipid_trait_Neale ORDER BY chr ASC"
    views.append(ukbbLipid)
    ukbbStatin = "SELECT rsid, chr, bp, p_value, beta FROM Lipid_UKB_statin_usage_Neale ORDER BY chr ASC"
    views.append(ukbbStatin)
   # return (views)

#conn, gene, cutoff, margin
#def dataframe():
    viewLen = len(views)
    i = 0
    while i < 5:
        df = pd.read_sql(views[i], conn) #columns=["rsid", "chr", "bp", "p_value", "beta"]
        # for i in gene:
        print(df)
        return(df)

def exitMethod():
    exitInput = int(input("\nTo exit program enter 1:\nRun a new query by entering 0:\n"))
    if exitInput == 1:
        sys.exit()
    elif exitInput == 0:
        main()
    else:
        exitMethod()

def main():
   
    while True:
        getGene()
        getMargin()
        getPValue()
        SQLViews()
        break
        #dataframe()
    exitMethod()
if __name__ == "__main__":
    main()