import sys
import pyodbc
import pandas as pd
import argparse
import numpy as np
import os

#server connection
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

#set up for command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--margin', type=int, metavar = '',  help='Base Pair margin interval', default=200000)
parser.add_argument('-p','--pval', type=float, metavar = '', help='Cutoff value of p-value', default=0.00005)
parser.add_argument('-i', '--input', type=str, metavar = '', help= 'Input file containing genes to run through script', required=True)
args = parser.parse_args()
pval = args.pval

#reading file input and then closing file
#cast file input to list
f = open(args.input,"r")
genes = f.read()
f.close()
genes = str(genes).split(",")
print(genes)

#this method connects to the database, specifically the hg19_ensembl table
#it grabs the genes out of the table that are input file, along with a few other  params
#this returns a dataframe of genes with their name, chr #, and adj bp start/end locations
def hg19Search(margin):
    geneList = []
    i = 0
    geneDF = pd.DataFrame(columns=["gene_name", "chr", "temp_start", "temp_end"])
    while i < len(genes):
        queryString = "SELECT gene_name, chr, MIN(chr_start) - %d AS temp_start, MAX(chr_end) + %d AS temp_end FROM HG19_ensembl WHERE gene_name = UPPER(\'%s\') COLLATE SQL_Latin1_General_CP1_CS_AS GROUP BY gene_name, chr" %(margin, margin, genes[i])
        geneList.append(queryString)
        geneFrame = pd.read_sql(geneList[i], conn)
        geneDF = geneDF.append(geneFrame)
        i += 1 
    return(geneDF)

#this methods builds the string for the "Where" segement of the sql search. 
#the method takes the parameters given to if from user, and then builds a custom
#string for each gene to be passed through each table
def sqlSearchString(geneDF, pval):
    stringList = []
    i = 0
    while i < len(genes):
        string = "chr = \'%d\' AND bp BETWEEN %d AND %d AND p_value < \'%f\'"  %(int(geneDF.iloc[i,1]), int(geneDF.iloc[i,2]), int(geneDF.iloc[i,3]), pval)
        stringList.append(string)
        i += 1
    return(stringList)

#This method parses through the gene_name column of the geneDataFrame
#it returns this column in a list format
def geneName(geneDF):
    i = 0
    geneNameList = []
    while i < len(genes):
        geneName = geneDF.iloc[i,0]
        geneNameList.append(geneName)
        i += 1
    return(geneNameList)

def tableLoad(stringList, geneNameList):
    #init counter vars
    j = 0
    x = 0
    while j < len(stringList):
            #list to hold connection strings
            views = []
            #SQL connection strings to get looped thru and then populate dataFrame
            bivariate = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Bivariate_Siewert_CGPM' AS study_name FROM Lipid_Bivariate_Siewert_CGPM WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(bivariate)
            chinese = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Chinese_exome_Lu_HMG' AS study_name FROM Lipid_Chinese_exome_Lu_HMG WHERE " + stringList[j] + " ORDER BY p_value ASC" 
            views.append(chinese)
            ehr = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_EHR_Hoffmann_NG' AS study_name FROM Lipid_EHR_Hoffmann_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ehr)
            surakka = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Engage_Surakka_NG' AS study_name FROM Lipid_Engage_Surakka_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(surakka)
            eastAsian = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Exome_Lu_East_Asian_NG' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(eastAsian)
            euroAsian = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Exome_Lu_European_and_East_Asian_NG' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(euroAsian)
            finnish = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Finnish_Locke' AS study_name FROM Lipid_Finnish_Locke WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(finnish)
            glgcExome = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_GLGC_Exome_Liu_NG' AS study_name FROM Lipid_GLGC_Exome_Liu_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcExome)
            glgcWiller = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_GLGC_Willer_NG' AS study_name FROM Lipid_GLGC_Willer_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcWiller)
            japanese = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Japanese_lipid_trait_Kanai_NG' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(japanese)
            korean = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Korean_Exome_Moon_SR' AS study_name FROM Lipid_Korean_Exome_Moon_SR WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(korean)
            mvp = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_MVP_Klarin_NG' AS study_name FROM Lipid_MVP_Klarin_NG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(mvp)
            spracklen = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Spracklen_Hum_Mol_Genetics' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(spracklen)
            teslovich = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Teslovich_Nature' AS study_name FROM Lipid_Teslovich_Nature WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(teslovich)
            uk10k = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UK10K_prins_SR' AS study_name FROM Lipid_UK10K_prins_SR WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(uk10k)
            ukbbHigh = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup' AS study_name FROM Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ukbbHigh)
            ukbbLipid = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_lipid_trait_Neal' AS study_name FROM Lipid_UKBB_lipid_trait_Neale WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ukbbLipid)
            ukbbStatin = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_statin_usage_Neale' AS study_name FROM Lipid_UKBB_statin_usage_Neale WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ukbbStatin)
        
            # #create DataFrame from sql connection string
            i = 0
            #list to hold the temp dfs for each table, will then get concat
            frames = []
            viewLen = len(views)
            while i < viewLen:
                tempdf = pd.read_sql(views[i], conn)
                frames.append(tempdf)
                i += 1
            df = pd.concat(frames) 
            df.insert(0, 'gene_name', geneNameList[j]) 
            path = 'E:/cross_ref/jgrot/output/%s.txt' %(geneNameList[j])
            df.to_csv(path, encoding = 'utf-8', sep='\t',index = False)
            df = df.iloc[0:0]
            j += 1  

if __name__ == "__main__":
    tableLoad(sqlSearchString(hg19Search(args.margin), args.pval), geneName(hg19Search(args.pval)))