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

#establish vars to hold file information
currDir = os.getcwd()
outputFolderName = 'Output'    
try:
    os.mkdir(outputFolderName)
except FileExistsError:
    newName = input("\nFile with the name 'Output' already exists in this directory.\nEnter new name for output file: ")
    outputFolderName = newName
    os.mkdir(outputFolderName)



#reading file input and then closing file
#cast file input to list
with open(args.input) as f:
    genes = f.read().splitlines()

capitalGenes = [gene.upper() for gene in genes]    
f.close()
print (capitalGenes)


#this method connects to the database, specifically the hg19_ensembl table
#it grabs the genes out of the table that are input file, along with a few other  params
#this returns a dataframe of genes with their name, chr #, and adj bp start/end locations
#it also creates two txt files, one containg all genes from input that are searched by the 
#script, and all the genes that are not.
def hg19Search(margin):
    geneList = []
    i = 0
    geneDF = pd.DataFrame(columns=["gene_name", "chr", "temp_start", "temp_end"])
    while i < len(genes):
        queryString = "SELECT gene_name, chr, MIN(chr_start) - %d AS temp_start, MAX(chr_end) + %d AS temp_end FROM hg19_ensembl WHERE gene_name = UPPER(\'%s\') COLLATE SQL_Latin1_General_CP1_CS_AS GROUP BY gene_name, chr" %(margin, margin, genes[i])
        geneList.append(queryString)
        geneFrame = pd.read_sql(geneList[i], conn)
        geneDF = geneDF.append(geneFrame)
        i += 1 
    xGenes = geneDF['chr']=='X'
    xframe = geneDF[xGenes]
    xlist = list(xframe['gene_name'])
    geneDF = geneDF[geneDF.chr != 'X']
    colGenes = list(geneDF['gene_name'])
    errList = list(set(capitalGenes) - set(colGenes) - set(xlist))

    #creation of the file holding all genes not searched
    genesNotSearched = open(currDir + '/' + outputFolderName + "/genesNotSearched.txt", 'w')
    genesNotSearched.write("Genes in input file and in X chr: \n")
    for x in xlist:
        genesNotSearched.write(x)
        genesNotSearched.write('\n')
    genesNotSearched.write("Genes in input file and not in database: \n")
    for err in errList:
        genesNotSearched.write(err)
        genesNotSearched.write('\n')
    genesNotSearched.close()

    #creation of the file holding all genes that are searched
    valGenFile = open(currDir + "/" + outputFolderName + "/validGeneFile.txt", "w")
    for gene in colGenes:
        valGenFile.write(gene)
        valGenFile.write('\n')
    valGenFile.close()
    return(geneDF)

#this methods builds the string for the "Where" segement of the sql search. 
#the method takes the parameters given to if from user, and then builds a custom
#string for each gene to be passed through each table
def sqlSearchString(geneDF, pval):
    stringList = []
    geneCount = len(geneDF)
    i = 0
    while i < geneCount:
            string = "chr = \'%d\' AND bp BETWEEN %d AND %d AND p_value < \'%f\'"  %(int(geneDF.iloc[i,1]), int(geneDF.iloc[i,2]), int(geneDF.iloc[i,3]), pval)
            stringList.append(string)
            i += 1
    return(stringList)

#This method parses through the gene_name column of the geneDataFrame
#it returns this column in a list format
def geneName(geneDF):
    i = 0
    geneCount = len(geneDF)
    geneNameList = []
    while i < geneCount:
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
            bivariateHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Bivariate_Siewert_CGPM_HDL' AS study_name FROM Lipid_Bivariate_Siewert_CGPM_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(bivariateHD)
            bivariateLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Bivariate_Siewert_CGPM_LDL' AS study_name FROM Lipid_Bivariate_Siewert_CGPM_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(bivariateLD)
            bivariateTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Bivariate_Siewert_CGPM_TG' AS study_name FROM Lipid_Bivariate_Siewert_CGPM_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(bivariateTG)
            bivariateTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Bivariate_Siewert_CGPM_TC' AS study_name FROM Lipid_Bivariate_Siewert_CGPM_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(bivariateTC)
            chineseHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Chinese_exome_Lu_HMG_HDL' AS study_name FROM Lipid_Chinese_exome_Lu_HMG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC" 
            views.append(chineseHD)
            chineseLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Chinese_exome_Lu_HMG_LDL' AS study_name FROM Lipid_Chinese_exome_Lu_HMG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC" 
            views.append(chineseLD)
            chineseTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Chinese_exome_Lu_HMG_TG' AS study_name FROM Lipid_Chinese_exome_Lu_HMG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC" 
            views.append(chineseTG)
            chineseTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Chinese_exome_Lu_HMG_TC' AS study_name FROM Lipid_Chinese_exome_Lu_HMG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC" 
            views.append(chineseTC)
            ehrHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_EHR_Hoffmann_NG_HDL' AS study_name FROM Lipid_EHR_Hoffmann_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ehrHD)
            ehrLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_EHR_Hoffmann_NG_LDL' AS study_name FROM Lipid_EHR_Hoffmann_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ehrLD)
            ehrTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_EHR_Hoffmann_NG_TC' AS study_name FROM Lipid_EHR_Hoffmann_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ehrTC)
            ehrTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_EHR_Hoffmann_NG_TG' AS study_name FROM Lipid_EHR_Hoffmann_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ehrTG)
            surakkaHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Engage_Surakka_NG_HDL' AS study_name FROM Lipid_Engage_Surakka_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(surakkaHD)
            surakkaLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Engage_Surakka_NG_LDL' AS study_name FROM Lipid_Engage_Surakka_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(surakkaLD)
            surakkaTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Engage_Surakka_NG_TG' AS study_name FROM Lipid_Engage_Surakka_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(surakkaTG)
            surakkaTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Engage_Surakka_NG_TC' AS study_name FROM Lipid_Engage_Surakka_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(surakkaTC)
            eastAsianHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_East_Asian_NG_HDL' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(eastAsianHD)
            eastAsianLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_East_Asian_NG_LDL' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(eastAsianLD)
            eastAsianTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_East_Asian_NG_TG' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(eastAsianTG)
            eastAsianTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_East_Asian_NG_TC' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(eastAsianTC)
            euroAsianHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_European_and_East_Asian_NG_HDL' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(euroAsianHD)
            euroAsianLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_European_and_East_Asian_NG_LDL' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(euroAsianLD)
            euroAsianTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_European_and_East_Asian_NG_TG' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(euroAsianTG)
            euroAsianTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Exome_Lu_European_and_East_Asian_NG_TC' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(euroAsianTC)
            finnishHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Finnish_Locke_HDL' AS study_name FROM Lipid_Finnish_Locke_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(finnishHD)
            finnishLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Finnish_Locke_LDL' AS study_name FROM Lipid_Finnish_Locke_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(finnishLD)
            finnishTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Finnish_Locke_TG' AS study_name FROM Lipid_Finnish_Locke_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(finnishTG)
            finnishTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Finnish_Locke_TC' AS study_name FROM Lipid_Finnish_Locke_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(finnishTC)
            glgcExomeHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Exome_Liu_NG_HDL' AS study_name FROM Lipid_GLGC_Exome_Liu_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcExomeHD)
            glgcExomeLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Exome_Liu_NG_LDL' AS study_name FROM Lipid_GLGC_Exome_Liu_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcExomeLD)
            glgcExomeTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Exome_Liu_NG_TG' AS study_name FROM Lipid_GLGC_Exome_Liu_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcExomeTG)
            glgcExomeTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Exome_Liu_NG_TC' AS study_name FROM Lipid_GLGC_Exome_Liu_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcExomeTC)
            glgcWillerHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Willer_NG_HDL' AS study_name FROM Lipid_GLGC_Willer_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcWillerHD)
            glgcWillerLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Willer_NG_LDL' AS study_name FROM Lipid_GLGC_Willer_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcWillerLD)
            glgcWillerTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Willer_NG_TG' AS study_name FROM Lipid_GLGC_Willer_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcWillerTG)
            glgcWillerTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_GLGC_Willer_NG_TC' AS study_name FROM Lipid_GLGC_Willer_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(glgcWillerTC)
            japaneseHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Japanese_lipid_trait_Kanai_NG_HDL' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(japaneseHD)
            japaneseLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Japanese_lipid_trait_Kanai_NG_LDL' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(japaneseLD)
            japaneseTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Japanese_lipid_trait_Kanai_NG_TG' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(japaneseTG)
            japaneseTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Japanese_lipid_trait_Kanai_NG_TC' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(japaneseTC)
            koreanHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Korean_Exome_Moon_SR_HDL' AS study_name FROM Lipid_Korean_Exome_Moon_SR_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(koreanHD)
            koreanLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Korean_Exome_Moon_SR_LDL' AS study_name FROM Lipid_Korean_Exome_Moon_SR_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(koreanLD)
            koreanTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Korean_Exome_Moon_SR_TG' AS study_name FROM Lipid_Korean_Exome_Moon_SR_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(koreanTG)
            koreanTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Korean_Exome_Moon_SR_TC' AS study_name FROM Lipid_Korean_Exome_Moon_SR_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(koreanTC)
            mvpHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_MVP_Klarin_NG_HDL' AS study_name FROM Lipid_MVP_Klarin_NG_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(mvpHD)
            mvpLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_MVP_Klarin_NG_LDL' AS study_name FROM Lipid_MVP_Klarin_NG_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(mvpLD)
            mvpTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_MVP_Klarin_NG_TG' AS study_name FROM Lipid_MVP_Klarin_NG_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(mvpTG)
            mvpTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_MVP_Klarin_NG_TC' AS study_name FROM Lipid_MVP_Klarin_NG_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(mvpTC)
            spracklenHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Spracklen_Hum_Mol_Genetics_HDL' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(spracklenHD)
            spracklenLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Spracklen_Hum_Mol_Genetics_LDL' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(spracklenLD)
            spracklenTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Spracklen_Hum_Mol_Genetics_TG' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(spracklenTG)
            spracklenTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Spracklen_Hum_Mol_Genetics_TC' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(spracklenTC)
            teslovichHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Teslovich_Nature_HDL' AS study_name FROM Lipid_Teslovich_Nature_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(teslovichHD)
            teslovichLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Teslovich_Nature_LDL' AS study_name FROM Lipid_Teslovich_Nature_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(teslovichLD)
            teslovichTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Teslovich_Nature_TG' AS study_name FROM Lipid_Teslovich_Nature_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(teslovichTG)
            teslovichTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_Teslovich_Nature_TC' AS study_name FROM Lipid_Teslovich_Nature_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(teslovichTC)
            uk10kHD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_UK10K_prins_SR_HDL' AS study_name FROM Lipid_UK10K_prins_SR_HDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(uk10kHD)
            uk10kLD = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_UK10K_prins_SR_LDL' AS study_name FROM Lipid_UK10K_prins_SR_LDL WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(uk10kLD)
            uk10kTG = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_UK10K_prins_SR_TG' AS study_name FROM Lipid_UK10K_prins_SR_TG WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(uk10kTG)
            uk10kTC = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_UK10K_prins_TC' AS study_name FROM Lipid_UK10K_prins_SR_TC WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(uk10kTC)
            ukbbHigh = "SELECT rsid, chr, bp, p_value, beta, 'Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup' AS study_name FROM Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ukbbHigh)
            ukbbLipid = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_lipid_trait_Neal' AS study_name FROM Lipid_UKBB_lipid_trait_Neale WHERE " + stringList[j] + " ORDER BY p_value ASC"
            views.append(ukbbLipid)

        
            # #create DataFrame from sql connection string
            i = 0
            #list to hold the temp dfs for each table, will then get concat
            frames = []
            viewLen = len(views)
            while i < viewLen:
                tempdf = pd.read_sql(views[i], conn)
                frames.append(tempdf)
                i += 1
            df = pd.concat(frames, ignore_index=True, sort=True) 
            df.insert(0, 'gene_name', geneNameList[j]) 
            path =  currDir + '/' + outputFolderName + '/%s.txt' %(geneNameList[j])
            df.to_csv(path, encoding = 'utf-8', sep='\t',index = False)
            df = df.iloc[0:0]
            j += 1  

if __name__ == "__main__":
    tableLoad(sqlSearchString(hg19Search(args.margin), args.pval), geneName(hg19Search(args.pval)))