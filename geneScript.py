import sys
import pyodbc
import pandas as pd
import argparse

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


f = open(sys.argv[1],"r")
genes = f.read()
f.close()

def error_messeage():
    print("""
          ***********************************
          *****Please enter valid input!*****
          ***********************************
          """)

def setMargin():
    while True:
        margin = input("""\n----------SET YOUR MARGIN----------\n
Press just 'ENTER' to use default value of 200,000.
Otherwise, specify your margin:\n""")
        try:
            #200000 is default value
            if margin == '':
                margin = 200000
                break
            #cast input to int... helps create valueError
            else:
                margin = int(margin)
                break
        except ValueError:
            error_messeage()
            continue        
    return(margin)

def setPValue():
    while True:
        cutoff = input("""\n----------SET YOUR P-VALUE---------\n
Press just 'ENTER' to use the default value of 0.05.
Otherwise, specify the cutoff for p-value:\n""")
        try:
            #.05 is default value
            if cutoff == '':
                cutoff = 0.05
                break
            #catch outliers
            elif (float(cutoff) < 0.0) or (float(cutoff) > 1.0):
                error_messeage()
                continue 
            #cast input to float... helps create valueError
            else:
                cutoff = float(cutoff)  
                break    
        except ValueError:
                error_messeage()
                continue
    return(cutoff)

def hg19Search(margin):

    queryString = """SELECT gene_name, chr, MIN(chr_start) - %d AS temp_start, MAX(chr_end) + %d AS temp_end FROM HG19_ensembl
                    WHERE gene_name = UPPER(\'%s\') COLLATE SQL_Latin1_General_CP1_CS_AS GROUP BY gene_name, chr""" %(margin, margin, genes)
    geneFrame = pd.read_sql(queryString, conn)
    print(geneFrame)
    return(geneFrame)

def sqlSearchString(geneFrame, cutoff):
    string = "chr = %d AND bp BETWEEN %d AND %d AND p_value < \'%f\'"  %(geneFrame['chr'], geneFrame['temp_start'], geneFrame['temp_end'], cutoff)
    return(string)

def tableLoad(string):
     #list that will hold connection strings
    views = []

    #SQL connection strings to get looped thru and then populate dataFrame
    bivariate = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Bivariate_Siewert_CGPM' AS study_name FROM Lipid_Bivariate_Siewert_CGPM WHERE " + string + " ORDER BY p_value ASC"
    views.append(bivariate)
    chinese = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Chinese_exome_Lu_HMG' AS study_name FROM Lipid_Chinese_exome_Lu_HMG WHERE " + string + " ORDER BY p_value ASC" 
    views.append(chinese)
    ehr = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_EHR_Hoffmann_NG' AS study_name FROM Lipid_EHR_Hoffmann_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(ehr)
    surakka = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Engage_Surakka_NG' AS study_name FROM Lipid_Engage_Surakka_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(surakka)
    eastAsian = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Exome_Lu_East_Asian_NG' AS study_name FROM Lipid_Exome_Lu_East_Asian_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(eastAsian)
    euroAsian = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Exome_Lu_European_and_East_Asian_NG' AS study_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(euroAsian)
    # finnish = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Finnish_Locke' AS study_name FROM Lipid_Finnish_Locke WHERE " + string + " ORDER BY p_value ASC"
    # views.append(finnish)
    glgcExome = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_GLGC_Exome_Liu_NG' AS study_name FROM Lipid_GLGC_Exome_Liu_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(glgcExome)
    glgcWiller = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_GLGC_Willer_NG' AS study_name FROM Lipid_GLGC_Willer_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(glgcWiller)
    japanese = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Japanese_lipid_trait_Kanai_NG' AS study_name FROM Lipid_Japanese_lipid_trait_Kanai_NG WHERE " + string + " ORDER BY p_value ASC"
    views.append(japanese)
    # korean = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Korean_Exome_Moon_SR' AS study_name FROM Lipid_Korean_Exome_Moon_SR WHERE " + string + " ORDER BY p_value ASC"
    # views.append(korean)
    # mvp = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_MVP_Klarin_NG' AS study_name FROM Lipid_MVP_Klarin_NG WHERE " + string + " ORDER BY p_value ASC"
    # views.append(mvp)
    spracklen = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Spracklen_Hum_Mol_Genetics' AS study_name FROM Lipid_Spracklen_Hum_Mol_Genetics WHERE " + string + " ORDER BY p_value ASC"
    views.append(spracklen)
    teslovich = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_Teslovich_Nature' AS study_name FROM Lipid_Teslovich_Nature WHERE " + string + " ORDER BY p_value ASC"
    views.append(teslovich)
    uk10k = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UK10K_prins_SR' AS study_name FROM Lipid_UK10K_prins_SR WHERE " + string + " ORDER BY p_value ASC"
    views.append(uk10k)
    ukbbHigh = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup' AS study_name FROM Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup WHERE " + string + " ORDER BY p_value ASC"
    views.append(ukbbHigh)
    # ukbbLipid = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKBB_lipid_trait_Neal' AS study_name FROM Lipid_UKBB_lipid_trait_Neale WHERE " + string + " ORDER BY p_value ASC"
    # views.append(ukbbLipid)
    # ukbbStatin = "SELECT rsid, chr, bp, p_value, beta, trait, 'Lipid_UKB_statin_usage_Neale' AS study_name FROM Lipid_UKB_statin_usage_Neale WHERE " + string + " ORDER BY p_value ASC"
    # views.append(ukbbStatin)
   
    #create DataFrame from sql connection strings
    viewLen = len(views)
    i = 0
    df = pd.DataFrame(columns=["rsid", "chr", "bp", "p_value", "beta", "trait", "study_name"])
    while i < viewLen:
        tempDF = pd.read_sql(views[i], conn)
        df = df.append(tempDF)
        i += 1 
        
    print(df)
    return(df)

def main():
    tableLoad(sqlSearchString(hg19Search(setMargin()), setPValue()))
if __name__ == "__main__":
    main()