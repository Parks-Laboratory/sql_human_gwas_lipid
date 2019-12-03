# Analysis_human_gwas_lipid_From_SQL_server
This python script is ran from the command line taking an input txt file, a base pair margin, and a pvalue cutoff. The margin and pvalue are optional arguments, so if not entered they will defer to default values. The script will take these arguments and locate all of the genes in the hg19 table. Once we get our information from hg19, it is cross referenced with all of the lipid gwas tables. The result is a pandas data frame, sorted by ascending pvalue, that is saved as a text file in an output folder. 

# INPUT ARGUMENTS 
In order to run this script naviagate to the directory containing script. From here, enter the following commands:
"py LIPID_GWAS_SCRIPT.py [-h] [-m] [-p] -i" Brackets indicate optional arguments and should not be included when running the script from the command line. 

 - LEGEND:
 
 


