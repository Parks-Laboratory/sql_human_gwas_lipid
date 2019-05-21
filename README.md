# sql_human_gwas_lipid
This python program is used for query of interest Human_GWAS lipid genes. Since the size of total tables in the SQL server are nearly 20GB or expected to be larger, this program is intended to find only necessary data with the efficient algorithm over dozen of multiple large size tables.
* You need to enter table you want to search, gene name, desired margin, and desired p-value.
* You can select multiple values separated by comma.
* Then the corresponding result will show up, and this is a txt file with all significant (defined by setup p-value) SNP associated with genes.
* The returned output is sorted by p-value.
* This result can be saved as csv and automatically opened.
