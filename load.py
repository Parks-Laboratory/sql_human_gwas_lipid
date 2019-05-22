import pyodbc
import csv
import os
from functools import reduce
import argparse


# "main" method which calls the functions and reads the txt file
def main():
    choice = input("what name you want to save this file in sql server ?\n")
    tablename = choice
    server = 'PARKSLAB'
    create = False
    path = "E:\\cross_ref\\load"
    database = "Human_GWAS"

    print (database)
    # connect to database

    cnxn = pyodbc.connect('DRIVER={SQL Server}' + \
                        ';SERVER=' + server + \
                        ';DATABASE=' + database + \
                        ';Trusted_Connection= Yes')

    print('connected to the database: %s successfully!' % database)
    cursor = cnxn.cursor()

    # if create is specified in the command line
    # if create:
    query = "create table {!s} ".format(tablename) + "(" + \
            "beta varchar (100)," \
            "bp varchar (100)," \
            "chr varchar (100)," \
            "gene_name varchar (100)," \
            "p_value varchar (100),"\
            "rsid varchar (100)," \
            "table_name varchar (100)," \
            "trait varchar (100));"


    print(query)
    cursor.execute(query)
    cursor.commit()
    print("table %s successfully created in database %s" % (tablename, database))



    # check if path exists
    if not os.path.isdir(path):
        print('path does not exist')
        exit(1)

        # change directory
    os.chdir(path)

    fileNames = []

    # for each of the files in the dir ending with txt, add to the list
    for file in os.listdir(path):
        if file.endswith(".txt"):
            fileNames.append(file)

    for fileName in fileNames:

        print("Reading from: " + str(fileName))

        with open(fileName, 'r') as txtFile:

            txtReader = csv.reader(txtFile, delimiter='\t')

            fileFormat = next(txtReader)

            errorCounter = 0

            index = 0
            # Resort each row to resemble the database format
            for rows in txtReader:

                index = index + 1

                list = []
                wanted = [0,1, 2, 3, 4, 5, 6, 7]

                listCounter = 0

                for cols in rows[0:]:
                    if listCounter in wanted:
                        list.append(cols)
                    listCounter += 1

                try:
                    query = "insert into dbo.{!s}".format(tablename) + \
                            " values ( "

                    for j in range(0, (len(list)-1)):
                        query += ("{!r}".format(list[j]) + ", ")

                    query += "{!r}".format(list[len(list) - 1])
                    query += ");"

                    cursor.execute(query)
                    cursor.commit()

                # write errmsg if file I/O exception
                except Exception as eex:
                    errorCounter += 1

                    if errorCounter == 1:
                        f = open("GC_{!r}_err.txt".format(tablename), "w")
                        f.write(str(list) + "\n")
                        print("ERROR in index " + str(index) + "!")
                    else:
                        print("ERROR in index " + str(index) + "!")
                        f.write(str(list) + "\n")

                else:
                    print("Insert " + str(index) + " was successful!")

    cursor.commit()
    print("File Read Done!" + str(fileName))
    cnxn.close()

main()