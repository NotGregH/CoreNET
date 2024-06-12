#!/usr/bin/env python

"""

This is a tool for Querying expression data
with a gene list and performing a pearson
correlation and generating a new sif file
with the edges being corr:r and a fourth column for p-value .  The tool
then queries the database for the interactions
that mee the p-Value threshold.  This tool is
part of the Networking tool for Virtual
Plant 2.0. http://coruzzilab.bio.nyu.edu/

Version Notes:
    - Introduced sql row factory
    - Added geneList as a default none param for format_expression_data
        - This allows the function to be run on the full expression set if desired
    - Added edgeList as a default none param to pull_interactions()
        - This allows for selective filtering of edgeTypes pulled from the database
    - Added pval as a default none param to run_corr()
        - This allows the user to decide if they would like to filter based on a pvalue
    - Cleaned up the code for readability
        - Added white space and cleaned up some naming
    - Still need to work on "chunking"
    - Need to work on changing the sif object to a dictionary type

"""

__author__ = 'Gregory Hamilton'
__version__ = '1.1.0'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"


####################################
# Modules
####################################

import argparse as argp
import sqlite3 as sql
import numpy as np
import scipy.stats as spst
from csv import writer
import time

####################################
# Functions
####################################

def format_expression_data(expressionFile, gene_list=None):

    """

    Pulls the expression data from a larger file for
     only genes in the gene_list and puts it in the format required
     to run pearson correlation.

    :param gene_list: List of genes to query
    :param expression_file: The file containing normalized
    expression data.

    :return: datadict  = The expression data in a dictionary object ( i.e.
        { gene1 : [exp array 1], ...,geneN : [Exp Array N]} )
        gene_arr = A Numpy array with every gene in the data_dict.
        line_count = The number genes being analyzed.
    """


    expressionData = []
    datadict = dict()
    line_count = 0

    with open(expressionFile, "rb", -1) as f:

        for line in f:

            if line.startswith("\t") or line.startswith(" "):

                pass

            else:


                line = line.strip().split("\t")

                if gene_list:

                    if line[0] in gene_list:

                        expressionData.append(line)
                        line_count += 1

                else:

                    expressionData.append(line)
                    line_count += 1

    column_count = len(expressionData[0]) - 1

    for i, row in enumerate(expressionData):

        for j, item in enumerate(row):

            expressionData[i][j] = item

    for i in range(0,len(expressionData)):

        datadict[expressionData[i][0]] = np.array(expressionData[i][1:])

    gene_arr = np.array([expressionData[i][0] for i in range(0, len(expressionData))])

    return datadict, gene_arr, line_count, column_count


def pearson_corr(x, y, xm, ym, ssxm, ssym, n):

    """

    Each param are the standard params required to preform correlation.
    :param x:
    :param y:
    :param xm:
    :param ym:
    :param ssxm:
    :param ssym:
    :param n:
    :return:
    """

    r_num = np.add.reduce(xm * ym)
    r_den = np.sqrt(ssxm * ssym)
    r = r_num / r_den
    r = max(min(r, 1.0), -1.0)
    df = n-2

    if abs(r) == 1.0:

        prob = 0.0

    else:

        t_squared = r*r * (df / ((1.0 - r) * (1.0 + r)))
        prob = spst.betai(0.5*df, 0.5, df / (df + t_squared))

    return round(r, 6), prob

#@profile
def run_corr(datadict, gene_arr, line_count, column_count, pval=None):

    """

    Runs pearson correlation on expression data
    from the format expression data function.

    :param datadict: datadict  = The expression data in a dictionary object ( i.e.
                     { gene1 : [exp array 1], ...,geneN : [Exp Array N]} )
    :param gene_arr: A Numpy array with every gene in the data_dict.
    :param line_count:  Int Object representing the number genes being analyzed.
    :param column_count: Int Object representing the number experiments being analyzed.

    :return: A sif object with Gene1 corr:p_value Gene2 r:Correlation_coefficient

    """

    # initializing all required objects.
    solution_arr = np.zeros((line_count, line_count, 2))
    m = np.zeros((line_count,column_count), dtype=np.float64)
    ssm = np.zeros((line_count,), dtype=np.float64)
    # filling the pearson correlation required objects

    for i in range(0, line_count):

        x = np.asfarray(datadict[gene_arr[i]])
        n = len(x)
        mx = x.mean()
        m[i] = x-mx
        ssm[i] = spst.ss(m[i])

    corr_sif = {}



    for i in range(0, line_count):

        for j in range(i+1, line_count):

            solution_arr[i, j] = pearson_corr(datadict[gene_arr[i]], datadict[gene_arr[j]],
                                           m[i], m[j], ssm[i], ssm[j], n)

            if pval:


                if (solution_arr[i, j, 1] < pval) & (gene_arr[i] != gene_arr[j]):
                    node1 = gene_arr[i].replace("\"", "")
                    node2 = gene_arr[j].replace("\"", "")
                    p = "pval:" + str(solution_arr[i, j, 1])
                    r = "corr:" + str(solution_arr[i, j, 0])
                    try:

                        corr_sif[node1][node2] = [r +"::"+ p ]

                    except KeyError:

                        corr_sif[node1] = {node2 : [r + "::" + p]}

            else:

                if (gene_arr[i] != gene_arr[j]):

                    node1 = gene_arr[i].replace("\"", "")
                    node2 = gene_arr[j].replace("\"", "")
                    p = "pval:" + str(solution_arr[i, j, 1])
                    r = "corr:" + str(solution_arr[i, j, 0])

                    try:

                        corr_sif[node1][node2] = [r +"::"+ p ]

                    except KeyError:

                        corr_sif[node1] = {node2 : [r + "::" + p]}

    return corr_sif


def correlation(expressionFile, gene_file=None, pval=None):

    if gene_file != None:

        gene_list = []

        with open(gene_file, 'rb', -1) as f:

            for line in f:

                gene_list.append(line.strip())

            f.close()

    else:

        gene_list = None

    datadict, gene_arr, line_count, column_count = format_expression_data(expressionFile, gene_list)

    corrSif = run_corr(datadict,gene_arr,line_count,column_count,pval)

    return corrSif


def sifUnion( corrSif , speciesSif ):

    """

    Creates a union of the known interactions Sif and the
    correlation sif object.

    :param corrSif: A dictionary type correlation sif object.
    :param speciesSif: A dictionary type known interaction sif object.

    :return unionSif:  A dictionary type sif object from the union of the two
    dictionary object put in.

    """

    unionSif = {}

    if len(corrSif) > len(speciesSif):

        unionSif = corrSif

        for i in speciesSif:

            try:

                if unionSif[i] != None:

                    for l in speciesSif[i]:

                        try:

                            if unionSif[i][l] != None:

                                dictValues = unionSif[i].get( l , [] )

                                for k in speciesSif[i][l]:

                                    dictValues.append(k)

                                unionSif[i][l] = list(set(dictValues))

                        except KeyError:

                            unionSif[i][l] = speciesSif[i][l]

            except KeyError:

                unionSif[i] = speciesSif[i]

    else:

        unionSif = speciesSif

        for i in corrSif:

            try:

                if unionSif[i] != None:

                    for l in corrSif[i]:

                        try:

                            if unionSif[i][l] != None:

                                dictValues = unionSif[i].get( l , [] )

                                for k in corrSif[i][l]:

                                    dictValues.append(k)

                                unionSif[i][l] = list(set(dictValues))

                        except KeyError:

                            unionSif[i][l] = corrSif[i][l]

            except KeyError:

                unionSif[i] = corrSif[i]

    return unionSif


def pull_interactions(corr_sif, database, edgeFile=None):

    """

    This is the function for pulling the interactions with zero hops.

    :param gene_list: A List object of genes.
    :param database: The database to pull the interactions from.

    :return: A list object set up like a sif file.

    """
    if edgeFile:

        with open(edgeFile, 'rb', -1) as f:

            edges = []

            for line in f:

                line = line.strip()

                edges.append(line)
    else:

        edges = None

    con = sql.connect(database)
    con.row_factory = sql.Row
    con.text_factory = str
    cur = con.cursor()
    ints = {}
    gene_list = []

    for i in corr_sif:

        ints[i.upper()] = {}
        gene_list.append(i.upper())

        for l in corr_sif[i]:

            ints[l.upper()] = {}

            gene_list.append(l.upper())

    gene_list = list(set(gene_list))
    

    for i in gene_list:

        cur.execute('SELECT e.edge_name, n2.node_name FROM '
                    'interactions i, nodes n1, nodes n2, edges e WHERE '
                    'n1.node_name = ? AND n1.node_id=i.node_1_id AND '
                    'i.node_2_id=n2.node_id AND e.edge_id=i.edge_id', (i,))


        for row in cur.fetchall():

            edge, node2 = row

            if edges:

                if edge in edges:

                    if "met:" not in edge and "miRNA:" not in edge:

                        try:

                            if ints[node2] != None:

                                try:

                                    if ints[i][node2] != None:

                                        dictValues = ints[i].get(node2,[])
                                        dictValues.append(edge)
                                        ints[i][node2] = list(set(dictValues))

                                except KeyError:

                                    ints[i][node2] = [edge]

                        except KeyError:

                            continue

                    if "met:" in edge or "miRNA:" in edge:

                        try:

                            if ints[i][node2] != None:

                                dictValues = ints[i].get(node2,[])
                                dictValues.append(edge)
                                ints[i][node2] = list(set(dictValues))

                        except KeyError:

                            ints[i][node2] = [edge]

            else:

                if "met:" not in edge and "miRNA:" not in edge:

                    try:

                        if ints[node2] != None:

                            try:

                                if ints[i][node2] != None:

                                    dictValues = ints[i].get(node2,[])
                                    dictValues.append(edge)
                                    ints[i][node2] = list(set(dictValues))

                            except KeyError:

                                ints[i][node2] = [edge]

                    except KeyError:

                        continue

                if "met:" in edge or "miRNA:" in edge:

                    try:

                        if ints[i][node2] != None:

                            dictValues = ints[i].get(node2,[])
                            dictValues.append(edge)
                            ints[i][node2] = list(set(dictValues))

                    except KeyError:

                        ints[i][node2] =  [edge]


        cur.execute('SELECT n1.node_name, e.edge_name FROM '
                    'interactions i, nodes n1, nodes n2, edges e WHERE '
                    'n2.node_name = ? AND n1.node_id=i.node_1_id AND '
                    'i.node_2_id=n2.node_id AND e.edge_id=i.edge_id', (i,))

        for row in cur.fetchall():

            node1, edge  = row

            if edges:

                if edge in edges:

                    if "met:" not in edge and "miRNA:" not in edge:

                        try:

                            if ints[node1] != None:

                                try:

                                    if ints[node1][i] != None:

                                        dictValues = ints[node1].get(i,[])
                                        dictValues.append(edge)
                                        ints[node1][i] = list(set(dictValues))

                                except KeyError:

                                    ints[node1][i] = [edge]

                        except KeyError:

                            continue

                    if "met:" in edge or "miRNA:" in edge:

                        try:

                            if ints[node1] != None:

                                try:

                                    if ints[node1][i] != None:

                                        dictValues = ints[node1].get(i,[])
                                        dictValues.append(edge)
                                        ints[node1][i] = list(set(dictValues))

                                except KeyError:

                                    ints[node1][i] = [edge]

                        except KeyError:

                            ints[node1] =  { i : [edge] }

            else:

                if "met:" not in edge and "miRNA:" not in edge:

                    try:

                        if ints[node1] != None:

                            try:

                                if ints[node1][i] != None:

                                    dictValues = ints[node1].get(i,[])
                                    dictValues.append(edge)
                                    ints[node1][i] = list(set(dictValues))

                            except KeyError:

                                ints[node1][i] = [edge]

                    except KeyError:

                        continue

                if "met:" in edge or "miRNA:" in edge:

                    try:

                        if ints[node1] != None:

                            try:

                                if ints[node1][i] !=  None:

                                    dictValues = ints[node1].get(i,[])
                                    dictValues.append(edge)
                                    ints[node1][i] = list(set(dictValues))

                            except KeyError:

                                ints[node1][i] = [ edge ]

                    except KeyError:

                        ints[node1] =  { i : [ edge ] }

    corr_sif = sifUnion(corr_sif,ints)

    return corr_sif



def saveSif( crossSif , outfile ):
    """
    This function takes a dictionary sif object and
    creates a sif files and a .tsv file.  The
    tab separated file contains the raw data.

    :param crossSif: The
    :param outfile: The name you would like
    the output sif file to have.

    :return Null
    """


    with open( outfile+".SIF" , "wb" , -1 ) as f:

        write = writer(f, delimiter='\t')
        write.writerow(["node_1","edge","node_2"])

        for i in crossSif:

            if len(crossSif[i]) > 0:

                for node2 in crossSif[i]:

                    for edge in crossSif[i][node2]:

                        if "::" in edge:

                            edge = edge.split("::")

                            if "corr:" in edge[0] and "corr:-" not in edge[0]:

                                write.writerow( ( i, "corr", node2, edge[0], edge[1]))

                            elif "corr:-" in edge[0]:

                                write.writerow( ( i, "negCorr", node2, edge[0], edge[1] ) )

                        else:

                            write.writerow( ( i, edge, node2 ) )
        f.close()

    return

#@profile
def main(gene_file, expression_file, database,  outfile, pval=None,edges=None):

    ## Profile Code
    #print gene_file
    #start = time.time()

    corr_sif = correlation(expression_file,gene_file,pval)

    ## Profile Code
    #corrEnd = time.time()
    #print "Correlation Over" , corrEnd - start
    #qStart = time.time()
    
    corr_sif = pull_interactions(corr_sif, database, edges)
    
    ## Profile Code
    #qEnd = time.time()
    #print "Query Over" , qEnd - qStart

    saveSif(corr_sif, outfile)

    ## Profile
    #end = time.time()
    #print "Finished", end - start, "\n\n"

    return

####################################
####################################

if __name__ == '__main__':
    parser = argp.ArgumentParser()
    parser.add_argument("-d", '--dbname', help='Name for the database to insert the data')
    parser.add_argument("-g", "--genes", help="The List of genes you would like to query against"
                                              " the expression data.", default=None)
    parser.add_argument("-e", "--expression_data", help="The expression data you'd like to query from")
    parser.add_argument("-l", "--edgeList", help="The edgeList you'd like to query for", default=None)
    parser.add_argument("-p", "--pval",help="The pvalue cutoff you'd like to use", default=None)
    parser.add_argument("-o", "--outfile", help="The output sif file.", default="tempCorr.sif")
    args = parser.parse_args()
    main(args.genes, args.expression_data, args.dbname,args.outfile,args.pval,args.edgeList)
