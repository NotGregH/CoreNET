#!/usr/bin/env python

"""

This is a tool for Querying expression data
with a gene list and performing a pearson
correlation and generating a new sif file
with the edges being corr:r and a fourth column for p-value .  The tool
then queries the the AgrisPred data for the interactions
that meet the pvalue threshold.  This tool is
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

                continue

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
    A pearson correlation tool created by Akshay Jain.
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

    return round(r, 8), prob

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


####################################
####################################

if __name__ == '__main__':
    print "Part of the Cross Species Network Tool"