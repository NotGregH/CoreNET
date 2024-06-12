#!/usr/bin/env python

"""
This is a tool for Querying network data
in the SQLite database associated with
the Networking tool for Virtual
Plant 2.0. http://coruzzilab.bio.nyu.edu/

Version Notes
- Can be run in the command line
- Gene lists that are uploaded must be singular for now.
- New interaction object set up
    - New Sif file creation based on the new dictionary type object
- some code clean up to make it easier to read
- Need to test out speed vs node with edge query

"""

__author__ = 'Gregory Hamilton'
__version__ = '1.2.0'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

###########
# Packages
###########

import sqlite3 as sql
from csv import writer
import argparse
import time


###########
# Functions
###########

#@profile
def pull_interactions_0(gene_list, database, edges=None):

    """

    This is the function for pulling the interactions with zero hops.

    :param gene_list: A List object of genes.
    :param database: The database to pull the interactions from.

    :return: A list object set up like a sif file.

    """

    con = sql.connect(database)
    con.row_factory = sql.Row
    con.text_factory = str
    cur = con.cursor()
    ints = {}

    for i in gene_list:

        ints[i.upper()] = {}

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

    return ints


##### requires further optimization
#@profile
def pull_interactions_1(gene_list, database, edges=None):

    """

    This is the function for pulling the interactions with one hop.

    :param gene_list: A List object of genes.
    :param database: The database to pull the interactions from.
    :return: A list object set up like a sif file.

    """

    con = sql.connect(database)
    con.text_factory = str
    con.row_factory = sql.Row
    cur = con.cursor()
    ints = {}

    for i in gene_list:

        ints[i.upper()] = {}

    for i in gene_list:

        cur.execute('SELECT e.edge_name, n2.node_name FROM '
                    'interactions i, nodes n1, nodes n2, edges e WHERE '
                    'n1.node_name = ? AND n1.node_id=i.node_1_id AND '
                    'i.node_2_id=n2.node_id AND e.edge_id=i.edge_id', (i,))

        for row in cur.fetchall():

            edge, node2 = row

            if edges:

                if edge in edges:

                    try:

                        if ints[i][node2] != None:

                            dictValues = ints[i].get(node2,[])
                            dictValues.append(edge)
                            ints[i][node2] = list(set(dictValues))

                    except KeyError:

                        ints[i][node2] = [edge]

            else:

                try:
                    if ints[i][node2] != None:

                        dictValues = ints[i].get(node2,[])
                        dictValues.append(edge)
                        ints[i][node2] = list(set(dictValues))

                except KeyError:

                    ints[i][node2] = [edge]

        cur.execute('SELECT n1.node_name, e.edge_name FROM '
                    'interactions i, nodes n1, nodes n2, edges e WHERE '
                    'n2.node_name = ? AND n2.node_id=i.node_2_id AND '
                    'i.node_1_id=n1.node_id AND e.edge_id=i.edge_id', (i,))

        for row in cur.fetchall():

            node1, edge = row

            if edges:

                if edge in edges:

                    try:

                        if  ints[node1] != None:

                            try:

                                if ints[node1][i] != None:

                                    dictValues = ints[node1].get(i,[])
                                    dictValues.append(edge)
                                    ints[node1][i] = list(set(dictValues))

                            except KeyError:

                                ints[node1][i] = [edge]

                    except KeyError:

                        ints[node1] = {i : [edge]}

            else:

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

                    ints[node1] = {i : [edge]}

    return ints


def create_sif(ints, sif_file_name):
    """
    Creates a sif file from the interaction objects.
    :param ints: The interactions pulled from the database.
    :param sif_file_name: The name you'd like for the sif file.f
    :return: N/A
    """

    with open(sif_file_name, 'wb') as f:

        write = writer(f, delimiter='\t')

        for i in ints:

            if len(ints[i]) > 0:

                for node2 in ints[i]:

                    for edge in ints[i][node2]:

                        write.writerow((i, edge, node2))

    return f.close()


#@profile
def main(genes, db, hops, outfile, edges=None):

    #Profile Code
    #start = time.time()

    gene_list = []
    edge_list = None

    with open(genes, 'rb', -1) as f:

        for line in f:

            line = line.strip()
            if ">" in line:

                continue

            if "\t" in line:

                line = line.split("\t")
                line = line[0]

            gene_list.append(line.upper())

    if edges:

        with open(edges, 'rb', -1) as f:

            edge_list = []

            for line in f:

                line = line.strip()

                edge_list.append(line)

    if hops == 0:

        ints = pull_interactions_0(gene_list, db, edge_list)

    if hops == 1:

        ints = pull_interactions_1(gene_list, db, edge_list)

    create_sif(ints, outfile)

    #Profile Code
    #end = time.time()
    #print genes
    #print end - start

    return

############
############

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--dbname', help='Name for the database to pull the data')
    parser.add_argument("-g", "--genes", help="The file containing the gene lists to search")
    parser.add_argument("-e", "--edges", help="The list of edges you'd like to keep", default=None)
    parser.add_argument("-j", "--jumps", type=int, help="The number of jumps, must be 0 or 1 .", default=0)
    parser.add_argument("-o", "--outfile", help="The output sif file.", default='TempOut.txt')
    args = parser.parse_args()
    main(args.genes, args.dbname, args.jumps, args.outfile, args.edges)

