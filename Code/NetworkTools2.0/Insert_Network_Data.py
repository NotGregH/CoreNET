#!/usr/bin/env python

"""
This is a tool for inserting network data
into the SQLite database associated with
the Networking tool for Virtual
Plant 2.0. http://coruzzilab.bio.nyu.edu/

Version Notes 1.2.0
 - Improved the scripts speed so it takes less than two minutes for full insert.
 - removed the reverse interaction table
 - add indexing after data insert for speed
 - Code cleanup
    - Removed debug code, fixed object naming for consistency
"""

__author__ = 'Gregory Hamilton'
__version__ = '1.2.0'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

###########
# Packages
###########

import sqlite3 as sql
from sqlite3 import IntegrityError
import argparse
import fileinput
import time

###########
# Functions
###########


def get_nodes(file_list):

    """
    Compiles a single list of nodes from the input
    sif files.

    :param file_list: The file in list form from read_large_file

    :return nodesList:A list of nodes.

    """
    nodesList = []

    for i in file_list:

        nodesList.append((i[0]))
        nodesList.append((i[2]))

    # filters out redundant nodes
    nodesList = list(set(nodesList))

    return nodesList


def get_edges(file_list):

    """
    Pulls the edges from the input sif files

    :param file_list: The sif file in list form from read_large_file

    :return edgesList: A list of edges

    """

    edgesList = []

    for i in file_list:

        if i[1] not in edgesList:

            edgesList.append((i[1]))

    return edgesList


def insert_nodes(nodes_list, database):

    """

    Inserts the nodes into the database.

    :param nodes_list: A list of nodes

    :param database: The database to insert the nodes.

    :return: VOID

    """

    con = sql.connect(database)
    cur = con.cursor()
    count = 0
    for i in nodes_list:

        try:

            cur.execute('INSERT INTO nodes VALUES (NULL, ?)',
                        (i,))
            count = count +1

        ### allows for ignoring redundant inserts
        except IntegrityError:

            continue

    con.commit()
    con.close()
    print "Nodes Added : " , count
    return


def insert_edges(edges_list, database):

    """

    Inserts the edges into the database

    :param edges_list: List of edges
    :param database: The database to insert the edges to.

    :return: VOID

    """

    con = sql.connect(database)
    cur = con.cursor()
    count = 0
    for i in edges_list:

        try:

            cur.execute('INSERT INTO  edges VALUES (NULL, ?)', (i,))
            count = count +1
        ### allows for ignoring redundant inserts
        except IntegrityError:

            continue

    con.commit()
    con.close()
    print "Edges Added :" , count
    return


def insert_ints(fileList, db):
    """

    Inserts the interactions into the database.
    This is not the fastest but it works.

    :param fileList: The fileList from the sif file. Each item
    is a line from the sif file.
    :param db: The database you'd like to insert the data to.

    :return: VOID
    """

    con = sql.connect(db)
    cur = con.cursor()
    count = 0
    for i in fileList:

        try:

            cur.execute('INSERT INTO interactions (id, node_1_id, node_2_id,'
                        'edge_id) VALUES (NULL, (SELECT node_id FROM nodes '
                        'WHERE node_name = ?), (SELECT node_id FROM nodes '
                        'WHERE node_name = ?), (SELECT edge_id FROM edges '
                        'WHERE edge_name = ?) )', (i[0], i[2], i[1] ))
            count = count +1
        ### allows for ignoring redundant inserts
        except IntegrityError:

            continue

    con.commit()
    con.close()
    print "Interactions Added : " , count
    return

#@profile
def main(input_files, database):
    ## Profiling code
    #start = time.time()
    ##
    con = sql.connect(database)
    cur = con.cursor()

    cur.execute('DROP INDEX IF EXISTS nodeids_x')
    cur.execute('DROP INDEX IF EXISTS nodeids2_x')
    cur.execute('DROP INDEX IF EXISTS edges_x')
    cur.execute('DROP INDEX IF EXISTS nodes_x')

    for f in input_files:

        file_list = []
        print f
        for line in fileinput.input(f):

            line = line.strip('\n').strip('\r').split("\t")
            line[0] = line[0].upper()
            line[2] = line[2].upper()
            line = tuple(line)
            file_list.append(line)

        nodes = get_nodes(file_list)
        edges = get_edges(file_list)
        insert_edges(edges, database)
        insert_nodes(nodes, database)
        insert_ints(file_list, database)

    ## Creating Indexes

    cur.execute('CREATE INDEX nodeids_x ON interactions (node_1_id)')
    cur.execute('CREATE INDEX nodeids2_x ON interactions (node_2_id)')
    cur.execute('CREATE INDEX edges_x ON edges (edge_name, edge_id)')
    cur.execute('CREATE INDEX nodes_x ON nodes (node_name, node_id)')
    con.commit()
    con.close()
    ## Profiling code
    #end = time.time()
    #print "Finished", end - start
    return


###########
###########

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dbname", help='Name for the database to insert the data', required=True)
    parser.add_argument("-i", "--inputfiles", nargs="+", help="List of the sif files to insert", required=True)
    args = parser.parse_args()
    main(args.inputfiles, args.dbname)