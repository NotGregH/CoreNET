#!/usr/bin/env python

"""
This contains functions for
the orthology part of the
cross species network
analysis package.

"""

__author__ = 'Gregory Hamilton'
__version__ = '0.0.5'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

############################
# Modules
############################

from csv import writer
import time

############################
# Functions
############################
def pullOrthologyData( orthologyDataFile , geneList1 ,
                       geneList2 ):

    """

    This function pulls data from an orthology data file in the OrthoMCL format.
    The gene lists are generated by pullGeneList() and geneList1 is the shorter list and geneList2
    is the longer list.  They must be from two separate species.

    :param orthologyDataFile: The Orthology datafile in the orthoMCL format
    :param geneList1: The longer geneList .
    :param geneList2: The shorter geneList.

    :return orthoDict: {Species 1 gene : { species 2 gene : orthoId , ... , species 2 gene n : orthoID } }
            revOrthoDict:  The same as orthoDict but with species 2 genes as the primary keys.

    """

    orthoFilter2 = []
    orthoFilter = []
    orthoDict = {}
    revOrthoDict = {}

    #start = time.time()
    with open( orthologyDataFile , "rb" , -1 ) as f:

        for line in f:

            line = line.strip()

            for gene_sp1 in geneList1:

                if gene_sp1 in line or gene_sp1 in line.upper() :

                    orthoFilter.append( line )

                    break

                else:

                    continue

        f.close()
    #end = time.time()
    #print "Filter 1" , end-start

    for line in orthoFilter:

        for gene_sp2 in geneList2:

            if str( line ).find( gene_sp2 ) != -1 or str( line ).upper().find( gene_sp2 ) != -1 :

                orthoFilter2.append( line )

                break

            else:

                continue

    for line in orthoFilter2:

        for gene_sp1 in geneList1:

            if line.find( gene_sp1 ) != -1 or line.upper().find( gene_sp1 ) != -1:

                for gene_sp2 in geneList2:

                    if line.find( gene_sp2 ) != -1 or line.upper().find( gene_sp2 ) != -1 :

                        try:

                            temp = orthoDict.get( gene_sp1 , [] )
                            temp.append(gene_sp2)
                            orthoDict[gene_sp1] = list(set(temp))

                        except KeyError:

                            orthoDict[gene_sp1] = [ gene_sp2 ]

    for line in orthoFilter2:

        for gene_sp2 in geneList2:

            if line.find(gene_sp2) != -1 or line.upper().find(gene_sp2) != -1:

                for gene_sp1 in geneList1:

                    if line.find(gene_sp1) != -1 or line.upper().find(gene_sp1) != -1:

                        try:

                            temp = revOrthoDict.get(gene_sp2, [])
                            temp.append(gene_sp1)
                            revOrthoDict[gene_sp2] = list(set(temp))

                        except KeyError:

                            revOrthoDict[gene_sp2] = [ gene_sp1 ]

    return orthoDict, revOrthoDict


def filterFlatOrthoFile(orthoFile, geneList1,
                            geneList2):

    """

    Filters a flat orthology file set up like a tab delimited sif file.
    (i.e. Gene1 orthoType Gene2)

    geneList1 will output to orthoDict.

    :param orthoFile: A orthology file structured like a sif file
    :param geneList1: Species 1 gene list
    :param geneList2: Species 2 gene list

    :return:

    """
    orthoDict = {}

    with open(orthoFile, "rb", -1) as f:
        for line in f:

            if "\t" in line:

                line = line.strip().split("\t")

            elif " " in line and "\t" not in line:

                line = line.strip().split(" ")

            gene_sp1 = line[0].upper()
            gene_sp2 = line[2].upper()

            try:

                if geneList1[gene_sp1] != None:

                    try:

                        if geneList2[gene_sp2] != None:

                            try:

                                temp = orthoDict.get( gene_sp1 , [])
                                temp.append(gene_sp2)
                                orthoDict[gene_sp1] = list(set(temp))

                            except KeyError:

                                orthoDict[ gene_sp1 ] =  [ gene_sp2 ]

                    except KeyError:

                        continue

            except KeyError:

                continue

    return orthoDict


def unionOrthoData( orthoDict1 , orthoDict2 ):

    """

    Takes two orthology dictionaries from the same species and merges them.

    :param orthoDict1:
    :param orthoDict2:

    :return orthoDict( 1 or 2):

    """

    if len(orthoDict1) >= len(orthoDict2):

        for i in orthoDict2:

            try:

                if orthoDict1[i] != None:

                    temp = orthoDict2[i] + orthoDict1[i]
                    orthoDict1[i] = list(set(temp))


            except KeyError:

                orthoDict1[i] = orthoDict2[i]

        return orthoDict1

    elif len(orthoDict1) < len(orthoDict2):

        for i in orthoDict1:

            try:

                if orthoDict2[i] != None:

                    temp = orthoDict2[i] + orthoDict1[i]
                    orthoDict2[i] = list(set(temp))


            except KeyError:

                orthoDict2[i] = orthoDict1[i]

        return orthoDict2

def saveOrtho(ortho, outfile):
        """
        This function takes a dictionary sif object and
        creates a sif files.

        :param crossSif: The
        :param outfile: The name you would like
        the output sif file to have.

        :return Null
        """

        with open(outfile, "wb", -1) as f:

            write = writer(f, delimiter='\t')

            for i in ortho:

                for node2 in ortho[i]:

                    write.writerow((i, "ortho", node2))

        return


#############################
#############################


if __name__ == '__main__':
    print "Not for commandline use"