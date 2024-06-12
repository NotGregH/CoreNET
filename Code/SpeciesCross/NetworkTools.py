#!/usr/bin/env python

"""
This tool contains functions
 for the networking
 part of the cross species
 network analysis package.

"""

__author__ = 'Gregory Hamilton'
__version__ = '0.0.5'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

############################
# Modules
############################

from Network_Correlation import format_expression_data,run_corr
from Network_Query import pull_interactions_0
from csv import writer

############################
# Functions
############################


def runCorrelation( expFile , pVal=None ):

    """

    Runs the correlation analysis for an expression data set.

    :param expFile1: Normalized Expression set from species one.

    :return corrSif: a correlation sif in list format.
    (i.e. { Gene1 : { Gene2 : [ ( corr:r::pval:p ) ] },..., GeneN : []}})
    """

    datadict, gene_arr, line_count, column_count = format_expression_data(expFile, None)
    corrSifDict = run_corr(datadict, gene_arr, line_count, column_count,pVal)

    return corrSifDict


def pullGeneList( expFile ):

    """

    Pulls a list of genes from the expression data set.

    :param expFile: The expression file.

    :return geneDict: A simple dictionary of the genes from the exp file.
        {gene1 : {}, gene2 : {}, .... , geneN : {}}

    """

    geneDict = {}

    with open(expFile,'rb',-1) as f:


        for line in f:

            if line.startswith("\t") or line.startswith(" "):

                continue

            else:

                line = line.strip().split('\t')[0]

                if line != "":

                    geneDict[line.upper()] = {}

    return geneDict


def filterSifFile( geneDict , sifFile ):

    """

    This filters the sif file for "0 hop" interactions based
     on the geneDict.  The interactions that will be kept must
     have genes from the geneDict in both node locations of the
     interaction.

    :param sifDict:  The geneDict created from pullGeneList().
    :param sifFile: The path to the sif file that will
    be filtered.

    :return sifDict: sifDict with the interactions added in.  The
    object is setup like a JSON object.
    (i.e. {Gene1 : {Gene2 : [edge1,..,edgeN]},..., GeneN : [edgeN]}})

    """
    sifDict = {}

    for i in geneDict:

        sifDict[i] = {}

    with open(sifFile,"rb",-1) as f:

        for line in f:

            line = line.strip().split("\t")
            line[0] = line[0].upper()
            line[2] = line[2].upper()

            if "met:" in line[1]:

                try:

                    if sifDict[line[0]] != None:

                        try:

                            if sifDict[line[0]][line[2]] != None:

                                dictValues = sifDict[line[0]].get(line[2], [])
                                dictValues.append( line[1] )
                                sifDict[line[0]][line[2]] = dictValues

                        except KeyError:

                            sifDict[line[0]][line[2]] = [line[1]]

                except KeyError:

                    try:

                        if sifDict[line[2]] != None:

                            try:

                                if sifDict[line[0]] != None :

                                    dictValues = sifDict[line[0]].get(line[2], [])
                                    dictValues.append( line[1] )
                                    sifDict[line[0]][line[2]] = dictValues

                            except KeyError:

                                sifDict[line[0]] = {line[2] : [line[1]]}

                    except KeyError:

                        pass

            else:

                try:

                    if sifDict[line[0]] != None:

                        try:

                            if sifDict[line[2]] != None:

                                try:

                                    if sifDict[line[0]][line[2]] != None:

                                        dictValues = sifDict[line[0]].get(line[2], [])
                                        dictValues.append( line[1] )
                                        sifDict[line[0]][line[2]] = dictValues

                                except KeyError:

                                    sifDict[line[0]][line[2]] = [line[1]]

                        except KeyError:

                            continue

                except KeyError:

                    continue

    ##### Removing empty dictionary keys
    removeList = []

    for i in sifDict:

        if len(sifDict[i]) == 0 :

            removeList.append(i)

    for i in removeList:

        del sifDict[i]

    return sifDict


def knownInteractions( geneList ,  interactionSif=None ,
                       interactionDB=None , edgeList= None ):

    """

    Pulls the known interactions using the geneList generated from
    the expression file.

    :param geneList: List of genes to query. Generated by pullGeneList().
    :param interactionDB: A database of known interactions created using the virtual plant
    networking tool pipeline.
    :param interactionSif: A known network sif file.

    :return speciesSif: A dictionary type "sif" object of the "0 hop" interactions
    pulled from the database.

    """
    speciesSif = {}

    if interactionDB != None and interactionSif == None:

        if edgeList:

            speciesSif = pull_interactions_0(geneList,interactionDB,edgeList)

        else:

            speciesSif = pull_interactions_0(geneList,interactionDB,None)

    if interactionSif != None and interactionDB == None:

        if len(interactionSif) == 1:

            speciesSif = filterSifFile(geneList,interactionSif[0])

        elif len(interactionSif) > 1:

            gate = 0

            for fileName in interactionSif:

                if gate == 0:

                    speciesSif = filterSifFile(geneList,fileName)
                    gate = gate + 1

                if gate >= 1 :

                    tempSif = filterSifFile(geneList,fileName)
                    speciesSif = sifUnion(tempSif,speciesSif)

    if interactionDB != None and interactionSif != None:
        if edgeList:

            speciesSif = pull_interactions_0(geneList,interactionDB,edgeList)

        else:

            speciesSif = pull_interactions_0(geneList,interactionDB,None)

        if len(interactionSif) == 1:

            tempSif = filterSifFile(geneList,interactionSif[0])
            speciesSif = sifUnion(tempSif,speciesSif)

        elif len(interactionSif) > 1:

            gate = 0

            for fileName in interactionSif:

                if gate == 0:

                    tempSif = filterSifFile(geneList,fileName)
                    speciesSif = sifUnion(tempSif,speciesSif)
                    gate = gate + 1

                if gate >= 1 :

                    tempSif = filterSifFile(geneList,fileName)
                    speciesSif = sifUnion(tempSif,speciesSif)

    return speciesSif


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
    if ".SIF" not in outfile:
        
        outfile = outfile +".SIF"

    with open( outfile , "wb" , -1 ) as f:

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

### Probably can delete ## is not used
def sifIntersect( corrSif , knownSif ):

    """

    Intersects the corrSif and knownSif object.

    :param corrSif:  A correlation Sif object from runCorrelation().
    :param knownSif: A knownSif object generated from knownInteractions().

    :return intersectSif: A dictionary sif object set up like a JSON object.
    ( { Gene1 : { Gene2 : [ edge1, corr:r::pval: p , edge3 ] } } )

    """

    intersectSif = {}

    for i in knownSif:

        intersectSif[i] = {}


    for i in knownSif:

        try:

            if corrSif[i] != None:

                for l in knownSif[i]:

                    try:

                        if corrSif[i][l] != None:

                            dictValues = intersectSif[i].get( l , [] )

                            for k in knownSif[i][l]:

                                dictValues.append(k)

                            dictValues.append(corrSif[i][l])
                            intersectSif[i][l] = dictValues

                    except KeyError:

                        ### To make sure we keep any self feedback interactions.
                        if l == i:

                            intersectSif[i][l] = knownSif[i][l]

                        else:

                            continue

        except KeyError:

            continue

    removeList = []

    for i in intersectSif:

        if len(intersectSif[i]) == 0 :

            removeList.append(i)

    for i in removeList:

        del intersectSif[i]

    return intersectSif

#############################
#############################

if __name__ == '__main__':
    print "Not for commandline use"