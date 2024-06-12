#!/usr/bin/env python

"""

This tool is meant for performing cross species
network analysis.

Version Notes

- Adding metabolic edge functionality
- Updating orthoDict type to replace the type of orthology with an int
- Changing correlation type to combine pval and corr into a single edge item.

"""

__author__ = 'Gregory Hamilton'
__version__ = '0.0.5'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

############################
# Modules
############################

import OrthologyTools
import CrossTools
import NetworkTools
import argparse
import sys
import os.path
import time

############################
# Functions
############################


def fileAndDatabaseCheck(expFile1, expFile2,sifFiles1=None,sifFiles2=None,db1=None,
                         db2=None,orthoMCL=None,orthoFile1=None,orthoFile2=None):
    """

    This function is for checking the path location to all of the files
    required for the cross species tool.

    All of the arguments are the file based input arguments.

    """

    if os.path.isfile(expFile1) == False:

        print "ERROR: Species 1 Expression File does not exist."
        print "Exiting ... "
        sys.exit()

    if os.path.isfile(expFile2) == False:

        print "ERROR: Species 2 Expression File does not exist."
        print "Exiting ... "
        sys.exit()

    if sifFiles1 != None:

        for i in sifFiles1:

            if os.path.isfile(i) == False:

                print "ERROR: Known Network file %d not found." %(i)
                print "Exiting ... "
                sys.exit()

    if sifFiles2 != None:

        for i in sifFiles2:

            if os.path.isfile(i) == False:

                print "ERROR: Known Network file %d not found." %(i)
                print "Exiting ... "
                sys.exit()

    if db1 != None:

        if os.path.isfile(db1) == False:

            print "ERROR: Known Network database %d not found." %(db1)
            print "Exiting ... "
            sys.exit()

    if db2 != None:

        if os.path.isfile(db2) == False:

            print "ERROR: Known Network database %d not found." %(db2)
            print "Exiting ... "
            sys.exit()

    if orthoMCL != None:

        if os.path.isfile(orthoMCL) == False:

            print "ERROR: OrthoMCL file %d not found." %(orthoMCL)
            print "Exiting ... "
            sys.exit()

    if orthoFile1 != None:
        for i in orthoFile1:

            if os.path.isfile(i) == False:

                print "ERROR: OrthoMCL file %d not found." %(i)
                print "Exiting ... "
                sys.exit()

    if orthoFile2 != None:
        for i in orthoFile2:

            if os.path.isfile(i) == False:

                print "ERROR: OrthoMCL file %d not found." %(i)
                print "Exiting ... "
                sys.exit()

    return


def speciesNetwork( expFile , pval , db = None ,
                    sif = None , edgeList = None ):
    """

    :param expFile:
    :param pval:
    :param db:
    :param sif:
    :param edgeList:

    :return:
    """

    sifDict = {}
    genes = NetworkTools.pullGeneList( expFile )
    corrSif = NetworkTools.runCorrelation( expFile , pval )

    if db != None and sif != None :

        sifDict = NetworkTools.knownInteractions( genes , sif )
        tempDict = NetworkTools.knownInteractions( genes ,
                                      None , db , edgeList )
        sifDict = NetworkTools.sifUnion(sifDict,tempDict)

    elif db != None and sif == None:

        sifDict = NetworkTools.knownInteractions( genes , sif ,
                                     db , edgeList )

    elif db == None and sif != None:

        sifDict = NetworkTools.knownInteractions( genes , sif )

    sifDict = NetworkTools.sifUnion( corrSif , sifDict )

    return sifDict, genes


def orthoData( geneDict1 , geneDict2 ,
               orthoMCL = None ,
               orthoSif1 = None ,
               orthoSif2 = None ):
    """

    :param geneDict1:
    :param geneDict2:
    :param orthoMCL:
    :param orthoSif:

    :return:
    """
    spec1ortho = {}
    spec2ortho = {}

    if orthoMCL != None:



        spec1ortho , spec2ortho = OrthologyTools.pullOrthologyData( orthoMCL ,
                                                                    geneDict1 ,
                                                                    geneDict2 )


        if orthoSif1 != None:

            for i in orthoSif1:

                temp1ortho = OrthologyTools.filterFlatOrthoFile( i ,
                                                                 geneDict1 ,
                                                                 geneDict2 )

                spec1ortho = OrthologyTools.unionOrthoData( spec1ortho , temp1ortho )

        if orthoSif2 != None:

            for i in orthoSif2:

                temp2ortho = OrthologyTools.filterFlatOrthoFile( i ,
                                                                 geneDict2 ,
                                                                 geneDict1 )
                spec2ortho = OrthologyTools.unionOrthoData( spec2ortho , temp2ortho )

    if orthoMCL == None:

        if len(orthoSif1) == 1:

            spec1ortho = OrthologyTools.filterFlatOrthoFile(orthoSif1[0],
                                                                  geneDict1,
                                                                  geneDict2)
        if len(orthoSif1) > 1:

            spec1ortho, spec2ortho = OrthologyTools.filterFlatOrthoFile(orthoSif1[0],
                                                                            geneDict1,
                                                                            geneDict2)
            for i in orthoSif1[1:]:

                temp1ortho, temp2ortho = OrthologyTools.filterFlatOrthoFile(i,
                                                                                geneDict1,
                                                                                geneDict2)
                spec1ortho = OrthologyTools.unionOrthoData(spec1ortho, temp1ortho)

        if len(orthoSif2) == 1:

            spec2ortho= OrthologyTools.filterFlatOrthoFile(orthoSif2[0],
                                                           geneDict2,
                                                           geneDict1)

        if len(orthoSif2) > 1 :

            spec2ortho = OrthologyTools.filterFlatOrthoFile(orthoSif2[0],
                                                            geneDict2,
                                                            geneDict1)

            for i in orthoSif2[1:]:

                temp2ortho, = OrthologyTools.filterFlatOrthoFile(i,
                                                                 geneDict2,
                                                                 geneDict1)


                spec2ortho = OrthologyTools.unionOrthoData(spec2ortho, temp2ortho)

    return spec1ortho, spec2ortho


def crossSpeciesNetworks(sifDict1 , sifDict2, ortho1, ortho2,
                         direction , outfile ):
    """

    :param sifDict1:
    :param sifDict2:
    :param ortho1:
    :param ortho2:
    :param direction:

    :return:
    """
    # Debug Code
    #NetworkTools.saveSif(sifDict1, "Ara_PreFilter.sif")
    #NetworkTools.saveSif(sifDict2, "Rice_PreFilter2.sif")
    ##

    sifDict1 = CrossTools.filterSif( sifDict1 , ortho1 )
    sifDict2 = CrossTools.filterSif( sifDict2 , ortho2 )

    #OrthologyTools.saveOrtho(ortho1 , "AraOrtho.sif")
    #OrthologyTools.saveOrtho(ortho2, "RiceOrtho.sif")

    # Debug Code
    #NetworkTools.saveSif(sifDict1, "Ara_Filter.sif")
    #NetworkTools.saveSif(sifDict2, "Rice_Filter.sif")
    ##

    crossSif1to2 = None
    crossSif2to1 = None

    if direction == "B":

            crossSif1to2 = CrossTools.strictXIntersectV2( sifDict1 , sifDict2 , ortho1 )
            crossSif2to1 = CrossTools.strictXIntersectV2( sifDict2 , sifDict1 , ortho2 )

    if direction == "1to2":

            crossSif1to2 = CrossTools.strictXIntersectV2( sifDict1 , sifDict2 , ortho1 )

    if direction == "2to1":

            crossSif2to1 = CrossTools.strictXIntersectV2( sifDict2 , sifDict1 , ortho2 )

    if len( outfile ) == 1:

        if crossSif1to2 != None:

            NetworkTools.saveSif( crossSif1to2 , outfile[0] )

        if crossSif2to1 != None:

            NetworkTools.saveSif( crossSif2to1 , outfile[0] )

    if len( outfile ) == 2:

        if crossSif1to2 != None:
            #print "file 1"

            NetworkTools.saveSif( crossSif1to2 , outfile[1] )

        if crossSif2to1 != None:

            #print "file 2"
            NetworkTools.saveSif( crossSif2to1 , outfile[0] )

    return


############################
############################


#@profile
def main( expFile1 , expFile2 ,  pVal , outfile , direct ,
          orthoMCL = None, ortho1=None , ortho2=None , db1=None , db2=None ,
          sif1=None , sif2=None, edgeList=None ):
    
    # Benchmark Code
    #print expFile1 , "\n", expFile2 , "\n"
    #start = time.time()

    fileAndDatabaseCheck(expFile1, expFile2, sif1 , sif2 , db1 ,
                         db2 , orthoMCL, ortho1 , ortho2 )

    # Benchmark Code
    #netStart = time.time()

    species1Sif, species1genes = speciesNetwork( expFile1 , pVal ,
                                                 db1 , sif1 , edgeList )
    # Benchmark Code
    #net2start = time.time()
    #print "Sp1", net2start - netStart , "\n"

    species2Sif, species2genes = speciesNetwork( expFile2 , pVal ,
                                                 db2 , sif2 , edgeList )

    # Benchmark Code
    #netEnd = time.time()
    #print "Sp2" , netEnd - net2start  , "\n"
    #print "Network Total Time" , netEnd - netStart , "\n"
    #orthoStart = time.time()

    species1Ortho, species2Ortho = orthoData( species1genes , species2genes ,
                                              orthoMCL , ortho1 , ortho2 )

    # Benchmark Code
    #orthoEnd = time.time()
    #print "Ortho Time", orthoEnd - orthoStart , "\n"
    #crossStart = time.time()

    crossSpeciesNetworks( species1Sif , species2Sif , species1Ortho ,
                          species2Ortho, direct , outfile )

    # Benchmark Code
    #crossEnd = time.time()
    #print "Cross Time", crossEnd - crossStart , "\n"
    #print"Total Time", crossEnd - start, "\n"  , "\n" , "\n"

    return

#############################
#############################


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    ##### Species 1 Arguments ######
    parser.add_argument("-e1","--expSpecies1",
                        help="The path to the species 1 tab delimited normalized "
                             "read count expression file.")
    parser.add_argument("-d1","--dbSpecies1",
                        help="The path to the species 1 interaction database.",
                        default=None)
    parser.add_argument("-s1","--sifSpecies1",
                        help="THe path(s) + name(s) to/of the species 1"
                             " sif file.", default=None, nargs="+")
    parser.add_argument("-o1","--orthoSpecies1",
                        help="THe path(s) + name(s) to/of the species 1 "
                             "orthology \"sif like\" file.", nargs="+", default=None)

    ##### Species 2 Arguments ######
    parser.add_argument("-e2","--expSpecies2",
                        help="The path to the species 2 tab delimited normalized"
                             "read count expression file.")
    parser.add_argument("-d2","--dbSpecies2",
                        help="The path to the species 2 interaction database.",
                        default=None)
    parser.add_argument("-s2","--sifSpecies2",
                        help="The path(s) + name(s) to/of the species 2 sif file(s)",
                        default=None, nargs="+")
    parser.add_argument("-o2","--orthoSpecies2",
                        help="The path(s) + name(s) to/of the species 2 "
                             "orthology \"sif like\" file(s).", nargs="+", default=None)

    ##### Universal Arguments ######
    parser.add_argument( "-p","--pValue",
                        help="The p-Value cutoff for the correlation.",
                        type=float,default=0.05 )
    parser.add_argument( "-oM" ,"--orthoMCL",
                        help="The path + name  to/of the orthoMCL file.",
                        default=None )
    #### Removing####################################################
    #parser.add_argument( "-Xp", "--Xpref",help="\"Strict (s) or loose (l)\" "
    #                                          "cross species network intersect",
    #                     default="s" )
    ######################################################
    parser.add_argument( "-d", "--direction",help=" Option for choosing if you want an output"
                                                 " based on moving from species 1 to species 2 (1to2),"
                                                 "vice versa(2to1) or both(B). ",
                         default="1to2" )
    parser.add_argument("-o","--outfile" ,
                         help=" The name for the output sif file. Note:"
                              " If you choose B for the direction option you must"
                              " provide two file names. (1to2 file name first then 2to1.",
                        nargs="+")

    args = parser.parse_args()

    # Error Checking for Arguments #

    if args.expSpecies1 == None or args.expSpecies2 == None :

        if args.expSpecies1 == None and args.expSpecies2 != None :

            print "ERROR: Species 1 expression file not specified."
            print "Exiting ...."
            sys.exit()

        elif args.expSpecies2 == None and args.expSpecies1 != None :

            print "ERROR: Species 2 expression file not specified."
            print "Exiting ...."
            sys.exit()

        elif args.expSpecies1 == None and args.expSpecies2 == None :

            print "ERROR: Both expression files have not been specified."
            print "Exiting ...."
            sys.exit()

    if args.dbSpecies1 == None and args.sifSpecies1 == None :

        print "WARNING: You have uploaded " \
              "no known network information" \
              " for species 1."

    if args.dbSpecies2 == None and args.sifSpecies2 == None :

        print "Warning: You have uploaded " \
              "no known network information" \
              " for species 2."

    if args.orthoSpecies1 == None and args.orthoSpecies2 == None and args.orthoMCL == None:

        print "ERROR: No Orthology data has been specified."
        print "Exiting ...."
        sys.exit()

    if args.outfile == None:

        print "ERROR: No output file specified."
        print "Exiting ...."
        sys.exit()

    if args.direction == "B":
        if len(args.outfile) == 1:

            print "ERROR: You specified two outputs but " \
                  "only provided one outfile name."
            print "Exiting ...."
            sys.exit()

    ### End of Error Check #####

    main( args.expSpecies1 , args.expSpecies2 ,args.pValue , args.outfile ,
          args.direction.upper() , args.orthoMCL ,args.orthoSpecies1,
          args.orthoSpecies2, args.dbSpecies1 , args.dbSpecies2 ,
          args.sifSpecies1 , args.sifSpecies2 )

