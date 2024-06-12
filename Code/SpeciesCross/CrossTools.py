#!/usr/bin/env python

"""

This tool contains functions for
 performing cross species
 network analysis as part of
 the CrossSpecies package.

"""

__author__ = 'Gregory Hamilton'
__version__ = '0.0.5'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

############################
# Functions
############################

def filterSif( sifDict , orthoDict ):

    """

    Reduces the sifDict to only interactions that have both
        nodes in the orthoDict.

    :param sifDict: A dictionary sif object.
    :param orthoDict: The orthology dictionary from
    filterOrthoFlatFile().

    :return filteredSif: A sif object reduced to only interactions
      that have both nodes in the orthoDict.

    """

    filteredSif = {}

    for i in sifDict:

        filteredSif[i] = {}

    for node1 in sifDict:

        for node2 in sifDict[node1]:

            for edge in sifDict[node1][node2]:

                if "met:" in edge or edge == "met":

                    try:

                        if orthoDict[node1] != None:

                            filteredSif[node1][node2] = sifDict[node1][node2]

                    except KeyError:

                        try:

                            if orthoDict[node2] != None:

                                filteredSif[node1][node2] = sifDict[node1][node2]

                        except KeyError:

                            continue
                if "met:" not in edge and edge != "met":

                    try:

                        if orthoDict[node1] != None and orthoDict[node2] != None:

                            filteredSif[node1][node2] = sifDict[node1][node2]

                    except KeyError:

                        continue
    removeList = []
    for i in filteredSif:

        if len(filteredSif[i]) < 1 :

            removeList.append(i)

    for i in removeList:

        del filteredSif[i]

    return filteredSif

def strictXIntersectV2(species1sif, species2sif,
                       orthoDict):

    """



    """
    conservedNet = {}

    for n1_sp1 in species1sif:

        for n2_sp1 in species1sif[n1_sp1]:

            for edge_sp1 in species1sif[n1_sp1][n2_sp1]:

                if "pp:" in edge_sp1 or edge_sp1 == "pp":

                    for n1_sp2 in orthoDict[n1_sp1]:

                        for n2_sp2 in orthoDict[n2_sp1]:

                            try:

                                for edge_sp2 in species2sif[n1_sp2][n2_sp2]:

                                    if edge_sp2 == "pp" or "pp:" in edge_sp2:

                                        try:

                                            cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                            cons.append(edge_sp1)
                                            conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                        except KeyError:

                                            conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                    else:

                                        continue

                            except KeyError:

                                try:

                                    for edge_sp2 in species2sif[n2_sp2][n1_sp2]:

                                        if edge_sp2 == "pp" or "pp:" in edge_sp2:

                                            try:

                                                cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                cons.append(edge_sp1)
                                                conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                            except KeyError:

                                                conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                        else:

                                            continue
                                except KeyError:

                                    continue

                elif "corr:" in edge_sp1 and "corr:-" not in edge_sp1:

                    for n1_sp2 in orthoDict[n1_sp1]:

                        for n2_sp2 in orthoDict[n2_sp1]:

                            try:

                                for edge_sp2 in species2sif[n1_sp2][n2_sp2]:

                                    if "corr:" in edge_sp2 and "corr:-" not in edge_sp2:

                                        try:

                                            cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                            cons.append(edge_sp1)
                                            conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                        except KeyError:

                                            conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                    else:

                                        continue

                            except KeyError:

                                try:

                                    for edge_sp2 in species2sif[n2_sp2][n1_sp2]:

                                        if "corr:" in edge_sp2 and "corr:-" not in edge_sp2:

                                            try:

                                                cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                cons.append(edge_sp1)
                                                conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                            except KeyError:

                                                conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                        else:

                                            continue

                                except KeyError:

                                    continue

                elif "corr:-" in edge_sp1:

                    for n1_sp2 in orthoDict[n1_sp1]:

                        for n2_sp2 in orthoDict[n2_sp1]:

                            try:

                                for edge_sp2 in species2sif[n1_sp2][n2_sp2]:

                                    if  "corr:-" in edge_sp2:

                                        try:

                                            cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                            cons.append(edge_sp1)
                                            conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                        except KeyError:

                                            conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                    else:

                                        continue

                            except KeyError:

                                try:

                                    for edge_sp2 in species2sif[n2_sp2][n1_sp2]:

                                        if "corr:-" in edge_sp2:

                                            try:

                                                cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                cons.append(edge_sp1)
                                                conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                            except KeyError:

                                                conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                        else:

                                            continue

                                except KeyError:

                                    continue


                elif "reg:" in edge_sp1 or edge_sp1 == "reg":

                    for n1_sp2 in orthoDict[n1_sp1]:

                        for n2_sp2 in orthoDict[n2_sp1]:

                            try:

                                for edge_sp2 in species2sif[n1_sp2][n2_sp2]:

                                    if "reg:" in edge_sp2 or "reg" == edge_sp2:

                                        try:

                                            cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                            cons.append(edge_sp1)
                                            conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                        except KeyError:

                                            conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                    else:

                                        continue

                            except KeyError:

                                continue

                elif "met:" in edge_sp1 or edge_sp1 == "met":

                    try:

                        for n1_sp2 in orthoDict[n1_sp1]:

                            try:

                                for edge_sp2 in species2sif[n1_sp2][n2_sp1]:

                                    if "met:" in edge_sp2:

                                        try:

                                            cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                            cons.append(edge_sp1)
                                            conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                        except KeyError:

                                            conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                    else:

                                        continue

                            except KeyError:

                                try:

                                    for edge_sp2 in species2sif[n2_sp1][n1_sp2]:

                                        if "met:" in edge_sp2:

                                            try:

                                                cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                cons.append(edge_sp1)
                                                conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                            except KeyError:

                                                conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                        else:

                                            continue

                                except KeyError:

                                    continue

                    except KeyError:

                        try:

                            for n2_sp2 in orthoDict[n2_sp1]:

                                try:

                                    for edge_sp2 in species2sif[n1_sp1][n2_sp2]:

                                        if "met:" in edge_sp2:

                                            try:

                                                cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                cons.append(edge_sp1)
                                                conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                            except KeyError:

                                                conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                        else:

                                            continue

                                except KeyError:

                                    try:

                                        for edge_sp2 in species2sif[n2_sp2][n1_sp1]:

                                            if "met:" in edge_sp2:

                                                try:

                                                    cons = conservedNet[n1_sp1].get(n2_sp1, [])
                                                    cons.append(edge_sp1)
                                                    conservedNet[n1_sp1][n2_sp1] = list(set(cons))

                                                except KeyError:

                                                    conservedNet[n1_sp1] = {n2_sp1: [edge_sp1]}

                                            else:

                                                continue

                                    except KeyError:

                                        continue

                        except KeyError:

                            continue

    return conservedNet


#############################
#############################


if __name__ == '__main__':
    print "Not for commandline use"
