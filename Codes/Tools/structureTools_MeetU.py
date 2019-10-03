#!/usr/bin/env python

import math, string



########################################
#         PARSING TOOLS
########################################



def parsePDBMultiChains(infile, charge = 1, chargeFromInfile = False, bfactor = False, CG = False) :
    """ purpose: to parse a pdb file (infile)
        input: PDB file
        output: a dico dPDB which contains for each atom of each residue of
        each chain, its corresponding 3D coordinates. Please take a look to
        the code to understand the structure of the dico.

    """

    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["reslist"] = []
    dPDB["chains"] = []
    
    # parcoure le PDB   
    for line in lines :
        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and ( (string.strip(line[17:20]) == "MET") or  (string.strip(line[17:20]) == "MSE") )) :
            chain = line[21]
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s"%(line[22:27]).strip()
            resnum = "%s"%(line[22:26]).strip()
            if not curres in dPDB[chain]["reslist"] : # first time we encounter it
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = string.strip(line[17:20])
                dPDB[chain][curres]["atomlist"] = []
                #dPDB[chain][curres]["atomlistTowrite"] = []
                alternateoccupancy = None #"%s"%(line[16:17])
                occupancy = "%s"%(line[16:17]) 
                if occupancy != " " :
                    alternateoccupancy = occupancy
                

            else: # this is not a new residue
                occupancy = "%s"%(line[16:17])

                if occupancy != " " and alternateoccupancy == None : # means we are in the first alternate location of that residue
                    alternateoccupancy = occupancy
            
            if CG : # means we are parsing a CG model so we have to treat the CSE atomtypes which can be redondant in terms of name the same res
                atomtype = "%s_%s"%(string.strip(line[6:11]), string.strip(line[12:16]))
            else:
                atomtype = string.strip(line[12:16])
            
            #if not atomtype in dPDB[chain][curres]["atomlist"] :
            if occupancy == alternateoccupancy  or occupancy == " " : # means this atom corresponds to the first rotamer found in the PDB for this residue

                #if CG :
                    #dPDB[chain][curres]["atomlistTowrite"].append(atomtype.split("_")[1]) # necessary for the writing later
                
                dPDB[chain][curres]["atomlist"].append(atomtype)
                dPDB[chain][curres][atomtype] = {}
                dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
                dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
                dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
                dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
                if bfactor == True :
                    dPDB[chain][curres][atomtype]["bfactor"] = float(line[60:67].strip())
                #if chargeFromInfile == True :
                #    dPDB[chain][curres][atomtype]["charge"] = float(line[60:67])
                #else :
                #    dPDB[chain][curres][atomtype]["charge"] = charge


            dPDB[chain][curres]["resnum"] = resnum
            #dPDB[chain][curres]["inser"] =  "%s"%(line[26:27])
 

    return dPDB


#################################################
#           WRITING TOOLS
#################################################


def writePDB(dPDB, filout = "out.pdb", bfactor = False) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each residue (one key per residue) in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"] :
            for atom in dPDB[chain][res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"],dPDB[chain][res]["bfactor"] ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"] ))
                    
    fout.close()


def initBfactor(dPDB):
    """purpose: initiation of the bfactor key for each residue
       input: a dico with the dPDB format
    """

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"]:
            dPDB[chain][res]["bfactor"] = 0


            
def generateFastPDB(x, y, z, res = "GEN", atomname = "X", atomid = 1, resid = 1, chain = " ", bfactor = ""):
    """ //// DEBUG FUNCTION ////
        purpose: creates a mini dico dPDB for one atom and its 3D coordinates.
        The idea is to call after the writePDB(my_mini_dico) in order to visualize
        with Pymol the coordinates of the corresponding atom.
        input: x, y, z (3D coordinates of the atom we want to visualize)
        output: a mini dPDB dico for one atom
        usage: my_mini_dico = generateFastPDB(xi, yi, zi) 

    """

    dPDB = {}
    dPDB["chains"] = [chain]
    dPDB[chain] = {}
    dPDB[chain]["reslist"] = [resid]
    dPDB[chain][resid] = {}
    dPDB[chain][resid]["atomlist"] = [atomname]
    dPDB[chain][resid][atomname] = {}
    dPDB[chain][resid][atomname]["id"] = atomid
    dPDB[chain][resid]["resname"] = res
    dPDB[chain][resid][atomname]["x"] = x
    dPDB[chain][resid][atomname]["y"] = y
    dPDB[chain][resid][atomname]["z"] = z
    if bfactor != "":
        dPDB[chain][resid][atomname]["bfactor"] = bfactor

    return dPDB



