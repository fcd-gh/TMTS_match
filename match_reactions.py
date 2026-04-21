 # -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:06:06 2024

"""

from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions,GetStereoisomerCount
from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters
from rdkit.Chem import AllChem
from _kinetic import iso_alkyl_C3,iso_alkyl_C4,iso_alkyl_C5,iso_alkyl_C6,bs_alkyl_CH,bs_alkyl_CC,habs_H,habs_C,recomb
#import subprocess as sub
import copy
#from anytree import Node, RenderTree, PreOrderIter
import ast
#import os
from rdkit import RDLogger  


def check_reaction(dic_reaction):
    
    
    okay=True
    reactants=dic_reaction["Target"][0]
    products=dic_reaction["Target"][1]
    if len(reactants)==0 or len(products)==0:
        okay=False
    for r in reactants:
        Rrad=Chem.MolFromSmiles(r)
        if Rrad==None or r=='':
            okay=False
    for p in products:
        Prad=Chem.MolFromSmiles(p)
        if Prad ==None or p=='':
            okay=False
    return okay
    

def findmatchHabs(typeofHabs,radreactant, molreactant, radproduct,molproduct):
    """
    Parameters
    ----------
    typeofBS: TYPE str
        DESCRIPTION. either "H" or "CH", define the type of radical abstracting the H (hydrogen = "H" or alkyl = "CH")
    radreactant : TYPE str
        DESCRIPTION. Smiles of the radical abstracting the H
    molreactant : TYPE str
        DESCRIPTION. Smiles of the molecule given the H
    radprodcut: TYPE str
        DESCRIPTION. Smiles of the radical product
    molproduct: TYPE str
        DESCRITPTION. Smiles of the molecular product
    Returns
    -------
    TMTSmodel: TYPE dic
            DESCRIPTION. dictionary cointaining the reaction description as key 
            (format "Rmiles>>Psmile1+Psmile2") and the list of Arrhenius parameters as 
            descriptions [A (s-1),n, Ea(cal/mol) and Lc]
            /!\ Return None if the reaction is not covered by ouyr database
    match: TYPE float
            DESCRIPTION. the fraction of model matching the target reaction
    """
    if typeofHabs=="H":
        kinetic=habs_H
    elif typeofHabs=="CH":
        kinetic=habs_C
    #print(kinetic)
    Rrad=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(radreactant)))
    natomRrad=Rrad.GetNumAtoms()
    Rmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(molreactant)))
    natomRmol=Rmol.GetNumAtoms()
    Prad=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(radproduct)))
    #natomPrad=Prad.GetNumAtoms()
    Pmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(molproduct)))
    #natomPmol=Pmol.GetNumAtoms()
    nmatch=0
    TMTSmodel={}
    wincase=""
    match=0
    for cas in kinetic:
        [reactants,products]=cas.split(">>")
        modelmolR,modelradR=reactants.split("+")
        modelradP,modelmolP=products.split("+")
        RmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(modelmolR)))
        RmodelRad=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(modelradR)))
        PmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(modelmolP)))
        PmodelRad=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(modelradP)))
        nRmol=len(Rmol.GetSubstructMatch(RmodelMol))
        nRrad=len(Rrad.GetSubstructMatch(RmodelRad))
        nPmol=len(Pmol.GetSubstructMatch(PmodelMol))
        nPrad=len(Prad.GetSubstructMatch(PmodelRad))
        ntotal=nRmol+nRrad
        #ntestP=len(Pmolradical.GetSubstructMatch(PalkylmodelMol))+len(Pmolalkene.GetSubstructMatch(PalkenemodelMol))
        if nRmol==nPrad and nRrad==nPmol:
            if nRmol>0 and nRrad>0 and nPmol>0 and nPrad>0 and ntotal>nmatch:
                match=ntotal/(natomRrad+natomRmol)
                nmatch=ntotal
                wincase=cas
    if wincase!='' and wincase not in TMTSmodel:
        TMTSmodel[wincase]=kinetic[wincase]
        return TMTSmodel,match,modelradP
    else:
        return None,0,None
    
def findmatchBS(typeofBS,radicalP,alkeneP,reactant):
    """
    Parameters
    ----------
    typeofBS: TYPE str
        DESCRIPTION. either "H" or "C", define the type of bond scissions (C-C) or (C-H)
    radicalP : TYPE str
        DESCRIPTION. Smiles of the radical product being formed in the B-scission
    alkeneP: TYPE str
        DESCRIPTION. Smiles of the product alkene
    reactant: TYPE str
        DESCRITPTION. Smiles of the alkyl radical reactant
    Returns
    -------
    TMTSmodel: TYPE dic
            DESCRIPTION. dictionary cointaining the reaction description as key 
            (format "Rmiles>>Psmile1+Psmile2") and the list of Arrhenius parameters as 
            descriptions [A (s-1),n, Ea(cal/mol) and Lc]
            /!\ Return None if the reaction is not covered by ouyr database
    match: TYPE float
            DESCRIPTION. the fraction of mo
    """
    Rmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(reactant)))
    Pmolradical=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(radicalP)))
    natomPradical=Pmolradical.GetNumAtoms()
    Pmolalkene=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(alkeneP)))
    natomPalkene=Pmolalkene.GetNumAtoms()
    natomP=natomPradical+natomPalkene
    if typeofBS=="H":
        kinetic=bs_alkyl_CH
    elif typeofBS=="C":
        kinetic=bs_alkyl_CC
    else:
        return None
    nmatch=0


    TMTSmodel={}
    wincase=""
    match=0
    for case in kinetic:
        [reactant,products]=case.split(">>")
        alkylPmodel,alkenePmodel=products.split("+")
        RmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(reactant)))
        PalkylmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(alkylPmodel)))
        PalkenemodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(alkenePmodel)))
        ntestR=len(Rmol.GetSubstructMatch(RmodelMol))
        if typeofBS=="H":
            ntestR+=1
        ntestP=len(Pmolradical.GetSubstructMatch(PalkylmodelMol))+len(Pmolalkene.GetSubstructMatch(PalkenemodelMol))
        if ntestR==ntestP:
            if ntestP>nmatch:
                match=ntestP/natomP
                nmatch=ntestP
                wincase=case
    if wincase not in TMTSmodel:
        TMTSmodel[wincase]=kinetic[wincase]
        #print(TMTSmodel)
        return TMTSmodel,match
    else:
        return None,0
    
def refiningarrheniusBS(TMTSmodel,Rsmiles):
    """

    Parameters
    ----------
    TMTSmodel : TYPE dic
        DESCRIPTION. dictionary cointaining the reaction description as key 
        (format "Rmiles>>Psmiles") and the list of Arrhenius parameters as 
        descriptions [A (s-1),n, Ea(cal/mol) and Lc]
    Rsmiles : TYPE str
        DESCRIPTION. Smiles of the reactant species, it can be canonical or not
    Returns
    -------
    TMTSmodel : TYPE. dic
    DESCRIPTION. corrected with a list of Arrhenius parameters after the correction factor is applied

    """       
    for model in TMTSmodel:
        Lcmodel=TMTSmodel[model][-1]
        fc=1/float(Lcmodel)
        Rmodel,Pmodels=model.split(">>")
        RmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rmodel)))
        nisoRmodel=CalcNumAtomStereoCenters(RmodelMol)+1
        m3= Chem.AddHs(RmodelMol)
        AllChem.EmbedMolecule(m3) 
        AllChem.MMFFOptimizeMolecule(m3)
        rdmolfiles.MolToPDBFile(m3,"C:/pwt/sym.pdb",flavor=16)
        #callplatonsym()
        extsymRmodel=tradsym(foundsym("C:/pwt/sym.lis"))
        RMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rsmiles)))
        nisoR=CalcNumAtomStereoCenters(RMol)+1
        m3R= Chem.AddHs(RMol)
        AllChem.EmbedMolecule(m3R) 
        AllChem.MMFFOptimizeMolecule(m3R)
        rdmolfiles.MolToPDBFile(m3R,"C:/pwt/sym.pdb",flavor=16)
        file = open("C:/pwt/sym.pdb","w+")
        textinput = file.read()
        textinputcorrige="COMPND    C:/pwt/base.xyz\nAUTHOR    GENERATED BY OPEN BABEL 3.1.0\n"+textinput
        file.write(textinputcorrige)
        #callplatonsym()
        extsymR=tradsym(foundsym("C:/pwt/sym.lis"))
        fc=extsymR*nisoRmodel/(extsymRmodel*nisoR)
        Lmodelfinal=fc*float(Lcmodel)
        TMTSmodel[model]=[str(fc*float(TMTSmodel[model][0])),TMTSmodel[model][1],TMTSmodel[model][2],str(Lmodelfinal)]
    return TMTSmodel


def matchreaction_recomb(Rsmiles,Psmile):
    """
    Parameters
    ----------
    Rsmiles : TYPE list of str
        DESCRIPTION. list of the smiles of the reactants
    Psmile : TYPE str
        DESCRIPTION. smile of the product

    Returns
    -------
    Arrhenius : TYPE list
        DESCRIPTION. list of dictionaries containing the reaction TMTS model with the 
        three Arrhenius paramters for the reaction and the statistical factor,
        as the reaction may happen by different pathways
        /!\ if no pathway is identified, then a None will be returned
    """
    Arrhenius={}
    kinetic=recomb
    PMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Psmile)))
    nP=len(PMol.GetAtoms())
    RMol1=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rsmiles[0])))
    nR1=len(RMol1.GetAtoms())
    RMol2=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rsmiles[1])))
    nR2=len(RMol2.GetAtoms())
    ntotal=nP+nR1+nR2
    nmatchb=0
    wincase=""
    for model in  kinetic:
        Rmodelssmi,Pmodelssmi=model.split(">>")
        Rmodelsmi1,Rmodelsmi2=Rmodelssmi.split("+")
        PmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Pmodelssmi)))
        RmodelMol1=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rmodelsmi1)))
        RmodelMol2=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rmodelsmi2)))
        ntestRtotal=max(len(RMol1.GetSubstructMatch(RmodelMol1))+len(RMol2.GetSubstructMatch(RmodelMol2)),len(RMol2.GetSubstructMatch(RmodelMol1))+len(RMol1.GetSubstructMatch(RmodelMol2)))
        ntestP=len(PMol.GetSubstructMatch(PmodelMol))
        nmatchtest=(ntestRtotal+ntestP)/ntotal
        if Rsmiles[1]=="[H]" or Rsmiles[0]=="[H]":
            ntestP+=1
        if nmatchtest>nmatchb and ntestP==ntestRtotal:
            wincase=copy.deepcopy(model)
            nmatchb=nmatchtest
    if nmatchb==1:
        Arrhenius[wincase]=kinetic[wincase]
    elif wincase!="":
        #print(wincase)
        Rmodelssmi,Pmodelssmi=wincase.split(">>")
        Rmodelsmi1,Rmodelsmi2=Rmodelssmi.split("+")
        if Rmodelsmi1==Rmodelsmi2 and Rsmiles[0]!=Rsmiles[1]:
            Arrhenius[wincase]=[str(float(kinetic[wincase][0])*2),kinetic[wincase][1],kinetic[wincase][2]]
        else:
            Arrhenius[wincase]=kinetic[wincase]    
    if Arrhenius=={}:
        Arrhenius=None
    else:
        Arrhenius=[Arrhenius]
    return Arrhenius

def matchreaction_Habs(Rsmiles,Psmiles):
    """
    Parameters
    ----------
    Rsmiles : TYPE list of str
        DESCRIPTION. list containning list of the smiles (str) reactant species, it can be canonical or not
    Psmiles : TYPE list of str
        DESCRIPTION. list containing a list of the smiles (str) of the product species,
        it can be canonical or not

    Returns
    -------
    Arrhenius : TYPE list
        DESCRIPTION. list of dictionaries containing the reaction TMTS model with the 
        three Arrhenius paramters for the reaction and the statistical factor,
        as the reaction may happen by different pathways
        /!\ if no pathway is identified, then a None will be returned
    """
    typeofHabs="CH"
    for R in Rsmiles:
        if "[" and "]" in R:
            radreactant=copy.deepcopy(R)
            if radreactant=="[H]":
                typeofHabs="H"
        else:
            molreactant=copy.deepcopy(R)
    for P in Psmiles:
        if "[" and "]" in P and P!="[HH]":
            radproduct=copy.deepcopy(P)
        elif P=="[HH]":
            molproduct=copy.deepcopy(P)
        else:
            molproduct=copy.deepcopy(P)
    #print(typeofHabs)
    TMTSmodel,match,modelradP=findmatchHabs(typeofHabs,radreactant,molreactant,radproduct,molproduct)
    #print(TMTSmodel)
    #print(match)
    Arrhenius=[]
    #Arrhenius.append(TMTSmodel)
    if match==1:
        Arrhenius.append(TMTSmodel)
    elif TMTSmodel!=None:
        #Arrhenius.append(refiningarrheniusBS(TMTSmodel,Rsmiles))
        fcorr=count_Heq(radproduct)/count_Heq(modelradP)
        for model in TMTSmodel:
            Arrhenius.append({model:[str(float(TMTSmodel[model][0])*fcorr),TMTSmodel[model][1],TMTSmodel[model][2]]})
            #print(TMTSmodel)
    else:
        Arrhenius=None
    return Arrhenius

def matchreaction_BS(Rsmiles,Psmiles):
    """
    Parameters
    ----------
    Rsmiles : TYPE str
        DESCRIPTION. Smiles of the reactant species, it can be canonical or not
    Psmiles : TYPE list of str
        DESCRIPTION. list containint thelist of the smiles of the product species,
        it can be canonical or not

    Returns
    -------
    Arrhenius : TYPE list
        DESCRIPTION. list of dictionaries containing the reaction TMTS model with the 
        three Arrhenius paramters for the reaction and the statistical factor,
        as the reaction may happen by different pathways
        /!\ if no pathway is identified, then a None will be returned
    """
    typeofbs="C"
    rad="[H]"
    alkene=""
    for P in Psmiles:
        if P=="[H]":
            typeofbs="H"
            rad=copy.deepcopy(P)
        if "=" in P:
            alkene=copy.deepcopy(P)
    if typeofbs!="H":
        for P in Psmiles:
            if P != alkene:
                rad=copy.deepcopy(P)

    TMTSmodel,match=findmatchBS(typeofbs,rad,alkene,Rsmiles)
    Arrhenius=[]
    if match==1:
        Arrhenius.append(TMTSmodel)
    elif TMTSmodel!=None:
        Arrhenius.append(refiningarrheniusBS(TMTSmodel,Rsmiles))
    else:
        Arrhenius=None
    return Arrhenius 

def matchreaction_iso(Rsmiles,Psmiles):
    """
    Parameters
    ----------
    Rsmiles : TYPE str
        DESCRIPTION. Smiles of the reactant species, it can be canonical or not
    Psmiles : TYPE str
        DESCRIPTION. Smiles of the product species, it can be canonical or not

    Returns
    -------
    Arrhenius : TYPE list
        DESCRIPTION. list of dictionaries containing the reaction TMTS model with the 
        three Arrhenius paramters for the reaction and the statistical factor,
        as the reaction may happen by different pathways
        /!\ if no pathway is identified, then a None will be returned
    """
    pathways,cs= cyclesize(Rsmiles,Psmiles)
    # print(pathways)
    # print(cs)
    #check if it was possible to find a pathway and a cyclesize:
    if len(cs)==0 or len(pathways)==0:
        return None
    Arrhenius=[]
    for i in range(0,len(cs)):
        cycle=cs[i]
        #print(cycle)
        path=pathways[i]
        TMTSmodel,match= searchmatch(Rsmiles,Psmiles,cycle)
        #print(TMTSmodel)
        if match==1:
            Arrhenius.append(TMTSmodel)
        elif TMTSmodel!=None:
            #here, we need to call the code for calculating the correction factor
            #for the statistical factor
            Arrhenius.append(refiningarrhenius(TMTSmodel,Rsmiles,Psmiles,path,cycle))
        else:
            Arrhenius=None
    return Arrhenius


def refiningarrhenius(TMTSmodel,Rsmiles,Psmiles,path,cycle):
    """

    Parameters
    ----------
    TMTSmodel : TYPE dic
        DESCRIPTION. dictionary cointaining the reaction description as key 
        (format "Rmiles>>Psmiles") and the list of Arrhenius parameters as 
        descriptions [A (s-1),n, Ea(cal/mol) and Lc]
    Rsmiles : TYPE str
        DESCRIPTION. Smiles of the reactant species, it can be canonical or not
    Psmiles : TYPE str
        DESCRIPTION.Smiles of the reactant species, it can be canonical or not

    Returns
    -------
    TMTSmodel : TYPE. dic
    DESCRIPTION. corrected with a list of Arrhenius parameters after the correction factor is applied

    """
    for model in TMTSmodel:
        Lcmodel=TMTSmodel[model][-1]
        fc=1/float(Lcmodel)
        Rsmiles=Rsmiles.replace("@","")
        Rmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rsmiles)))
        nisoR=CalcNumAtomStereoCenters(Rmol)+1
        opts = StereoEnumerationOptions(unique=True)
        isomers = tuple(EnumerateStereoisomers(Rmol, options=opts))
        isomerssmi=[]
        for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
            isomerssmi.append(smi)
        isomerssmi=cleanisomers(isomerssmi)
        NR=len(isomerssmi)
        isomersextsym=[]
        for isomer in isomerssmi:
            moliso=Chem.MolFromSmiles(isomer)
            m3= Chem.AddHs(moliso)
            AllChem.EmbedMolecule(m3) 
            AllChem.MMFFOptimizeMolecule(m3)
            rdmolfiles.MolToPDBFile(m3,"C:/pwt/sym.pdb",flavor=16)
            file = open("C:/pwt/sym.pdb","w+")
            textinput = file.read()
            textinputcorrige="COMPND    C:/pwt/base.xyz\nAUTHOR    GENERATED BY OPEN BABEL 3.1.0\n"+textinput
            file.write(textinputcorrige)
            #callplatonsym()
            extsymiso=tradsym(foundsym("C:/pwt/sym.lis"))
            isomersextsym.append(copy.deepcopy(extsymiso))
        if len(isomersextsym)==1:
            sym=isomersextsym[0]
            LR=[sym/nisoR]
        else:
            LR=[sym/2 for sym in isomersextsym]
        # now we will create a dictionary of TS based on the structures:
        TS=TSdic(path,Rmol)
        subsTS=[TS[i] for i in TS]
        uniquecombofsubs=createTS(subsTS)
        NTS=len(uniquecombofsubs)
        LTSs=[]
        #print(uniquecombofsubs)
        for ts in uniquecombofsubs:
            nisoTS=2
            symTS=1
            #verifyif nisoTS=1
            if cycle <=4 or cycle ==6:
                if samesubextremes(ts)==True:
                    nisoTS=1
            invinvTS=middleCinversion(planinversion(ts))
            if invinvTS==ts:
                symTS=2
            LTS=nisoTS/symTS
            LTSs.append(copy.deepcopy(LTS))
        Lmodel=[]
        i=NR
        for LTS in LTSs:
            j=i%NR
            Lmodeli=LTS*LR[j]
            i+=1
            Lmodel.append(copy.deepcopy(Lmodeli))
        Lmodelfinal=sum(Lmodel)/NR
        fc=fc*Lmodelfinal
        TMTSmodel[model]=[str(fc*float(TMTSmodel[model][0])),TMTSmodel[model][1],TMTSmodel[model][2],str(Lmodelfinal)]
        
    return TMTSmodel

def symfromsmi(smiles):
    """
    

    Parameters
    ----------
    smiles : TYPE str
        DESCRIPTION. smiles of the molecule

    Returns
    -------
    sym: TYPE float
        DESCRIPTION. integer with the symmetry value

    """
    moliso=Chem.MolFromSmiles(smiles)
    m3= Chem.AddHs(moliso)
    AllChem.EmbedMolecule(m3) 
    AllChem.MMFFOptimizeMolecule(m3)
    rdmolfiles.MolToPDBFile(m3,"C:/pwt/sym.pdb",flavor=16)
    file = open("C:/pwt/sym.pdb","r")
    textinput = file.read()
    file.close()
    file2=open("C:/pwt/sym.pdb","w")
    textinputcorrige="COMPND    C:/pwt/base.xyz\nAUTHOR    GENERATED BY OPEN BABEL 3.1.0\n"+textinput
    textinputcorrige=textinputcorrige.replace(" H "," F ")
    file2.write(textinputcorrige)
    file2.close()
    # callplatonsym()
    symiso=tradsym(foundsym("C:/pwt/sym.lis"))
    
# def callplatonsym():
#     newpath = "C:\pwt"
    # os.chdir(newpath)
    # proc = sub.Popen(["platon", "-o", "sym.pdb\n"],
    #                  stdin=sub.PIPE, stdout=sub.PIPE, shell=True, text=True)
    # output2, errors = proc.communicate(input="NONSYM\nexit\n")           

def foundsym(filename):
    symmetry = "C1"
    # symfile = open(filename, "r")
    # text = symfile.read()
    # lines = text.splitlines()
    # for i in range(0, len(lines)):
    #     if lines[i][0:6] == "Resd #":
    #         data = lines[i+2].split()
    #         tol = float(data[9])
    #         if tol == 0.1:
    #             symmetry = str(data[6])
    # symfile.close()
    return symmetry

def tradsym(symmetry):
    extsym = 1
    if symmetry == "C1" or symmetry == "Cs":
        extsym = 1
    elif symmetry == "C2" or symmetry == "C2v":
        extsym = 2
    elif symmetry == "C3v":
        extsym = 3
    elif symmetry == "D2h":
        extsym = 4
    elif symmetry == "D3h" or symmetry == "D3d":
        extsym = 6
    elif symmetry == "D5h":
        extsym = 10
    elif symmetry == "Td":
        extsym = 12
    elif symmetry == "Oh":
        extsym = 24
    return extsym
        
def cleanisomers(listofiso):
    """
    

    Parameters
    ----------
    listofiso : TYPE list
        DESCRIPTION. list of smiles for isomers, chiral centers indicated with @ or @@

    Returns
    -------
    listofisoclean : TYPE list
        DESCRIPTION. list of smiles for isomers, chiral centers indicated with @ or @@, no
                    double structure
    """
    already=[]
    cleaniso=[]
    for iso in listofiso:
        if iso not in already:
            cleaniso.append(iso)
            already.append(iso)
            already.append(inverseiso(iso))
    return cleaniso

def inverseiso(isomer):
    """
    

    Parameters
    ----------
    isomer : TYPE str
        DESCRIPTION.

    Returns
    -------
    isomer : TYPE str
        DESCRIPTION. @ and @@ inverted

    """
    isomer=isomer.replace("@@","XXX")
    isomer=isomer.replace("@","XO")
    isomer=isomer.replace("XXX","@")
    isomer=isomer.replace("XO","@@")
    
    return isomer

def searchmatch(Rsmi,Psmi,C):
    """
    
    Parameters
    ----------
    Rsmi : TYPE str
        DESCRIPTION. SMILES of reactant
    Psmi : TYPE str
        DESCRIPTION. SMILES of the product
    C : TYPE integer
        DESCRIPTION. size of the TScycle

    Returns
    -------
     TMTSmodel : TYPE dic
         DESCRIPTION. dictionary cointaining the reaction description as key 
         (format "Rmiles>>Psmiles") and the list of Arrhenius parameters as 
         descriptions [A (s-1),n, Ea(cal/mol) and Lc]
         /!\ Return None if the reaction is not covered by ouyr database
    """
    Rmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rsmi)))
    natomR=Rmol.GetNumAtoms()
    Pmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Psmi)))
    natomP=Pmol.GetNumAtoms()
    if C==3:
        kinetic=iso_alkyl_C3
    elif C==4:
        kinetic=iso_alkyl_C4
    elif C==5:
        kinetic=iso_alkyl_C5
    elif C==6:
        kinetic=iso_alkyl_C6
    else:
        return None,0
    nmatch=0
    TMTSmodel={}
    wincase=""
    match=0
    for case in kinetic:
        [Rmodel,Pmodel]=case.split(">>")
        RmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Rmodel)))
        PmodelMol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Pmodel)))
        ntestR=len(Rmol.GetSubstructMatch(RmodelMol))
        ntestP=len(Pmol.GetSubstructMatch(PmodelMol))
        if ntestR==ntestP:
            if ntestR>nmatch:
                match=ntestR/natomR
                nmatch=ntestR
                wincase=case
    if wincase not in TMTSmodel:
        TMTSmodel[wincase]=kinetic[wincase]
        return TMTSmodel,match
    else:
        return None,0

def cyclesize(Rsmiles,Psmiles):
    """
    Parameters
    ----------
    Rsmiles : TYPE str
        DESCRIPTION. Smiles of the reactant species, it can be canonical or not
    Psmiles : TYPE str
        DESCRIPTION. Smiles of the product species, it can be canonical or not

    Returns
    -------
    path : TYPE tuple
        DESCRIPTION. tuple containing the indexes of atoms being part of the TS cycle
    cyclesize : TYPE integer
        DESCRIPTION. size of the TS cycle

    """
    Rsmilescanonical=Chem.MolToSmiles(Chem.MolFromSmiles(Rsmiles))
    Psmilescanonical=Chem.MolToSmiles(Chem.MolFromSmiles(Psmiles))
    RMol=Chem.MolFromSmiles(Rsmilescanonical)
    PMol=Chem.MolFromSmiles(Psmilescanonical)
    alkane=Chem.MolFromSmiles(equivalentalkane(Rsmilescanonical))
    ridxR=radatomindex(Rsmilescanonical)
    ridxP=radatomindex(Psmilescanonical)
    path=[]
    cyclesize=[]
    rpositionsR,dicR=equivalencelist(Rsmilescanonical)
    rpositionsP,dicP=equivalencelist(Psmilescanonical)
    for atomR in rpositionsR:
        for atomP in rpositionsP:
            if atomP!=atomR:
                path.append(Chem.rdmolops.GetShortestPath(alkane,atomR,atomP))
    path=cleanpath(path,dicR)
    cyclesize=[len(path[i])+1 for i in range(0,len(path))]
    return path,cyclesize

def cleanpath(path,dicR):
    """
    
    Parameters
    ----------
    path : TYPE list
        DESCRIPTION. list of tuples containing the atoms indexes of a isomerisation TS cycle
    dicR : TYPE dictionary
        DESCRIPTION. tom indexes (keys) with their respective cannonical classification

    Returns
    -------
    path: TYPE list of tuples
        DESCRIPTION. list of tuples containing the atoms indexes of a isomerisation TS cycle,
        but this times cleaned from equivalent ones

    """
    cleanpath=[]
    aux=[]
    for p in path:
        equivalent=[dicR[p[i]] for i in range(0,len(p))]
        if equivalent not in aux:
            cleanpath.append(p)
            aux.append(equivalent)
    return cleanpath

def equivalencelist(radsmiles):
    """
    Parameters
    ----------
    radsmiles : str
        DESCRIPTION. smiles of a radical alkyl

    Returns
    -------
    alleq : list 
        DESCRIPTION. atom indexes of the equivalent positions for the radical center
    dic: dictionary
        DESCRIPTION. atom indexes (keys) with their respective cannonical classification
    """
    nonradsmiles=equivalentalkane(radsmiles)
    Mol=Chem.MolFromSmiles(nonradsmiles)
    listofequivalents=list(Chem.CanonicalRankAtoms(Mol, breakTies=False))
    radMol=Chem.MolFromSmiles(radsmiles)
    radicalidx=radatomindex(radsmiles)
    canradical=listofequivalents[radicalidx]
    alleq=[i for i in range(0,len(listofequivalents)) if listofequivalents[i]==canradical]
    dic={}
    for i in range(0,len(listofequivalents)):
        dic[i]=listofequivalents[i] 
    return alleq,dic 

def radatomindex(radsmiles):
    """

    Parameters
    ----------
    radsmiles : string
        DESCRIPTION. simplified molecular input line entry specification
        it can be canonical or not
        /!\ it must contain only ONE radical site, otherise the function will 
            return None
    Returns
    -------
    indxrad: integer
        DESCRIPTION. index of the atom in the alkane equivalent that contains
        the radical site
    """
    indxradR=-1
    alkane=equivalentalkane(radsmiles)
    Rmol=Chem.MolFromSmiles(radsmiles)
    Amol=Chem.MolFromSmiles(alkane)   
    for atom in Amol.GetAtoms():
        atom.SetNumRadicalElectrons(1)
        if len(Rmol.GetSubstructMatch(Amol))==len(Rmol.GetAtoms()):
            indxradR=atom.GetIdx()
        atom.SetNumRadicalElectrons(0)
    return indxradR

# def radicalatomindex(smiles):
#     """

#     Parameters
#     ----------
#     smiles : string
#         DESCRIPTION. simplified molecular input line entry specification
#         it can be canonical or not
#         /!\ it must contain only ONE radical site, otherise the function will 
#             return None
#     Returns
#     -------
#     radical : integer
#         DESCRIPTION. index of the atom that contains the radical
#         if the species contain moe than 1 radicla site, it will return None
#     """
#     mol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
#     r=Chem.MolToSmiles(mol)
#     radical=[]
#     for atom in mol.GetAtoms():
#         n=atom.GetNumRadicalElectrons()
#         if n>0 :
#             radical.append(atom.GetIdx())
#     if len(radical)==1:
#         return radical[0]
#     else:
#         radical=None
#         return radical
    
def equivalentalkane(radicalsmiles):
    """
    
    Parameters
    ----------
    radicalsmiles : str
        DESCRIPTION.simplified molecular input line entry specification of a radical
        species, it can be canonical or not

    Returns
    -------
    newsmiles2 : str
        DESCRIPTION. simplified molecular input line entry specification
        of the equivalent non-radical species
    
    #important observation: as the Mol in RDkits are C++ objects,
    # it is hard to create a copy of them in python,
    # so I decide to create functions that take smiles all the time 
    # it can be done with xyz coordinates or others, but smiles are 
    # simpler to me
    """
    
    radicalmol=Chem.MolFromSmiles(radicalsmiles)
    radicalmol=Chem.AddHs(radicalmol)
    for atom in radicalmol.GetAtoms():
        # print("numeroatom", atom.GetIdx())
        # print("radical electron", atom.GetNumRadicalElectrons())
        if atom.GetNumRadicalElectrons()!=0:
            #print('esse entrou')
            atom.SetNumRadicalElectrons(0) 
            atom.SetFormalCharge(0)
    sminew=Chem.MolToSmiles(radicalmol)
    newsmiles=Chem.MolToSmiles(Chem.MolFromSmiles(sminew))
    
    return newsmiles

def count_Heq(radsmiles):
    """
    Parameters
    ----------
    radsmiles : str
        DESCRIPTION.simplified molecular input line entry specification of a radical
        species, it can be canonical or not

    Returns
    -------
    nH : int
        DESCRIPTION. number of H in the equivalent position
        
    """
    nH=1
    alkaneeq=equivalentalkane(radsmiles)
    Molcomplete=Chem.AddHs(Chem.MolFromSmiles(alkaneeq))
    indrad=radatomindex(radsmiles)
    for bond in Molcomplete.GetBonds():
        if bond.GetBeginAtomIdx()==indrad:
            if Molcomplete.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomicNum()==1:
                Hpos=int(bond.GetEndAtomIdx())
            
    #Radcomplete=Chem.AddHs(radMol)
    # for atom in Molcomplete.GetAtoms():
    #     if atom.GetAtomicNum()==1:
    #         print('found an H')
    #         with Chem.RWMol(Molcomplete) as mw:
    #             mw.RemoveAtom(atom.GetIdx())
    #             moltest=Chem.MolFromSmiles(Chem.MolToSmiles(mw))
    #         print(Chem.MolToSmiles(Chem.RemoveHs(moltest)))
    #         print(Chem.MolToSmiles(radMol))
    #         if len(moltest.GetSubstructMatch(radMol))==len(radMol.GetAtoms()):
    #             print('found THE H')
    #             Hpos=int(atom.GetIdx())
    
    
    listofequivalents=list(Chem.CanonicalRankAtoms(Molcomplete, breakTies=False))
    nH=listofequivalents.count(listofequivalents[Hpos])   
    
    return nH

def getnonradicalfrag(smiles):
    """
    

    Parameters
    ----------
    smiles : TYPE str
        DESCRIPTION. smiles containning a fragmented radical

    Returns
    -------
    smi : TYPE str
        DESCRIPTION. smiles of the non radical fragment

    """
    
    listofsmiles=smiles.split(".")
    for smi in listofsmiles:
        if "[" not in smi and "]" not in smi:
            return smi

def TSdic(path,R1):
    """
    
    Parameters
    ----------
    path : TYPE  tuple
        DESCRIPTION. indexes of the atoms in the TS cycle
    R1 : TYPE Mol object, from RDkit
        DESCRIPTION. radical reactant

    Returns
    -------
    TSdic : TYPE dic
        DESCRIPTION. dictionary containing the substituents of each atom on the cycle, follow the format:
            TSDic={ind:[0,0],id:["CCC",0]}..

    """
    TSdic={}
    for a in path: 
        bonds=[]
        fragsa=[]
        atom=R1.GetAtomWithIdx(a)
        #print("atom ",a)
        for x in atom.GetNeighbors():
            if x.GetIdx() not in path:
                #print("bonded to",x.GetIdx(),"with bond number", R1.GetBondBetweenAtoms(atom.GetIdx(),x.GetIdx()).GetIdx())
                bonds.append(R1.GetBondBetweenAtoms(atom.GetIdx(),x.GetIdx()).GetIdx())
        if len(bonds)>0:
            frag=Chem.FragmentOnSomeBonds(R1,tuple(bonds),addDummies=False)
            for chem in frag:
                chemfrag=Chem.MolToSmiles(chem)
                newchemfrag=getnonradicalfrag(chemfrag)
                fragsa.append(copy.deepcopy(newchemfrag))
        else:
            fragsa=[0,0]
        if len(fragsa)==2:
            TSdic[a]=fragsa
        else:
            fragsa.append(0)
            TSdic[a]=fragsa
    return TSdic
    
def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    # for bond in mol.GetBonds():
    #     bond.SetBondMapNum(bond.GetIdx())
    return mol


def creer_solutions(liste,pivot,arbre,longueur):
    if pivot < longueur:
        filsg = Node(liste[pivot],arbre)
        filsd = Node(list(reversed(liste[pivot])),arbre,)
        pivot += 1
        creer_solutions(liste,pivot,filsg,longueur)
        creer_solutions(liste,pivot,filsd,longueur)

def planinversion(subslist):
    newlist=[list(reversed(i)) for i in subslist]
    return newlist

def middleCinversion(sublist):
    return list(reversed(sublist))

def createTS(liste):
    liste.insert(0,0)
    pivot=1
    root = Node(liste[0])
    longueur=len(liste)
    creer_solutions(liste, pivot, root, longueur)
    cnodes = list(PreOrderIter(root, filter_=lambda node: node.is_leaf))
    combinaisons = [str(i)[7:-2].split(i.separator) for i in cnodes]
    combinaisons=[comb[1:] for comb in combinaisons]
    c=[]
    for comb in combinaisons:
       comb=[ast.literal_eval(i) for i in comb]
       c.append(copy.deepcopy(comb))
    uniquesubs=[]
    for comb in c:
        if len(comb)<=3:
            combinvplan=planinversion(comb)
            combinvC=middleCinversion(comb)
            doubleinv=planinversion(combinvC)
            if (combinvplan not in uniquesubs) and (combinvC not in uniquesubs) and (comb not in uniquesubs) and (doubleinv not in uniquesubs):
                uniquesubs.append(copy.deepcopy(comb))
        elif len(comb)>=4:
            combinvC=middleCinversion(comb)
            if combinvC not in uniquesubs and (comb not in uniquesubs):
                uniquesubs.append(copy.deepcopy(comb))
    return uniquesubs

def samesubextremes(subsTS):
    """
    

    Parameters
    ----------
    subsTS : TYPE list
        DESCRIPTION.

    Returns
    -------
    same type Boolean
    

    """
    p1=subsTS[0]
    p2=subsTS[-1]
    if p1==list(reversed(p1)) and p2==list(reversed(p2)):
        return True
    else:
        return False

def testforaromaticatom(speciessmiles):
    """
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an species

    Returns
    -------
    aromatic : TYPE boolean
        DESCRIPTION. True if there is an atom in a aromatic ring,
                    False if no aromatic atom is found

    """
    aromatic=False
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    for atom in rdkitmol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic=True  
    
    return aromatic

def testforheteroatom(speciessmiles):
    """
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an species

    Returns
    -------
    heteroatom : TYPE boolean
        DESCRIPTION. True if there is an heteroatom (atom different than C and H),
                    False if no heteroatom is found

    """
    heteroatom=False
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    for atom in rdkitmol.GetAtoms():
        if atom.GetAtomicNum()!=6 and atom.GetAtomicNum()!=1:
            heteroatom=True  
    
    return heteroatom

def testcharge(speciessmiles):
    """
    
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an specie

    Returns
    -------
    charged : TYPE Boolean
        DESCRIPTION. True if the charge of the specie is different than 0
                    False if the charge is 0

    """
    charged=False
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    for atom in rdkitmol.GetAtoms():
        if atom.GetFormalCharge()!=0:
            charged=True
    return charged

def testforcycle(speciessmiles):
    """
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an species

    Returns
    -------
    cycle : TYPE boolean
        DESCRIPTION. True if there is a ring,
                    False if no ring is found

    """
    cycle=False
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    for bond in rdkitmol.GetBonds():
        if bond.IsInRing() :
            cycle=True  
    
    return cycle

def testforbonds(speciessmiles):
    """
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an species

    Returns
    -------
    bonds : TYPE boolean
        DESCRIPTION. True if there is a triple, dative or aromatic bond
                    False if the bonds are only single and doubles

    """
    bonds=False
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    for bond in rdkitmol.GetBonds():
        if str(bond.GetBondType())!="DOUBLE" and str(bond.GetBondType())!="SINGLE":
            bonds=True  
    return bonds

def classreaction(rside,pside):
    """

    Parameters
    ----------
    rside : TYPE list
        DESCRIPTION.  list of stringscontaining the type os species (it can be "alkyl","alkene" and "alkane")
    pside : TYPE list
        DESCRIPTION. list of stringscontaining the type os species (it can be "alkyl","alkene" and "alkane")

    Returns
    -------
    reactiontype : TYPE str
        DESCRIPTION. string definin the type of reaction "isoCH","BSCH","1/BSCH","recombCH","1/recombCH" and "HabsCH"

    """
    typeofreaction=None
    if None in pside or None in rside:
        typeofreaction=None
    elif pside==rside==["alkyl"]:
        typeofreaction="isoCH"
    elif sorted(pside)==['alkene', 'alkyl'] and rside==["alkyl"]:
        typeofreaction="BSCH"
    elif sorted(rside)==['alkene', 'alkyl'] and pside==["alkyl"]:
        typeofreaction="1/BSCH"
    elif rside==["alkyl","alkyl"] and pside==["alkane"]:
        typeofreaction="recombCH"
    elif pside==["alkyl","alkyl"] and rside==["alkane"]:
        typeofreaction="1/recombCH"
    elif sorted(pside)==['alkane', 'alkyl'] and sorted(rside)==['alkane', 'alkyl']:
        typeofreaction="HabsCH"
        
    return typeofreaction


def classspecie(speciessmiles):
    """
    Parameters
    ----------
    speciessmiles : TYPE str
        DESCRIPTION. smiles of an species

    Returns
    -------
    speciestype : TYPE str or None
        DESCRIPTION. the type of species, can be classified as:
                    "alkyl" = only C and H atoms, one radical point, no double bonds, no cycle
                    "alkane"= only C and H atoms, no radical point, no double bond, no cycle
                    "alkene"=only C and H atom, no radical point, one double bond, no cycle
                    
    """
    typeofspecie=None
    rdkitmol=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(speciessmiles)))
    radical=0
    doublebond=0
    testar=testforaromaticatom(speciessmiles)
    testhet=testforheteroatom(speciessmiles)
    testcha=testcharge(speciessmiles)
    testcycle=testforcycle(speciessmiles)
    testbonds=testforbonds(speciessmiles)
    notincluded=True
    if testar==False and testhet==False and testcha==False and testcycle==False and testbonds==False:
        for b in rdkitmol.GetBonds():
            if str(b.GetBondType())=="DOUBLE":
                doublebond+=1
        for atom in rdkitmol.GetAtoms():
            if atom.GetNumRadicalElectrons()==1:
                radical+=1
        if doublebond==1 and radical==0:
            typeofspecie="alkene"
        elif doublebond==0 and radical==0:
            typeofspecie="alkane"
        elif doublebond==0 and radical==1:
            typeofspecie="alkyl" 
    else:
        typeofspecie=None 
    return typeofspecie

def cleanandclassify(dicofreaction):
    """
     Parameters
    ----------
    dicofreaction : TYPE dic
        DESCRIPTION. a dictionary of the the reactions, the key being the reaction as written in the mechanism,
        linked to a list of list of reactants followed by a list of products in the format:
            dic ={"C8H18-1=>R4CH3+R34C7H15":[["C(C)(C)(C)CC(C)C"],["[CH3]","[C](C)(C)CC(C)C"]]}
            /!\ it only includes reactions with smiles equivalents

    Returns
    -------
    newdicofreaction : TYPE dic
        DESCRIPTION. similar to the dictionary of the reactions, but the reactions not covered by our reaction 
        classes (isomerisations, Beta-scissions, H-abstraction and recombinations of alkyl radicals and alakenes)
        are excluded. Also, the type of the reaction is added as a string in the list:
            dic ={"C8H18-1=>R4CH3+R34C7H15":[["C(C)(C)(C)CC(C)C"],["[CH3]","[C](C)(C)CC(C)C"],"1/recomb"],}
            possible types of reactions listed as:
                "isoCH"=isomerisations of alkyl radicals
                "BSCH"=betascissions of alkyl radicals
                "1/BSCH"=inverse of BS, addition to double bond
                "recombCH"=recombination of two alkyl radicals
                "1/recombCH"=inverse of recombination of twoalkyl radicals, also known as unimolecular initiation
                "HabsCH"=H-abstraction of alkane by alkyl/H radicals
    """
    #0th thing: test if len (Rs) or len(Ps)>2
    #1st thing: test if there are heteroatoms in the R or P list
    #second thing: test if there are aromatics in R or P list 
    #reflection: maybe classify the type of molecule and radical?Fabi
    
    newdicofreaction={}
    for reaction in dicofreaction:
        excluded=False
        if len(dicofreaction[reaction][0])>2 or len(dicofreaction[reaction][1])>2:
            excluded=True
        else:
            rside=[]
            pside=[]
            for R in dicofreaction[reaction][0]:
                rside.append(classspecie(R))
            for P in dicofreaction[reaction][1]:
                pside.append(classspecie(P))
            reactiontype=classreaction(rside,pside)
            if reactiontype!=None:
                newdicofreaction[reaction]=[dicofreaction[reaction][0],dicofreaction[reaction][1],reactiontype]
    return newdicofreaction

def addnewkinetic(dicofreactions):
    """
    Parameters
    ----------
    dicofreactions : TYPE dic
        DESCRIPTION. dictionaryof reactions as the ones from cleanandclassify function, in the format:
            {'C8H18-1=>R4CH3+R34C7H15': [['C(C)(C)(C)CC(C)C'], ['[CH3]', '[C](C)(C)CC(C)C'], '1/recombCH'],..}

    Returns
    -------
    updateddic : TYPE dic
        DESCRIPTION. dictionary containint the reactions as key, linked with a list containing the reactants smiles,
        the product smiles, the type of reaction, and a list with the Arrhenius data. 

    """
    updateddic={}
    niso=0
    nBS=0
    nrecomb=0
    nhabs=0
    for reaction in dicofreactions:
        typeofreaction=dicofreactions[reaction][2]
        if typeofreaction=="isoCH":
            arrhenius=matchreaction_iso(dicofreactions[reaction][0][0],dicofreactions[reaction][1][0])
            if arrhenius!=None:
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
                niso+=1
        elif typeofreaction=="BSCH":
            arrhenius=matchreaction_BS(dicofreactions[reaction][0][0],dicofreactions[reaction][1])
            if arrhenius!=None:
                nBS+=1
                for model in arrhenius:
                    for m in model:
                        model[m]=[model[m][0],model[m][1],str(float(model[m][2])*1000)]
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
        elif typeofreaction=="1/BSCH":
            arrhenius=matchreaction_BS(dicofreactions[reaction][1][0],dicofreactions[reaction][0])
            if arrhenius!=None:
                nBS+=1
                for model in arrhenius:
                    for m in model:
                        model[m]=[model[m][0],model[m][1],str(float(model[m][2])*1000)]
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
        elif typeofreaction =="recombCH":
            arrhenius=matchreaction_recomb(dicofreactions[reaction][0],dicofreactions[reaction][1][0])
            if arrhenius!=None:
                nrecomb+=1
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
        elif typeofreaction =="1/recombCH":
            arrhenius=matchreaction_recomb(dicofreactions[reaction][1],dicofreactions[reaction][0][0])
            if arrhenius!=None:
                nrecomb+=1
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
        elif typeofreaction =="HabsCH":
            nhabs+=1
            arrhenius=matchreaction_Habs(dicofreactions[reaction][0],dicofreactions[reaction][1])
            if arrhenius!=None:
                for model in arrhenius:
                    for m in model:
                        model[m]=[model[m][0],model[m][1],str(float(model[m][2])*1000)]
                updateddic[reaction]=dicofreactions[reaction]
                updateddic[reaction].append(arrhenius)
    # print("niso:",niso)
    # print("nbs:",nBS)
    # print("nrecomb:", nrecomb)
    # print("nhabs",nhabs)
    return updateddic

RDLogger.DisableLog('rdApp.*') #this is only for avoiding RdKit warnings
# examples of how to use those functions:
# v
#arr=matchreaction_BS("CC(C)(C)CC([CH2])CC(C)(C)C",["[H]","CC(C)(C)CC(=C)CC(C)(C)C"])
#arr=matchreaction_iso("C[CH]C","CC[CH2]")
#arr=matchreaction_recomb(["CC(C)(C)CC([CH2])CC(C)(C)C","[H]"],'CC(C)(C)CC(C)CC(C)(C)C')
# important comment: make sure to check units before using.
