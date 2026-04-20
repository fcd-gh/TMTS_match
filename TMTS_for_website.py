# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:36:24 2026

@author: citrango1
"""

# import subprocess as sub
# import copy
# from anytree import Node, RenderTree, PreOrderIter
# import ast,copy
# import os
# from rdkit import *
# from rdkit import Chem
# from rdkit.Chem import rdmolfiles
# from rdkit.Chem import Draw
# from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions,GetStereoisomerCount
# from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters
import streamlit as st
from streamlit_ketcher import st_ketcher

################################################################################
#### Rate constants estimations using TMTS approach#############################
################################################################################

##### First step: title and infos#########################################
#st.image('Image1.JPG',width=850)

st.title('Estimation of rate constants for common reactions in pyrolysis and combustions')
st.subheader('Method available for isomerisations and beta-scissions, of alkyl radicals.'+
            ' H abstractions will be available soon')
st.markdown("***F.C. Destro, R. Fournet, R. Bounaceur, P.A.Glaude, B. Sirjean***")
st.markdown("***Université de Lorraine, CNRS, LRGP, F-54000 Nancy, France***")
st.write("----------------------------------------------------------")

st.write('For more information you dowload the related papers:\n '+
        ' [Isomerisations](https://linkinghub.elsevier.com/retrieve/pii/S0010218024004413)\n'+
        ' [Beta-scissions](https://linkinghub.elsevier.com/retrieve/pii/S0010218025004821)')

st.write("----------------------------------------------------------")
#st.image('Image2.JPG',width=500)

################################################################################
####################Second step: retrieving the smiles##########################
################################################################################

SMILES_REACTION = 'C[CH-]C>>CC[CH2-]' # Défaut reaction

DEFAULT_REACTION = SMILES_REACTION

reaction= st.text_input("Insert the SMILE notation of the reaction or draw it", DEFAULT_REACTION)

 
REQUESTED_REACTION = st_ketcher(reaction)
st.markdown(f"Smile code for reaction: ``{REQUESTED_REACTION}``")

# reaction=reaction.replace('-','')
# reactants,products=reaction.split(">>")
# reactants_list=reactants.split(".")
# products_list=products.split(".")

# dic_reaction={"Target":[reactants_list,products_list]}

################################################################################
####################Third step: matching with model TS##########################
################################################################################


from match_reactions import cleanandclassify, matchreaction_BS, matchreaction_recomb, matchreaction_Habs, matchreaction_iso,addnewkinetic

#detailed_dic_reaction=cleanandclassify(dic_reaction)
#st.write('Read the reaction in the format of a dictionary:',detailed_dic_reaction)
#arr_dic=addnewkinetic(detailed_dic_reaction)






