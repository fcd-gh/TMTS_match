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

st.title('Estimation of rate constants for common reactions in pyrolysis and combustions :fire:')
st.subheader('Method available for isomerisations and beta-scissions, of alkyl radicals.'+
            ' H-abstractions will be available soon')
st.markdown("***F.C. Destro, R. Fournet, R. Bounaceur, P.A.Glaude, B. Sirjean***")
st.markdown("***Université de Lorraine, CNRS, LRGP, F-54000 Nancy, France***")
st.write("----------------------------------------------------------")

st.write('For more information, you can download the related papers:\n ')
st.write(':newspaper: [Isomerisations](https://linkinghub.elsevier.com/retrieve/pii/S0010218024004413)\n')
st.write(':newspaper: [Beta-scissions](https://linkinghub.elsevier.com/retrieve/pii/S0010218025004821)')

st.write("----------------------------------------------------------")
#st.image('Image2.JPG',width=500)

################################################################################
####################Second step: retrieving the smiles##########################
################################################################################

# SMILES_REACTION = 'C[CH-]C>>CC[CH2-]' # Défaut reaction

# DEFAULT_REACTION = SMILES_REACTION
st.write(' You can either enter your reaction by inserting the SMILES of reactants and products \n in the following boxes')
col1, col2, col3, col4, col5 = st.columns([2,2,1,2,2], gap="small")
with col1:
    react1=st.text_input('Reactant 1')
with col2:
    react2=st.text_input('Reactant 2 (optional)')
with col3:
    st.write(':arrow_right:')
with col4:
    prod1=st.text_input('Product 1')
with col5:
    prod2=st.text_input('Product 2 (optional)')
reactants_list=[r for r in [react1,react2] if r!='']
products_list=[p for p in [prod1,prod2] if p!='']

st.write ('Or you can draw it in the following box:')
reaction = st_ketcher()
st.markdown(f"Smile code for reaction: ``{reaction}``")
if '-' in reaction:
    reaction=reaction.replace('-','')
if '>>' in reaction:
    reactants,products=reaction.split(">>")
    reactants_list=reactants.split(".")
    products_list=products.split(".")
st.write('Reactant(s)',reactants_list)
st.write('Product(s)',products_list)
    
st.write(dic_reaction={"Target":[reactants_list,products_list]})

################################################################################
####################Third step: matching with model TS##########################
################################################################################


from match_reactions import cleanandclassify, matchreaction_BS, matchreaction_recomb, matchreaction_Habs, matchreaction_iso,addnewkinetic

#detailed_dic_reaction=cleanandclassify(dic_reaction)
#st.write('Read the reaction in the format of a dictionary:',detailed_dic_reaction)
#arr_dic=addnewkinetic(detailed_dic_reaction)






