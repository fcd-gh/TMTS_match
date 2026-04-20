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

    



################################################################################
####################Third step: matching with model TS #########################
#after the button is pressed!#
################################################################################


from match_reactions import cleanandclassify, addnewkinetic,check_reaction

run_match=st.button("Find the rates for my reaction")
dic_reaction={"Target":[reactants_list,products_list]}
reaction_ok=check_reaction(dic_reaction)


if run_match and reaction_ok:
    detailed_dic_reaction=cleanandclassify(dic_reaction)
    arr_dic=addnewkinetic(detailed_dic_reaction)
    if len(arr_dic)>0:
        model=list(arr_dic["Target"][3][0].keys())[0]
        A,n,Ea=arr_dic["Target"][3][0][model][0],arr_dic["Target"][3][0][model][1],arr_dic["Target"][3][0][model][2]
        #st.write('Read the reaction in the format of a dictionary:',detailed_dic_reaction)
        st.write('The Arrhenius parameters for this reaction fitted between 500 and 2000K are (in cal, mol,s units):')
        res1,res2=st.columns(2)
        with res1:
            st.write('Pre-exponential parameter:')
            st.write('n parameter:')    
            st.write('Activation energy:')
        with res2:
            st.write(f'{A}')
            st.write(f'{n}')    
            st.write(f'{Ea}')
    else:
        st.write("We are sorry, but this reaction is not covered by our database")
        why_not_covered=st.button('Want to know why?')
        if why_not_covered==True:
            st.write('This part of the code is still under construction :construction:')
elif run_match and reaction_ok==False:
    st.write('We could not find the SMILES for your reaction :sob:')
    





