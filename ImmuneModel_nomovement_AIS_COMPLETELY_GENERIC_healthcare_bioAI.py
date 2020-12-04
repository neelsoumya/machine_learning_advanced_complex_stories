# -*- coding: utf-8 -*-
"""
@author: Soumya Banerjee
"""

#################################################################################################
#
# Description -  COMPLETELY GENERIC SIMULATOR SYNTHETIC DATA
#               Stochastic simulator to simulate T-cell negative
#                       selection in the thymus (no movement)
#               FOR AIS data
#               NOTE: has withincell turnstile model i.e. if 2 or more peptides match on a single
#                    interaction (TCR-pMHC), then negatively select
#               generic takes in En, num synpases and num mtecs as arguments,
#                 will do turnstile/threshold model, checkauto and also fixed TCR
#                       'yescheckauto'/'nocheckauto' 'turnstile'/'threshold' fixed TCR sequence 'NNLFF'/'notcrseq'
#
# Installation requirements
#       Install using
#       pip install mesa
#
# Usage:
#      from root of repo
#       python ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI  tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI 1000 False False 100 False   50 1000 2 nocheckauto turnstile notcrseq
#
#       python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI  tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI 1000 False False 100 False   50 1000 2 nocheckauto turnstile notcrseq
#
#       python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py 'tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 'tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 10 False False 100 False   70 1000 100 'yescheckauto' 'turnstile' '1001010010100101001010010100101001010010100101001010010100101001010010100101001010010100101001010010' > temp.out
#       python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py 'tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 'tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 10 False False 100 False   70 1000 100 'nocheckauto' 'turnstile' 'notcrseq' > temp.out
#       python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py 'tcell_survival_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 'tcell_matches_ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI' 10 False False 100 False   70 1000 100 'nocheckauto' 'threshold' 'notcrseq' > temp.out
#       nohup python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC_healthcare_bioAI.py.py 'm_ImmuneModel_nomovement_SYNTHETIC_turnstilewithincell_yescor_noauto_tcrandomunif_1e4tcell_en21_syn1000_tec174_c100iter' 's_ImmuneModel_nomovement_SYNTHETIC_turnstilewithincell_yescorr_noauto_tcrandomunif_1e4tcell_en21_syn1000_tec174_c100iter' 10000 False False 100 False   70 1000 174  > /dev/null
#
#       auto is True    b_auto_TCR_sequence_given = False pull from auto TCR pool
#       nohup python3 ImmuneModel_nomovement_AIS_COMPLETELY_GENERIC.py.py 'autopool_ImmuneModel_nomovement_SYNTHETIC_turnstilewithincell_AIRE_5aa_1e5tcell_c100iter_1' 'autopool_ImmuneModel_nomovement_SYNTHETIC_turnstilewithincell_AIRE_5aa_1e5tcell_c100iter' 100000 False False 100 True   70 1000 174  > /dev/null
#
# Testing: from root of repo
#       ./automated_doc_testing.sh
#
#
# Author  - Soumya Banerjee
# Website - https://sites.google.com/site/neelsoumya/
#
################################################################################################


####################################################################
# Import required libraries
####################################################################

import matplotlib.pyplot as plt
import numpy as np
import time
import pandas as pd
import os
import csv
import sys
import pdb # for debugging
# from functools import lru_cache
# import stochpy # for stochastic simulation algorithms
# from custom_memory_utilities import * # memory profiling, freeing memory
# from scipy.stats import genextreme
# import logging
# import unittest
# import coverage
# import pickle
# from joblib import Parallel, delayed
# import copy
# import warnings
# import pyproj
# import seaborn as sns
# from abc import ABCMeta, abstractmethod
from functions_lib_new import *


####################################################################
# Parameters for running simulation MODULE
####################################################################
# TODO: read model parameters from an model initialization (text) file
#       filename passed in from function arguments in main
# Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/27/read-input-model-parameters-from-file

# from model_input_parameters import *
i_num_tcells         = 1   # number of T-cells tested in each Monte-Carlo round
#i_num_mtecs          = 174
#i_num_interactions   = 192 # 4 days residence time (96 hours) * 0.5 hours for each contact (http://www.sciencedirect.com/science/article/pii/S0960982200008708) Fig. 1 
# Current Biology, Volume 10, Issue 24, 14 December 2000, Pages R923–R933, Jérôme Delona, Ronald N Germain
# Fig. 1 caption mentions 30 minutes
i_num_montecarlo     = 1 # number of times simulations to be repeated
b_auto_TCR_sequence_given = False # True, if want to try a given auto sequence. False pick randomly from offline file


#######################################################
# Required model parameters and variables MODULE
#######################################################

# Array of amino acid sequences
# from https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
str_single_aa = ['1', '0']
f_aa_mice_proteome_frequencies = [0.5, 0.5]

i_single_aa   = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                 13, 14, 15, 16, 17, 18, 19, 20]
# Array of amino acid sequences that are hydrophobic
#str_hydrophobic_aa = ['V', 'I', 'L', 'F', 'W', 'Y', 'M']
# Frequency of aa for Listeria monocytogenes
f_pathogen_frequencies = [0.5, 0.5]

# change length of TCR and peptide here
i_length_TCR = 74 #100  # Length of T-cell receptor
i_length_peptide = 74 #100 # length of mTEC presented peptide

#i_num_peptide_present = 1000 # each cell presents 1000 peptides from bag
# if TCR matches 5 or more, then kill thymocyte
#i_matching_threshold_within_peptide = 70 # number of letters/amino acids/symbols that must match for negatively selection in a single peptide
i_matching_threshold = 2 # if 2 or more peptides match on a single interaction (TCR-pMHC), then negatively select
#f_negative_selection_threshold = -21.0 # -21.0 from paper with 5 (?) aa
f_positive_selection_threshold = -18.5

tcell_death_probability = 0.1 # Tcell natural death or exit from thymus probability
mtec_death_probability  = 0.1 # mTEC natural death or exit from thymus probability

list_tcell_count   = []   # list to keep track of T-cell counts at each time step
list_tcell_matches = [] # list to keep track of number of positive matches TCR-pMHC
list_all_tcell_sequence_escape = [] # list of all TCR sequences that escaped thymic selection
list_all_tcell_sequence_killed = [] # list of all TCR sequences that are killed/deleted during thymic selection
list_all_tcell_sequence_escape_BUTAUTO = [] # list of all TCR sequences that escaped thymic selection BUT ARE AUTO


######################################################
# Directory parameters
######################################################
#DATA_DIR = "data/copy0_filter_aire" # copy >0 and AIRE and AIRE-dep only
#DATA_DIR = "data" # all data all genes copy > 0
#DATA_DIR = "data/copy_EQ_0_filter_aire_zeroinf" # copy >= 0 and AIRE and AIRE-dep only and zero inflated
#DATA_DIR = "data/copy_EQ_0_filter_aire_zeroinf_transcript" 
#   copy >= 0 and AIRE and AIRE-dep only and zero inflated and transcript abundance
#DATA_DIR = "data/copy0_filter_NOTaire" # copy >0 and NOT AIRE (AIRE-)
DATA_DIR = "synthetic_data" # synthetic data generated
# TODO: make generic output directory that can be used to run on different machines
str_output_directory = "/scratch/banerjees/agent_based_model_mesa/gambit"

######################################################
# Performance parameters
######################################################
#LRU_CACHE_SIZE = 10000


######################################################
# Load the MJ Miyazawa-Jernigan interaction matrix
######################################################
from static_protein_files import matrix_mj

# matrix_JJernigan(TCR, TEC_peptide) static dict or static method
# matrix_JJernigan(TCR, TEC_peptide) will return an energy
# thresholds for positive and negative selection
# pick from Kosmrj paper
# negative_threshold = 10
# positive_threshold = 2

#matrix_mj = dict or numpy array
#matrix_mj = np.array([ [1,2,3], [4,5,6] ])
#matrix_mj[0,1]

#####################################################################
# LOADING EXPERIMENTAL DATA MODULE
#####################################################################

#######################################################
# Load experimental data
#   generated offline using a bioinformatics pipeline
#   computationally predicts peptides that are likely
#   to be MHC-bound from gene expression data
#######################################################
#PEPTIDE_FILE = os.path.join(DATA_DIR, 
#            'final_list_peptides_combined_abundance_mtec_12028_R41_norm.csv')
#PEPTIDE_DF = pd.read_csv(PEPTIDE_FILE)#, header=None)
#PEPTIDE_DF.columns = ['peptide', 'abundance']
# Initialize bag of peptides and probability of selection
#str_peptide_bag = PEPTIDE_DF.peptide
#f_prob_peptides = PEPTIDE_DF.abundance

# List to store all data (each element of this list is a pandas dataframe)
global list_pd_df_data #= []
# print("Loading synthetic data ... \n")
#
# for temp_file in os.listdir(DATA_DIR):
#     if temp_file.endswith(".csv"):
#         # Load data in pandas dataframe and store in list
#         #print(temp_file)
#
#         # get full path for data file
#         temp_str_peptide_file = os.path.join(DATA_DIR, temp_file)
#
#         # read in pandas dataframe
#         temp_peptide_df = pd.read_csv(temp_str_peptide_file, header=None)
#         #temp_peptide_df.columns = ['peptide', 'abundance']
#
#         #########################################################################
#         # TODO: unit test for this
#         # change length of TCR and peptide here temp_peptide_df.peptide modify to have 3 - 7 aa
#         # possibly modulate for hydrophilic also
#         #temp_peptide_df.peptide = [pep[2:7] for pep in temp_peptide_df.peptide]
#         temp_peptide_df.peptide = [ pep for pep in temp_peptide_df[0] ]
#         #########################################################################
#
#         # store as one member of list
#         list_pd_df_data.append(temp_peptide_df)
#
        
#print(list_pd_df_data[0].peptide)
#print(list_pd_df_data[1].peptide)
        

#####################################################################
# AUXILIARY FUNCTIONS MODULE
#####################################################################

#######################################################
# Other required auxiliary functions
#######################################################

# TODO: optimize cache size later
# @lru_cache(maxsize=LRU_CACHE_SIZE)
def notmatch_two_receptors_within_threshold(tcell_receptor,
                                            peptide,
                                            i_matching_threshold_within_peptide):
    """
    takes two receptors and returns
        True  if they do NOT match within threshold
        False if they        match within threshold
    expects global i_matching_threshold
    """
    # print(sum([int(i==j) for (i,j) in zip(tcell_receptor,peptide)]) < i_matching_threshold)

    ######################
    # DEBUG code follows
    ######################
    # b_temp = sum([int(i==j) for (i,j) in zip(tcell_receptor,peptide)]) < i_matching_threshold
    # if (b_temp == False):
    #    pdb.set_trace()

    # peptide
    # 'KPWLANSKM'
    # tcell_receptor
    # 'KPQLAKHHM'
    # ind = (PEPTIDE_DF.peptide == 'KPWLANSKM')
    # PEPTIDE_DF[ind].head()
    #     peptide  abundance
    # 30427  KPWLANSKM   0.000004

    # print("within function",f_negative_selection_threshold,"\n")
    # print( sum([ matrix_mj[i,j] for (i, j) in zip(tcell_receptor, peptide) ]) )

    # NOTE: all energies negative, so now it becomes all energies must be >
    # return (sum([ matrix_mj[i,j] for (i, j) in zip(tcell_receptor, peptide) ]) > f_negative_selection_threshold)

    # Simple sequential matching model for TESTING
    # TODO: check if needsto call string_array_elementwise_matching_generic
    return (sum([int(i == j) for (i, j) in zip(tcell_receptor, peptide)]) < i_matching_threshold_within_peptide)



def receptor_matching_function_turnstile(tcell_receptor, array_peptide, i_matching_threshold_within_peptide, i_matching_threshold):
    """ Function that takes in T cell receptor and peptide sequence
        and produces matching decision
        True: matched too close (kill)
        False: survive

        :param tcell_receptor: 
        :type tcell_receptor: Str
        :param array_peptide: 
        :type array_peptide: array of Str

    Returns:
        a = ['D', 'I', 'P']
        b = ['A', 'I', 'P']
        receptor_matching_function_turnstile(a, b)
        receptor_matching_function_turnstile(a, [b, b])
        OT-1 receptor
        a = ['G', 'L', 'E', 'Q', 'L'] # 'G', 'L', 'E', 'Q', 'L', 'E'
        Ovalbumin
        b = ['I', 'N', 'F', 'E', 'K'] # 'S', 'I', 'N', 'F', 'E', 'K', 'L'
        from http://www.nature.com/nmeth/journal/v3/n3/fig_tab/nmeth858_F6.html
        receptor_matching_function_turnstile(a, [b, b])
        
        a = ['D', 'I', 'L', 'L', 'L']
        b = ['D', 'I', 'L', 'L', 'L']
        c = ['D', 'I', 'P', 'P', 'P']
        receptor_matching_function_turnstile(a, [b, b], 2)
        2        
        receptor_matching_function_turnstile(a, [b, b, b], 2)
        3
        receptor_matching_function_turnstile(a, [b, b, b], 2)
        4
        receptor_matching_function_turnstile(a, [b, b, b, c], 2)
        3
        
        Returns: bool
        True: matched too close (kill)
        False: survive

    """

    return(
        sum(
            [ not(notmatch_two_receptors_within_threshold(tcell_receptor, peptide, i_matching_threshold_within_peptide))
                for peptide in array_peptide ]
        ) >= i_matching_threshold
    )



def receptor_matching_function_threshold(tcell_receptor, array_peptide, i_matching_threshold_within_peptide):
    """ Function that takes in T cell receptor and peptide sequence
        and produces matching decision
        Requires global i_matching_threshold

        :param tcell_receptor:
        :type tcell_receptor: Str
        :param array_peptide:
        :type array_peptide: array of Str
        :param f_negative_selection_threshold:
        :type f_negative_selection_threshold: Float

    Returns:
        a = ['D', 'I', 'P']
        b = ['A', 'I', 'P']
        receptor_matching_function(a, b, -21.0)
        False
        receptor_matching_function(a, [b, b], -21.0)
        False
        OT-1 receptor
        a = ['G', 'L', 'E', 'Q', 'L'] # 'G', 'L', 'E', 'Q', 'L', 'E'
        Ovalbumin
        b = ['I', 'N', 'F', 'E', 'K'] # 'S', 'I', 'N', 'F', 'E', 'K', 'L'
        from http://www.nature.com/nmeth/journal/v3/n3/fig_tab/nmeth858_F6.html
        receptor_matching_function(a, [b, b], -21.0)
        2
        False

        a = ['D', 'I', 'P', 'L', 'L']
        b = ['A', 'I', 'P', 'L', 'L']
        receptor_matching_function(a, [b,b], -21.0)
        True


        Returns: bool
        True: matched too close (kill)
        False: survive

    """

    # TODO: fix for case when array_peptide has only one peptide
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/29/receptor_matching_function-returns-false
    # TODO: reimplement using map
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/28/new-map-function-to-operate-on-peptide-and
    return (
        not (
            all([notmatch_two_receptors_within_threshold(tcell_receptor, peptide, i_matching_threshold_within_peptide) for
                 peptide in array_peptide])
        )
    )



###########################################################
# Load list of auto sequences from precomputed file
###########################################################

# load from list_all_auto_tcr_offline.csv
#   where pool of auto TCR have been precomputated offline
# auto_tcr_file = pd.read_csv('list_all_auto_tcr_offline.csv', header=None)
# auto_tcr_file = pd.read_csv('list_all_auto_tcr_offline_AIRE_5aa.csv', header=None)
# auto_tcr_file = pd.read_csv('list_all_auto_tcr_offline_AIRE_5aa_1e6.csv', header=None)
#auto_tcr_file = pd.read_csv('list_all_auto_tcr_offline_ALLGENES_5aa_1e6.csv', header=None)


def get_auto_TCR_sequence(i_num_mtecs):
    """ return a TCR sequence which will match some peptide sequence
        within that threshold
        Needs access to global str_peptide_bag

        :param: None
        :type: None

    """

    # TODO: later pick randomly
    # now picks last element of peptide bag (sorted from R script in
    #   descending order); so this is least abundant and hence most
    #   likely to be auto
    # TODO: get_auto_TCR_sequence() will need to select TCR sequence that is less than -ve threshold but has lowest probability
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/35/get_auto_tcr_sequence-will-need-to-select
    # select peptide (TEC) sequence with lowest prob.; then find a complementary sequence that exceeds -ve selection threshold

    # Choose the bag for the first cell
    # choose cell randomly
    # i_temp_cell_index = 0


    # ==============================================================================
    #     i_temp_cell_index = np.random.choice( range(0,i_num_mtecs-1) )
    #     str_peptide_bag = list_pd_df_data[i_temp_cell_index].peptide
    #
    #     # Choose the peptide with the lowest probability from this bag
    #     # CAUTION - assumes that bag is sorted by relative abundance
    #     # TODO: choose from this bag with 1 - abundance probability
    #     #   i.e. preferrably choose low abundance peptides
    #     str_peptide = str_peptide_bag[len(str_peptide_bag) - 1]
    #     print("original peptide", str_peptide)
    #
    #     ##################################################################
    #     # pick something within some threshold/edit distance of peptide
    #     ##################################################################
    #
    #     # initialize TCR sequence
    #     str_TCR_sequence = ''
    #
    #     #while ( (notmatch_two_receptors_within_threshold(tcell_receptor=str_TCR_sequence, peptide=str_peptide) == True):
    #     #    # mutate peptide pointwise
    #     #    i_temp_site = np.random.choice(range(0,len(str_peptide)-1))
    #     #    str_peptide[i_temp_site] = np.random.choice(str_single_aa)
    #
    #     for aa in str_peptide:
    #
    #         # if not exceeded negative selection threshold yet
    #         if ( notmatch_two_receptors_within_threshold(str_TCR_sequence, str_peptide, f_negative_selection_threshold) ):
    #             # mutate pointwise
    #
    #             # first minimum energy of interaction with this aa
    #             min_value = min([ matrix_mj[aa, temp_aa]   for temp_aa in str_single_aa ])
    #
    #             print(min_value)
    #
    #             # do reverse lookup in dict to find that aa which has minimum energy
    #             #[j for (i,j) in matrix_mj[aa, temp_aa].items() if matrix_mj[i,j] == min_value][0]
    #             # find entry in MJ matrix that has the relevant first aa and its interaction with
    #             # another aa has the min. energy; if yes then return that aa
    #             substitution_aa = [j  for ((i,j),k) in matrix_mj.items()
    #                                 if matrix_mj[i,j] == min_value and i == aa ][0]
    #
    #             # Append to string
    #             # CAUTION: this is probably a reference/pointer
    #             print("sub", aa, "for", substitution_aa)
    #             str_TCR_sequence = str_TCR_sequence + substitution_aa
    #
    #         else:
    #             # if matched with self already, then done
    #             # you have auto and now just keep on adding aa
    #             # cannot break since then you will have TCR of length smaller than peptide
    #             print("reached negative selection threshold")
    #             str_TCR_sequence = str_TCR_sequence + aa
    #
    #
    # ==============================================================================

    # #########################################################
    # # TODO: Alternative model of auto TCR generator
    # #       1) generate offline
    # #       2) get a list of all low prob. peptides from all cells
    # array_peptide_lowprob = [ str_peptide_bag[len(str_peptide_bag) - 1] ,
    #                     str_peptide_bag[len(str_peptide_bag) - 2]
    #                     ]
    #
    # # keep on generating random TCR sequences until you get a match
    # while (1==1):
    #
    #     temp_str_TCR_sequence = np.random.choice(str_single_aa, i_length_TCR)
    #
    #     b_decision = receptor_matching_function(temp_str_TCR_sequence, array_peptide_lowprob)
    #     # pdb.set_trace()
    #
    #     if b_decision:
    #         # record # Matched too close and recognizes pathogen
    #         print("Auto sequence. TCR sequence is: ", temp_str_TCR_sequence)
    #         # if matches, then make this the TCR sequence that is returned
    #         str_TCR_sequence = temp_str_TCR_sequence
    #         break
    #     else:
    #         print("Try another auto sequence")
    #
    # #########################################################


    # randomly pick one TCR in this auto pool
    i_temp_num_elements = auto_tcr_file.shape[0]
    i_temp_index = np.random.choice(range(0, i_temp_num_elements))
    temp_tcr = auto_tcr_file.iloc[i_temp_index,]

    # put into an array
    if not b_auto_TCR_sequence_given:
        str_TCR_sequence = [a for a in temp_tcr]
    else:
        # TODO: analyze new offline auto TCR data (AIRE and 5aa)
        #   and get new TCR sequence which is pathologically relevant also
        # TODO: load auto TCR sequence (fixed) from file
        str_TCR_sequence = ['V', 'Q', 'H', 'C', 'F', 'S', 'G', 'C', 'R']

    # final negative selection status of putative auto TCR sequence
    # TODO: if not matched within threshold, repeat same procedure

    # print("Does final TCR sequence match peptide",
    #      not (notmatch_two_receptors_within_threshold(str_TCR_sequence, str_peptide, f_negative_selection_threshold) )
    #      )

    # TODO: also check if positive_selection()
    #print("Does this auto TCR sequence positively match at least one self-peptide?",
    #      positive_selection(str_TCR_sequence, i_num_mtecs)
    #      )

    print("modified auto TCR: ", str_TCR_sequence)
    return (str_TCR_sequence)


def get_complementary_to_TCR_sequence(str_survived_TCR_sequence, f_negative_selection_threshold):
    """ Given a TCR sequence (which has for example survived thymic selection),
        generate one or a set of peptides which are complementary to this
        i.e. a sequence of peptides which will match this TCR beyond the threshold

        :param str_survived_TCR_sequence: Str
        :type str_survived_TCR_sequence: Str
    """

    # TODO: later pick randomly
    # now picks last element of peptide bag (sorted from R script in
    #   descending order); so this is least abundant and hence most
    #   likely to be auto
    # TODO: get_auto_TCR_sequence() will need to select TCR sequence that is less than -ve threshold but has lowest probability
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/35/get_auto_tcr_sequence-will-need-to-select
    # select peptide (TEC) sequence with lowest prob.; then find a complementary sequence that exceeds -ve selection threshold

    str_survived_TCR_sequence
    print("original TCR sequence that has survived negative selection", str_survived_TCR_sequence)

    ##################################################################
    # pick something within some threshold/edit distance of peptide
    ##################################################################

    # initialize TCR sequence
    str_peptide_sequence = ''

    # while ( (notmatch_two_receptors_within_threshold(tcell_receptor=str_TCR_sequence, peptide=str_peptide) == True):
    #    # mutate peptide pointwise
    #    i_temp_site = np.random.choice(range(0,len(str_peptide)-1))
    #    str_peptide[i_temp_site] = np.random.choice(str_single_aa)

    for aa in str_survived_TCR_sequence:

        # if not exceeded negative selection threshold yet
        if (notmatch_two_receptors_within_threshold(str_peptide_sequence, str_survived_TCR_sequence,
                                                    f_negative_selection_threshold)):
            # mutate pointwise until you exceed threshold

            # first minimum energy of interaction with this aa
            min_value = min([matrix_mj[aa, temp_aa] for temp_aa in str_single_aa])

            print(min_value)

            # do reverse lookup in dict to find that aa which has minimum energy
            # [j for (i,j) in matrix_mj[aa, temp_aa].items() if matrix_mj[i,j] == min_value][0]
            # find entry in MJ matrix that has the relevant first aa and its interaction with
            # another aa has the min. energy; if yes then return that aa
            substitution_aa = [j for ((i, j), k) in matrix_mj.items()
                               if matrix_mj[i, j] == min_value and i == aa][0]

            # Append to string
            # CAUTION: this is probably a reference/pointer
            print("sub", aa, "with", substitution_aa)
            str_peptide_sequence = str_peptide_sequence + substitution_aa

        else:
            # if matched with self already, then done
            # you have auto and now just keep on adding aa
            # cannot break since then you will have complementary of length smaller than peptide
            print("reached negative selection threshold")
            str_peptide_sequence = str_peptide_sequence + aa

    # final negative selection status of putative auto TCR sequence
    print("Does final TCR sequence match peptide",
          not (notmatch_two_receptors_within_threshold(str_peptide_sequence, str_survived_TCR_sequence,
                                                       f_negative_selection_threshold))
          )

    # TODO: also check if positive_selection()
    print("complementary peptide:", str_peptide_sequence)
    return (str_peptide_sequence)


def positive_selection(tcell_receptor, i_num_mtecs):
    """ Models positive selection in thymus cortex
            Given a TCR sequence, it checks if it positively matches a peptide

    :param tcell_receptor:
    :type tcell_receptor: array of Str


    Returns:

    """

    # TODO: get mouse MHC 1 genes H2kb etc and predict/get MHC 1 proteins
    # get protein from UniProt

    # randomly choose cell
    i_temp_cell_index = np.random.choice(range(0, i_num_mtecs - 1))
    list_all_peptides_inthat_cell = list_pd_df_data[i_temp_cell_index].peptide

    # go through all peptides in that cell until one that matches above
    #   positive selection threshold
    for peptide in list_all_peptides_inthat_cell:
        # any match then break out and return True
        if (sum([matrix_mj[i, j] for (i, j) in zip(tcell_receptor, peptide)]) < f_positive_selection_threshold):
            print("TCR positively matched in cell index: ", i_temp_cell_index,
                  " with peptide", peptide)
            return True

    # if no match, return False
    print("TCR does not match positive")
    return False


def multiple_from_bag_of_peptides(i_tec_id, b_generate_mtec_random_aa, i_num_peptide_present):
    """ Return multiple strings of peptides for mTEC from real data
        (called only during initialization)
        Needs the unique id for an mTEC to index into the unique bag of
        peptides for each mTEC.
        Requires global i_num_peptide_present and list_pd_df_data

        :param i_tec_id: ID of TEC that has been selected
        :type i_tec_id: Int
        :param b_generate_mtec_random_aa: If False, generate peptides from
                amino acids with all same frequency, If True then from
                mice proteome frequencies
        :type b_generate_mtec_random_aa: Bool
        :param i_num_peptide_present: number of synapses
        :type i_num_peptide_present: Int

        :rtype: String
    """

    if not (b_generate_mtec_random_aa):
        # if not random idd, then use real data
        # Get peptides list and abundance for the unique id or unique bag for this mTEC
        #   index into global list of pandas dataframe
    
        #str_peptide_bag = list_pd_df_data[i_tec_id].peptide
        #f_prob_peptides = list_pd_df_data[i_tec_id].abundance
        str_peptide_bag = list_pd_df_data[i_tec_id]

        # Draw from this bag of peptides according to probabilities (normalized abundances) of peptides
        #str_array_peptide = np.random.choice(str_peptide_bag, i_num_peptide_present, p=f_prob_peptides)
        str_array_peptide = np.random.choice(str_peptide_bag[0], i_num_peptide_present)
    else:
        # randomly iid at each site from mice proteome frequencies
        str_array_peptide = [np.random.choice(str_single_aa, i_length_peptide, p=f_aa_mice_proteome_frequencies)
                             for _ in range(1, i_num_peptide_present)]

    return (str_array_peptide)


def return_TCR_survives_mult_rounds():
    """ Analyze data about TCRs that survive multiple rounds
            of thymic selection
    """

    str_survived_tcr = ''

    return str_survived_tcr


def analyze_aa_frequencies(list_tcell_sequence):
    """
    """

    # call return_TCR_survives_mult_rounds() repeatedly
    # store in     array_str_survived_tcr
    # analyze site specific probablities/frequencies
    # analyze non-site specific probablities

    # [for tcr in list_tcell_sequence]
    import scipy.stats
    return (scipy.stats.itemfreq(list_tcell_sequence))


def challenge_survivetcellpool_pathogen(list_tcr_escape_thymus,
                                        str_misc_output_file_name,
                                        i_matching_threshold_within_peptide,
                                        i_matching_threshold,
                                        i_num_peptide_present,
                                        str_turnstile_threshold):
    """ Given a TCR pool that has survived thymic selection,
            challenge this with a pathogen
    """

    # TODO: advanced challenge experiment (take pathogen/virus sequence -> protein -> NetMHC) and challengt to surving Tcells
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/51/advanced-challenge-experiment-take
    #   get sequences from virus/pathogen databases
    #   https://www.fludb.org/brc/sssearch.spg?method=ShowCleanInputPage&decorator=corona&preSelectDB=true
    #   https://www.fludb.org/brc/proteinSequence.spg?proteinGi=225403261&decorator=corona&range=2-10&context=SS_760020346912
    #   https://www.fludb.org/brc/viprDetails.spg?ncbiProteinId=ACN89743&decorator=corona


    generate_pathogen_peptide_pool = \
        [np.random.choice(str_single_aa, i_length_peptide, p=f_pathogen_frequencies)
         for _ in range(1, i_num_peptide_present)]

    # pdb.set_trace()
    # check if generate_pathogen_peptide_pool is array of peptides

    # open misc filename handles (NOTE: append mode here)
    # CAUTION: assumes that file already opened
    with open(str_misc_output_file_name, 'a', newline='\n') as outfile:
        wr = csv.writer(outfile)  # , quoting=csv.QUOTE_ALL)
        item = "Challenge negatively selected TCR pool with pathogen"
        wr.writerow(item)

        # for each TCR that has escaped, test if it matches with pathogen
        for tcr in list_tcr_escape_thymus:
            if str_turnstile_threshold == 'turnstile':
                b_decision = receptor_matching_function_turnstile(tcr, generate_pathogen_peptide_pool,
                                                                  i_matching_threshold_within_peptide, i_matching_threshold)
            elif str_turnstile_threshold == 'threshold':
                b_decision = receptor_matching_function_threshold(tcr, generate_pathogen_peptide_pool,
                                                                  i_matching_threshold_within_peptide)

            # pdb.set_trace()

            if b_decision:
                # record # Matched too close and recognizes pathogen
                print("Selected TCR recognizes pathogen. TCR sequence is: ", tcr)

                # save all output to misc file
                item = "Selected TCR recognizes pathogen. TCR sequence is: "
                wr.writerow(item)
                item = tcr
                wr.writerow(item)
            else:
                # record
                pass

    outfile.close()

    return (generate_pathogen_peptide_pool)


def degeneracy_pathogen_against_survivetcellpool_pathogen(list_generate_pathogen_peptide_pool,
                                                        list_tcr_escape_thymus,
                                                        str_misc_output_file_name,
                                                        i_matching_threshold_within_peptide,
                                                        i_matching_threshold,
                                                        i_num_peptide_present,
                                                        str_turnstile_threshold,
                                                        str_degen_diseasepep_againstTCR_filename):
    """ Given a pathogen peptide pool that has been generated in challenge_survivetcellpool_pathogen()
            test the degeneracy of each non-self peptide against the TCR pool that has survived thymic selection,
    """

    # pdb.set_trace()
    # check if generate_pathogen_peptide_pool is array of peptides

    # open misc filename handles (NOTE: append mode here),
    # CAUTION: assumes that file already opened

    list_degeneracy_nonself_against_TCR_repertoire = []

    with open(str_misc_output_file_name, 'a', newline='\n') as outfile:
        wr = csv.writer(outfile)  # , quoting=csv.QUOTE_ALL)
        item = "Degeneracy of non-self against TCR repertoire"
        wr.writerow(item)

        # for each non-self, test it against all TCRs
        for temp_disease_peptide in list_generate_pathogen_peptide_pool:


            #########################################################################
            # check additional how many degenerate interactions does this TCR have?
            #########################################################################

            i_degeneracy_escape = count_degeneracy_function(tcell_receptor=temp_disease_peptide,
                                                            array_peptide=list_tcr_escape_thymus,
                                                            i_matching_threshold_within_peptide=i_matching_threshold_within_peptide)

            # if i_degeneracy_escape < 3:
            #    print("Printing low degeneracy TCR sequence: ", tcell.receptor)
            #    print("Degeneracy is: ", i_degeneracy_escape)

            # append to list of degeneracy ALL INCOMING RANDOM THYMOCYTES
            list_degeneracy_nonself_against_TCR_repertoire.append(str(i_degeneracy_escape))


            # pdb.set_trace()

        # Quantify fraction detected
        item = "Detection fraction:"
        wr.writerow(item)

        f_detect_frac = sum([int(x) > 0 for x in list_degeneracy_nonself_against_TCR_repertoire]) / len(
            list_degeneracy_nonself_against_TCR_repertoire)
        item = str(f_detect_frac)
        wr.writerow(item)


    # close file handle
    outfile.close()

    # save list of disease peptides that are recognized to file
    # only if there is a list of disease peptides that are recognized
    if not (list_degeneracy_nonself_against_TCR_repertoire == []):
        fn_save_list_to_file(list_input=list_degeneracy_nonself_against_TCR_repertoire,
                             str_filename=str_degen_diseasepep_againstTCR_filename)



def challenge_survivetcellpool_diseasepeptide(list_tcr_escape_thymus,
                                              str_misc_output_file_name,
                                              i_matching_threshold_within_peptide,
                                              i_matching_threshold,
                                              str_turnstile_threshold):
    """ Given a TCR pool that has survived thymic selection,
            challenge this with a disease causing peptide
            to observe if this Tcell pool may react with these peptides in the periphery
    """

    # TODO: advanced challenge experiment (take peptide sequence -> protein -> NetMHC) and challenge to surving Tcells
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/73/run-for-21-hydroxylse-and-ins1-and-ins2-in
    #   get sequences from peptide databases
    #   http://immunet.cn/bdb/index.php
    # also from Toy transcript data (21 hrdroxylse and Ins1/2)

    # generate_disease_peptide_pool = \
    #    [np.random.choice(str_single_aa, i_length_peptide, p=f_pathogen_frequencies)
    #     for _ in range(1, i_num_peptide_present)]
    # generate_disease_peptide_pool = ['V', 'A', 'R', 'V', 'W'] #['V', 'P', 'V', 'A', 'R', 'V', 'W', 'L']
    # the following peptide is from transcript for Ins1 (ENSMUST00000036952) run through NETMHC4.0 (strong binder)

    # generate_disease_peptide_pool = [['H', 'L', 'V', 'E', 'A'], ['H', 'L', 'V', 'E', 'A']]#['G','P','H','L','V','E','A','L']

    # the following peptide is from transcript for Cyp21a1 (ENSMUST00000025223) run through NETMHC4.0 (strong binder)
    # generate_disease_peptide_pool = ['H', 'L', 'V', 'E', 'A']#['G','P','H','L','V','E','A','L']

    # peptides for Cyp21a1 in toy data ENSMUST00000025223
    # gives protein http://www.uniprot.org/uniprot/A0A0R4J048.fasta
    # which gives the following strong binders from NetMHC 4.0
    generate_disease_peptide_pool = [['G', 'I', 'N', 'L', 'P'], ['P', 'G', 'S', 'Q', 'L'],
                                     ['P', 'R', 'T', 'P', 'S'], ['W', 'A', 'V', 'A', 'F'],
                                     ['P', 'L', 'L', 'R', 'F'], ['T', 'T', 'A', 'T', 'T'],
                                     ['L', 'P', 'S', 'K', 'F'], ['Y', 'A', 'G', 'I', 'N'],
                                     ['N', 'L', 'P', 'I', 'Y'], ['L', 'A', 'P', 'G', 'F'],
                                     ['Q', 'K', 'W', 'V', 'D'], ['I', 'P', 'P', 'F', 'Q'],
                                     ['P', 'L', 'A', 'P', 'G']
                                    ]

    # randomly iid at each site from mice proteome frequencies
    #generate_disease_peptide_pool = [np.random.choice(str_single_aa, i_length_peptide, p=f_aa_mice_proteome_frequencies)
    #                                    for _ in range(1, i_num_peptide_present)]
    # [['Y','A','G','I','N','L','P','I']
    # LGPGSQLL
    # KNPRTPSF
    # LSWAVAFL
    # IIPLLRFL
    # TETTATTL
    # WELPSKFW
    # QPYAGINL
    # QPNLPIYL
    # PPLAPGFL
    # LIQKWVDF
    # LPIPPFQV
    # LPPLAPGF

    # pdb.set_trace()
    # check if generate_disease_peptide_pool is array of peptides

    # open misc filename handles (NOTE: append mode here)
    # CAUTION: assumes that file already opened
    with open(str_misc_output_file_name, 'a', newline='\n') as outfile:
        wr = csv.writer(outfile)  # , quoting=csv.QUOTE_ALL)
        item = "Challenge negatively selected TCR pool with disease peptide"
        wr.writerow(item)

        # for each TCR that has escaped, test if it matches with pathogen
        for tcr in list_tcr_escape_thymus:
            if str_turnstile_threshold == 'turnstile':
                b_decision = receptor_matching_function_turnstile(tcr, generate_disease_peptide_pool,
                                                                  i_matching_threshold_within_peptide, i_matching_threshold)
            elif str_turnstile_threshold == 'threshold':
                b_decision = receptor_matching_function_threshold(tcr, generate_disease_peptide_pool,
                                                                  i_matching_threshold_within_peptide)
            # pdb.set_trace()

            if b_decision:
                # record # Matched too close and recognizes pathogen
                print("Selected TCR recognizes disease peptides. TCR sequence is: ", tcr)

                # save all output to misc file
                item = "Selected TCR recognizes disease peptides. TCR sequence is: "
                wr.writerow(item)
                item = tcr
                wr.writerow(item)
            else:
                # record
                pass

    outfile.close()



def count_degeneracy_function(tcell_receptor, array_peptide, i_matching_threshold_within_peptide):
    """ Function that takes in T cell receptor and peptide sequence
        and produces sum of number of degenerate interactions
        Requires global i_matching_threshold

    Args:
        peptide:

    Returns:
        a = ['D', 'I', 'P']
        b = ['A', 'I', 'P']
        count_degeneracy_function(a,b, -21.0)
        2
        count_degeneracy_function(a,[b,b], -21.0)
        count_degeneracy_function(a,[b,b,b], -21.0)
        
        a = ['C', 'M', 'C', 'M']
        b = ['C', 'M', 'C', 'M']
        [  (notmatch_two_receptors_within_threshold(a, peptide)) for peptide in [b,b] ]
        [False, False]
        [  int(notmatch_two_receptors_within_threshold(a, peptide)) for peptide in [b,b] ]
        [0, 0]
        count_degeneracy_function(a, [b,b], -21.0)
        2

        Returns: bool
        True: matched too close (kill)
        False: survive

    """

    return ( sum([  int(not(notmatch_two_receptors_within_threshold(tcell_receptor, peptide, i_matching_threshold_within_peptide))) for peptide in array_peptide])
    )


#####################################################################
# CLASS MODULE
#####################################################################

#######################################################
# Class definitions of all agents
#######################################################



class ImmuneCellAgent():
    """A T-cell agent with fixed initial receptor."""

    def __init__(self, unique_id, b_generate_tcr_random_aa, b_auto_flag, i_num_mtecs, list_tcr_sequence_fixed):
        """ Constructor for T-cell/thymocyte

        Args:
            unique_id:

        Returns:

        """
        self.i_cell_type = 1  # T-cell
        self.unique_id = unique_id

        # If auto flag True, then get auto TCR sequence
        if b_auto_flag:
            self.receptor = get_auto_TCR_sequence()
            # else get normal random sequence
        else:
            if list_tcr_sequence_fixed == 'notcrseq':
                # if no fixed TCR sequence supplied, then generate one
                self.receptor = self.bag_of_TCR(b_generate_tcr_random_aa, i_num_mtecs)
            else:
                # if fixed TCR sequence supplied, then use it
                self.receptor = list_tcr_sequence_fixed
                print(list_tcr_sequence_fixed)

    def bag_of_TCR(self, b_generate_tcr_random_aa, i_num_mtecs):
        """ Return one string of TCR
        pick randomly given 20^N combinations
        (called only during initialization)

        :return: string
        """

        # TODO: get better logic here on T-cell receptor repertoire
        # TODO: testing are these equal probabilities?
        # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/14/unit-testing-for-bag_of_tcr
        # TODO: run ABM with bag_of_TCR() with freuquency from human proteome
        #   refactor bag_of_TCR() to take two options (no frequency/uniform or from human proteome)
        #   make global function that both Tcell and TEC can call
        #   Null model: compare gene sequencing model to frequency from human proteome model
        # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/32/run-abm-with-bag_of_tcr-with-freuquency

        if (not (b_generate_tcr_random_aa)):
            # generate TCR sequences completely at random (uniform)
            str_TCR_sequence = np.random.choice(str_single_aa, i_length_TCR)  # , p=[0.5, 0.1, 0.1, 0.3])
        else:
            # generate TCR sequences at random but with probability of aa following mice proteome frequencies
            str_TCR_sequence = np.random.choice(str_single_aa, i_length_TCR,
                                                p=f_aa_mice_proteome_frequencies)

        # str_TCR_sequence = 'YANDSVRKK'

        # TODO: also check if positive_selection()
        # while(1 == 1):
        #    str_TCR_sequence = np.random.choice(str_single_aa, i_length_TCR)
        #    if (positive_selection(str_TCR_sequence)):
        #        break
        #print("Does this random TCR sequence positively match at least one self-peptide?",
        #      positive_selection(str_TCR_sequence, i_num_mtecs)
        #      )

        print("TCR sequence is: ", str_TCR_sequence, "\n")

        return (str_TCR_sequence)

    def sensing(self, model):
        """ Sensing module for T-cell

        Returns:

        """

        # pass

    def apoptosis(self):
        """ Death function for T-cell/thymocyte
        """

        # pass

    def natural_death(self, model):
        """ Natural death of thymocytes or exit from compartment

        :return:
        """

        # if np.random.uniform(0,1) < tcell_death_probability:
        #        # self.apoptosis(model)
        #        pass


class EpithelialCellAgent():
    """A thymic epithelial cell agent with fixed initial receptor."""

    def __init__(self, unique_id, b_generate_mtec_random_aa, i_num_peptide_present):
        """ Constructor for TEC

        Args:
            unique_id:

        Returns:

        """

        self.i_cell_type = 2  # Thymic Epithelial-cell
        self.unique_id = unique_id  # a unique id for each mTEC
        # self.receptor = self.bag_of_peptides()
        # self.receptor = self.multiple_from_bag_of_peptides()
        # Redraw from bag of peptides for this mTEC (given by unique id) and reset presented peptides
        self.redraw_from_bag(self.unique_id, b_generate_mtec_random_aa, i_num_peptide_present)

    def redraw_from_bag(self, i_tec_id, b_generate_mtec_random_aa, i_num_peptide_present):
        """ Redraw from bag of peptides for this mTEC (given by unique id)
                and reset presented peptides
            Primary mTEC receptor setting function
        """

        self.receptor = multiple_from_bag_of_peptides(i_tec_id, b_generate_mtec_random_aa, i_num_peptide_present)

    def natural_death(self, model):
        """ Natural death or exit from compartment of TECs

        :return:
        """

        # if np.random.uniform(0,1) < mtec_death_probability:
        # inspired by https://github.com/projectmesa/mesa-examples.git
        # ==============================================================================
        #             x, y = self.pos
        #             this_cell = model.grid[y][x]
        #             # find ImmuneCellAgent
        #             epithelialcell = [obj for obj in this_cell if isinstance(obj, EpithelialCellAgent)]
        #             # Kill the immunecell
        #             model.grid[y][x].remove(epithelialcell)
        #             model.schedule.remove(epithelialcell)
        # ==============================================================================
        #   pass


class CompartmentThymus():
    """ Compartment for complete simulation (thymus)
    """

    def __init(self):
        """ Constructor for thymus compartment
        """
        pass

    def step(self, model):
        """ Function to implement logic for each step
        """
        pass


class BiologicalModel():
    """A biological model with some number of agents.
        Contains both T-cells and mTECs.
    """

    def __init__(self, N, M):
        """ Constructor for model. Initialize simulation

        Args:
            N: number of T-cells
            M: number of mTECs

        Returns:

        """
        ##########################################
        # Initialize
        ##########################################
        self.running = True
        self.num_agents_tcells = N
        self.num_agents_epithelial = M
        self.num_tcells_remaining = N
        self.num_mtecs_remaining = M
        self.num_tcells_matched = 0

        ##########################################
        # Create agents and add them to grid
        ##########################################
        # TODO: Use multiple cores using joblib and multiprocessing
        # for i in range(self.num_agents_tcells):
        #    a = ImmuneCellAgent(i)

        # if b_init_mtec:
        # TODO: Use multiple cores using joblib and multiprocessing
        # for i in range(self.num_agents_epithelial):
        #    a = EpithelialCellAgent(i)

    def step(self):
        """ Step logic for model

        Returns: Advance the model one step
        """

        # TODO: speed this up
        # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/17/speed-up-get_total_tcells
        # Also may not need pandas dataframe; use just a list/dict
        # NOTE CAUTION this assumes that mTECs do not die or exit compartment

        # list_tcell_count.append( self.schedule.get_agent_count() - self.num_agents_epithelial )

    def birth(self):
        """ Birth of new T-cells during the simulation (entry into compartment)

        Returns:

        """
        # TODO: implement birth of new T-cells during the simulation (entry into compartment)

        # pass



if __name__ == "__main__":

    ####################################################################
    # Get command line arguments
    ####################################################################
    # TODO: use argparse like sample_path = args.get('--sample')
    # TODO: make generic output directory that can be used to run on different machines
    # TODO: append to other filenames
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/72/make-generic-output-directory-that-can-be
    #str_output_directory = "/scratch/banerjees/agent_based_model_mesa/gambit"
    str_output_filename_survival = str(sys.argv[1]) # 'tcell_survival.csv'
    str_output_filename_matched  = str(sys.argv[2]) # 'tcell_matches.csv'
    i_num_independent_tcells     = int(sys.argv[3]) # number of independent T-cells to be generated and tested
    #b_generate_tcr_random_aa     = bool(sys.argv[4]) # If False, generate TCR sequences completely at random (uniform)
        # If True, generate TCR sequences at random but with probability of aa following mice proteome frequencies
    
    if   (sys.argv[4] == "False"):
        b_generate_tcr_random_aa = False
    elif (sys.argv[4] == "True"):
        b_generate_tcr_random_aa = True
        
    #b_generate_mtec_random_aa    = bool(sys.argv[5]) # If False, generate mTEC peptide sequences completely at random (uniform)
        # If True, generate mTEC peptide sequences at random but with probability of aa following mice proteome frequencies
    if   (sys.argv[5] == "False"):
        b_generate_mtec_random_aa = False
    elif (sys.argv[5] == "True"):
        b_generate_mtec_random_aa = True
        
    i_num_interactions = int(sys.argv[6]) # number of interactions a single T-cell has with mTECs
        # during its residence in the thymus
    
    #b_auto_flag if True, then generate T-cells with auto sequence
    # If False else generate TCR with random sequences
    if   (sys.argv[7] == "False"):
        b_auto_flag = False
    elif (sys.argv[7] == "True"):
        b_auto_flag = True

    # negative selection threshold
    i_matching_threshold_within_peptide = int(sys.argv[8])
    # print(f_negative_selection_threshold)

    i_num_peptide_present = int(sys.argv[9])
    print(i_num_peptide_present)

    i_num_mtecs = int(sys.argv[10])
    print(i_num_mtecs)

    str_checkauto = str(sys.argv[11])

    str_turnstile_threshold = str(sys.argv[12])

    str_tcr_sequence_fixed_FLAG = (sys.argv[13]) #",".join(list(sys.argv[11]))
    if (str_tcr_sequence_fixed_FLAG == 'notcrseq'):
        # NO TCR sequence supplied
        list_tcr_sequence_fixed = 'notcrseq'
    else:
        # TCR sequence supplied so parse it
        list_tcr_sequence_fixed = [a for a in str_tcr_sequence_fixed_FLAG]

    print(list_tcr_sequence_fixed)


    print("Loading synthetic data ... \n")

    # pick directory based on number of mTECs
    #if i_num_mtecs == 174:
    #    DATA_DIR = "synthetic_data"  # synthetic data generated

    #if i_num_mtecs == 10:
    #    DATA_DIR = "synthetic_data/synthetic_data_10mtec"  # synthetic data generated

    if i_num_mtecs == 100:
        DATA_DIR = "data"  # synthetic data generated
		
    print(i_num_mtecs)	

    #if i_num_mtecs == 1000:
    #    DATA_DIR = "synthetic_data/synthetic_data_1000mtec"  # synthetic data generated

    global list_pd_df_data  # = []

    list_pd_df_data = []

    i_temp_counter_formtecs = 0 	
    for temp_file in os.listdir(DATA_DIR):
        if temp_file.endswith(".CSV"):
            # Load data in pandas dataframe and store in list
            # print(temp_file)

            # get full path for data file
            temp_str_peptide_file = os.path.join(DATA_DIR, temp_file)

            # read in pandas dataframe
            temp_peptide_df = pd.read_csv(temp_str_peptide_file, header=None)
            # temp_peptide_df.columns = ['peptide', 'abundance']

            #########################################################################
            # TODO: unit test for this
            # change length of TCR and peptide here temp_peptide_df.peptide modify to have 3 - 7 aa
            # possibly modulate for hydrophilic also
            # temp_peptide_df.peptide = [pep[2:7] for pep in temp_peptide_df.peptide]
            temp_peptide_df.peptide = [pep for pep in temp_peptide_df[0]]
            #########################################################################

            i_temp_counter_formtecs = i_temp_counter_formtecs + 1			
            # store as one member of list
            list_pd_df_data.append(temp_peptide_df)


    print(i_temp_counter_formtecs)
    # pdb.set_trace()	
	
    # Create GIANT LIST of all peptides        
    print("Generating GIANT CELL with all synthetic peptides ... \n")        
    # [x for x in list_pd_df_data[0].iloc[0:,0]  ]        
    # [x for i_count in range(0,174)   for x in list_pd_df_data[i_count].iloc[0:,0]  ]
    # [x for i_count in range(0,174)   for x in list_pd_df_data[i_count].iloc[0:,0]  ]       
    list_giant_ALL_COMBINED_peptides = [x for i_count in range(0,i_num_mtecs) 
                                          for x in list_pd_df_data[i_count].iloc[0:,0]  ]


    # arg_parser = argparse.ArgumentParser(description="argument_parser")
    # arg_parser.add_argument("--neg_seln_thresh",
    #                        help="negative selection threshold (Float)",
    #                        type=float)
    # arguments = arg_parser.parse_args()
    # print(arguments.neg_seln_thresh)

    # str_input_param_filename = str(sys.argv[9]) # model_input_parameters.py
    # from eval(model_input_parameters) import *
    # pdb.set_trace() # debug command

    ####################################################################
    # Print parameters of run
    ####################################################################

    print("\nParameters of run")
    print("*****************************************************")
    print("Number of peptides presented on each mTEC: ", i_num_peptide_present)
    print("Auto T-cell: ", b_auto_flag, "\n", 
          "Negative selection Threshold: ", i_matching_threshold_within_peptide, "\n",
          "TCR generated from proteome freq.?: ", b_generate_tcr_random_aa, "\n",
          "mTEC generated from proteome freq.?:", b_generate_mtec_random_aa, "\n",
          "Number of independent T-cells generated:", i_num_independent_tcells, "\n",
          "Number of interactions a single T-cell has with mTECs:", i_num_interactions, "\n",
          "Negative selection energy threshold", i_matching_threshold_within_peptide, "\n",
          "Number of immunological synapses", i_num_peptide_present, "\n",
          "Number of peptides that must match on a single interaction (TCR-pMHC) to negatively select:", i_matching_threshold, "\n",
          "Number of thymic epithelial cells", i_num_mtecs, "\n"
          )
    print("Number of Monte-Carlo runs:", i_num_montecarlo)          
    print("*****************************************************\n")


    #####################################################################
    # STOCHASTIC SIMULATOR MODULE
    #####################################################################

    for _ in range(1,i_num_independent_tcells + 1):
        ####################################################################
        # Create model and initialize all agents
        ####################################################################
    
        print("Initializing stochastic simulator ... \n")
        first_model = BiologicalModel(i_num_tcells, i_num_mtecs) # b_init_mtec
    
        # Create a T-cell/thymocyte
        tcell = ImmuneCellAgent(unique_id = 1,
                                b_generate_tcr_random_aa = b_generate_tcr_random_aa, 
                                b_auto_flag = b_auto_flag,
                                i_num_mtecs=i_num_mtecs,
                                list_tcr_sequence_fixed=list_tcr_sequence_fixed)

        ####################################################################
        # Initialize variables to hold metrics about simulation
        ####################################################################
    
        # Initialize list with the number of Tcells remaining
        list_tcell_count.append(first_model.num_tcells_remaining)
        # Initialize list with number of Tcells recognized    
        list_tcell_matches.append(first_model.num_tcells_matched)
        # Initialize list of mTECs that have been initialized
        list_index_TEC_initialized = []
    
        ####################################################################
        # Run simulation
        ####################################################################
    
        # do Monte-Carlo runs
        for i_temp_counter_montecarlo in range(1,i_num_montecarlo + 1):
            
            #start = time.clock()
        
            # TODO: can potentially be made parallel (joblib or multiprocessing)
            for _ in range(1,i_num_interactions + 1):
                
                #first_model.step() # step function for model
                
                # Choose an mTEC randomly (i_random_index goes from 1 to number of mTECs)
                i_random_index = np.random.choice(range(1,i_num_mtecs))
                
                # Initialize that mTEC with a new set of presented peptides
                mtec = EpithelialCellAgent(unique_id = i_random_index, 
                                           b_generate_mtec_random_aa = b_generate_mtec_random_aa,
                                           i_num_peptide_present=i_num_peptide_present)
                
                #if i_random_index not in list_index_TEC_initialized:
                    # if this particular mTEC never called before, then initialize
                    #mtec = EpithelialCellAgent(i_random_index)
                    # add index to list of mTECs that have been initialized 
                    #list_index_TEC_initialized.append(i_random_index)
                #else:
                    # if this mTEC already exists, then just call the member function
                    # to redraw from unique bag of peptides and reset receptor
                    # object.redraw_from_bag() need object given i_random_index
                    
                
                # explictly call Tcell-mTEC interaction
                if str_turnstile_threshold == 'turnstile':
                    b_decision = receptor_matching_function_turnstile(tcell_receptor=tcell.receptor,
                                                                  array_peptide=mtec.receptor,
                                                                  i_matching_threshold_within_peptide=i_matching_threshold_within_peptide,
                                                                  i_matching_threshold=i_matching_threshold)
                elif str_turnstile_threshold == 'threshold':
                    b_decision = receptor_matching_function_threshold(tcell_receptor=tcell.receptor,
                                                                      array_peptide=mtec.receptor,
                                                                      i_matching_threshold_within_peptide=i_matching_threshold_within_peptide)

                if b_decision:
                    # Matched too close to mTEC (self), hence needs to be killed
                
                    #==============================================================================
                    #                     print("TCR matches with mTEC peptides within specified threshold. Parameters are \n") 
                    #                     print("Matched with mTEC cell index:", i_random_index)
                    #                     print("Monte-Carlo run number: ", i_temp_counter_montecarlo)
                    #                     
                    #==============================================================================
                    #print("TCR: ", tcell.receptor)
                    #print("mTEC peptides: ", mtec.receptor)
        
                    # update number of T-cells remaining
                    # decrement counter only if > 1
                    if (first_model.num_tcells_remaining >= 1):
                        first_model.num_tcells_remaining = first_model.num_tcells_remaining - 1
                                    
                    # update model variable with number of Tcells recognized
                    first_model.num_tcells_matched = first_model.num_tcells_matched + 1

                    # save TCR sequence that is killed before breaking
                    list_all_tcell_sequence_killed.append(tcell.receptor)

                    # DEBUGGING command
                    #pdb.set_trace()
                    
    
        
                # update list with the number of Tcells remaining
                list_tcell_count.append(first_model.num_tcells_remaining)
                
                # update list with number of Tcells recognized    
                list_tcell_matches.append(first_model.num_tcells_matched)
                
                
                if b_decision:
                    # Matched too close to mTEC (self), hence needs to be killed
                    # this particular T cell killed; hence no need to have more interactions
                    # break out of interaction loop
                    break

        
            # end for loop for number of interactions
        
            #end = time.clock()
        
            ####################################################################
            # Diagnostics
            ####################################################################
        
            #print("model execution time: ", end - start, "seconds")
            #print(first_model.num_agents_tcells) # public member
            #print(first_model.running)
            
        # end Monte-Carlo for loop    
    
        ####################################################################
        # Plotting
        ####################################################################
    
        # TODO: change this to new data structure so as to accommodate epithelial_receptor
        # TODO: now receptor is a char array, change
        # list_receptor_dist = [a.receptor for a in first_model.schedule.agents]
        #plt.hist(list_receptor_dist)#, bins=30)
        #plt.title("Distribution of Tcells")
        #plt.xlabel("Receptors")
        #plt.show()
        #plt.savefig("hist_receptor_distribution.png")
    
    
        ####################################################################
        # Save data to disk as pandas dataframe or pickle or csv
        ####################################################################
        # Save to csv file
        #df_tcell_survive.to_csv('tcell_survival.csv')
    
        #==============================================================================
        #         with open(str_output_filename_survival + '.csv', 'w', newline='\n') as outfile:
        #             wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        #             wr.writerow(list_tcell_count)
        #             
        #         with open(str_output_filename_matched + '.csv', 'w', newline='\n') as outfile:
        #             wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        #             wr.writerow(list_tcell_matches)
        #             
        #==============================================================================
    
        ####################################################################
        # Summary statistics of model output
        # Analyze model output
        ####################################################################
    
        # df_tcell_survive = first_model.datacollector.get_model_vars_dataframe()
        # df_tcell_survive.plot()
        # #plt.show()
    
        #==============================================================================
        #         plt.figure()
        #         plt.plot( (list_tcell_count) )
        #         plt.title("T-cell survival of negative selection in thymus")
        #         plt.xlabel('Number of interactions')
        #         plt.ylabel('Number of T-cells')
        #         plt.legend(['Surviving T-cells'])
        #         plt.savefig(str_output_filename_survival + ".png")
        #         plt.close()
        #         
        #         plt.figure()
        #         plt.plot( (list_tcell_matches) )
        #         plt.title("TCR-pMHC matching during negative selection in thymus")
        #         plt.xlabel('Number of interactions')
        #         plt.ylabel('Number of times T-cell(s) matches closely with pMHC')
        #         plt.legend(['Number of TCR-pMHC matches'])
        #         plt.savefig(str_output_filename_matched + ".png")
        #         plt.close()
        #==============================================================================
        
        
        # print(df_tcell_survive.describe())
        #pdb.set_trace()
    
        # other complex queries
        # from pandasql import *
        # temp_df_2 = sqldf("select * "
        #                       "from df_traffic_times_file"
        #                       )
        # # More complex queries in SQL like format from my website
        # ind = (df_traffic_times_file.from_place == ‘Boston’ )
        # temp = df_traffic_times_file[ind]
        # temp_2 = sqldf("select * "
        #                "from temp where duration_in_traffic > 6200"
        #                )
        # print("Most peaky departure time (Boston)\n”,temp_2.departure_time), #temp_2.departure_time.dt.dayofweek)



        ####################################################################
        # Look at interesting properties:
        #   1) escape probability of auto Tcell
        #   2) how does that scale with the number of TECs and thymocytes
        ####################################################################

        if (first_model.num_tcells_matched == 0):
            # Then TCR has escaped thymus
            print("TCR has escaped thymic negative selection: ", tcell.receptor)

            # append to list of TCR that escape thymus
            list_all_tcell_sequence_escape.append(tcell.receptor)

            if str_checkauto == 'yescheckauto':
                #####################################
                # check additional is this AUTO?
                #####################################
                if str_turnstile_threshold == 'turnstile':
                    b_escape_AUTO = receptor_matching_function_turnstile(tcell_receptor=tcell.receptor,
                                                            array_peptide=list_giant_ALL_COMBINED_peptides,
                                                            i_matching_threshold_within_peptide=i_matching_threshold_within_peptide,
                                                            i_matching_threshold=i_matching_threshold)
                elif str_turnstile_threshold == 'threshold':
                    b_escape_AUTO = receptor_matching_function_threshold(tcell_receptor=tcell.receptor,
                                                                         array_peptide=list_giant_ALL_COMBINED_peptides,
                                                                         i_matching_threshold_within_peptide=i_matching_threshold_within_peptide)

                if b_escape_AUTO:
                    # is actually AUTO but not killed
                    print("is actually AUTO but not killed and escaped:")
                    print(tcell.receptor)

                    # append to list of TCR that escape thymus but are AUTO
                    list_all_tcell_sequence_escape_BUTAUTO.append(tcell.receptor)

    # end for of multiple runs of different TCR each time


    #####################################################################
    # look at multiple simulations and analyze output
    #####################################################################
    print("All TCR that escaped: ", list_all_tcell_sequence_escape)
    aa_freq = analyze_aa_frequencies(list_all_tcell_sequence_escape)
    print("Amino acid frequencies of TCR that escaped: ", aa_freq)
    #pdb.set_trace() # DEBUG code

    #####################################################################
    # serialize output and store to disk for re-analysis offline
    #####################################################################
    # import pickle
    # pickle.dump(aa_freq, open(filename, 'wb'))
    # pickle.load(open(filename, 'rb'))


    #####################################################################
    # PRINTING MODULE
    #####################################################################

    #####################################################################
    # write all sequences of TCR that escape to csv file
    #   and aa frequencies of TCRs that escaped
    #####################################################################
    print("Saving output of all analysis to disk ..")
    with open("tcr_escape_sequence_" + str_output_filename_survival + ".csv", 'w', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        for item in list_all_tcell_sequence_escape:
            wr.writerow(item)
    outfile.close()

    with open("aa_freq_tcrescape_" + str_output_filename_survival + ".csv", 'w', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        for item in [(i,j) for (i,j) in zip(aa_freq[:,0], aa_freq[:,1]) ]:
            wr.writerow(item)
    outfile.close()

    if str_checkauto == 'yescheckauto':
        # if checkauto then save list of escape but auto to file
        with open("escapeBUTAUTO_" + str_output_filename_survival + ".csv", 'w', newline='\n') as outfile:
            wr = csv.writer(outfile)  # , quoting=csv.QUOTE_ALL)
            for item in list_all_tcell_sequence_escape_BUTAUTO:
                wr.writerow(item)
        outfile.close()

    with open("tcr_killed_sequence_" + str_output_filename_survival + ".csv", 'w', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        for item in list_all_tcell_sequence_killed:
            wr.writerow(item)
    outfile.close()

    #####################################################################
    # Analysis of fraction of T-cells that escaped negative selection
    #####################################################################
    print("\n Fraction of T-cells that escaped negative selection \n", 
          len(list_all_tcell_sequence_escape)/i_num_independent_tcells
          )   
    f_fraction_tcells_survived = \
        len(list_all_tcell_sequence_escape)/i_num_independent_tcells

    #####################################################################
    # save these metrics to a file
    #####################################################################
    str_misc_output_file = "misc_modeloutput_" + str_output_filename_survival + ".csv"
    with open(str_misc_output_file, 'w', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        item = "Fraction of T-cells that escaped negative selection"
        wr.writerow(item)
        item = [f_fraction_tcells_survived, f_fraction_tcells_survived]
        wr.writerow(item)
    outfile.close()
    
    
    ####################################################################
    # Testing
    ####################################################################
    # 1) unit tests in test folder 

    ####################################################################
    # Validation
    #   1) which peptides are most likely to be auto?
    #   2) do these make sense clinically/biologically?
    #   3) which TCRs actually survive thymic selection?
    #   4) do these correspond to viruses/pathogens?
    #   5) Does the aa frequency of escaping TCR match Chakrabarti et al. papers
    #           (done in benchmarking code)
    ####################################################################

    #  TCR with receptor YANDSVRKK survives at -ve -44, all 174 cell data
    #   complementary sequence for corresponding peptide is 
    #get_complementary_to_TCR_sequence('YANDSVRKK')
    # LLFFFLLLL
    #notmatch_two_receptors_within_threshold('LLFFFLLLL', 'YANDSVRKK')
    # True
    # TESTING: bug get_complementary_to_TCR_sequence('YANDSVRKK')
    # Issue: https://bitbucket.org/neelsoumya/agent_based_model_mesa/issues/40/bug-get_complementary_to_tcr_sequence
    
    # this fuzzy matches to murine hepatitis virus in 
    # LFVFLLLLP
    # https://www.fludb.org/brc/sssearch.spg?method=ShowCleanInputPage&decorator=corona&preSelectDB=true
    # https://www.fludb.org/brc/proteinSequence.spg?proteinGi=225403261&decorator=corona&range=2-10&context=SS_760020346912
    # https://www.fludb.org/brc/viprDetails.spg?ncbiProteinId=ACN89743&decorator=corona

    # However viral sequence LFVFLLLLP does not match TCR within threshold    
    #notmatch_two_receptors_within_threshold('LFVFLLLLP', 'YANDSVRKK')

    ########################################################################
    # Benchmarking
    #   Does the aa frequency of escaping TCR match Chakrabarti et al. papers
    #   (done in benchmarking code)
    ########################################################################

    # os.system("nohup python3 benchmark_ImmuneModel_nomovement.py 'bench_nocor_noauto_tcrandomprot_1e6tcell_c2inter' 'bench_nocorr_noauto_tcrandomprot_1e6tcell_c2inter' 1000000 True True 2 > /dev/null")
    # and compare with correlations on real data
    # os.system("nohup python3 benchmark_ImmuneModel_nomovement.py 'bench_yescor_noauto_tcrandomuniform_1e6tcell_c2inter' 'bench_yescorr_noauto_tcrandomuniform_1e6tcell_c2inter' 1000000 False False 2 > /dev/null")

    #####################################################################
    # CHALLENGE MODULE
    #####################################################################

    ####################################################################
    # Challenge Tcell that survives negative selection:
    #   1) with pathogen (virus or bacteria)
    #   2) sequences from virus databases
    #       https://www.fludb.org/brc/viprDetails.spg?ncbiProteinId=ACN89743&decorator=corona
    #   3) Disease specific peptides (form Ins1/2 or 21 hydroxylase transcripts) or some peptide database
    #               like http://immunet.cn/bdb/index.php
    ####################################################################

    # str_survived_tcr = return_TCR_survives_mult_rounds()
    # str_virus_peptide = 'LFVFLLLLP'
    # notmatch_two_receptors_within_threshold(str_survived_tcr, str_virus_peptide)
    #print("Challenge T-cells that survive negative selection with pathogen ..")
    #challenge_survivetcellpool_pathogen(list_all_tcell_sequence_escape, str_misc_output_file,
    #                                    f_negative_selection_threshold, i_num_peptide_present, str_turnstile_threshold)

    #print("Challenge T-cells that survive negative selection with disease causing peptides ..")
    #challenge_survivetcellpool_diseasepeptide(list_all_tcell_sequence_escape,
    #                                          str_misc_output_file,
    #                                          f_negative_selection_threshold,
    #                                          i_matching_threshold,
    #                                          str_turnstile_threshold)

    list_generate_pathogen_peptide_pool = \
        challenge_survivetcellpool_pathogen(list_tcr_escape_thymus=list_all_tcell_sequence_escape,
                                            str_misc_output_file_name=str_misc_output_file,
                                            i_matching_threshold_within_peptide=i_matching_threshold_within_peptide,
                                            i_matching_threshold=i_matching_threshold,
                                            i_num_peptide_present=i_num_peptide_present,
                                            str_turnstile_threshold=str_turnstile_threshold)
    ####################################################################
    # TODO: 
    # do any of these escaping sequences match with list of offline precomputed auto TCR
    #   get percentage auto metric
    ####################################################################
    #[list_all_tcell_sequence_escape in get_auto_TCR_sequence()]

    #######################################################
    # Degeneracy of non-self against all TCR repertoire
    #######################################################
    str_degen_diseasepep_againstTCR_filename = "degen_diseasepep_against_TCRrepertoire_" + str_output_filename_survival + ".csv"
    degeneracy_pathogen_against_survivetcellpool_pathogen(list_generate_pathogen_peptide_pool=list_generate_pathogen_peptide_pool,
                                                         list_tcr_escape_thymus=list_all_tcell_sequence_escape,
                                                         str_misc_output_file_name=str_misc_output_file,
                                                         i_matching_threshold_within_peptide=i_matching_threshold_within_peptide,
                                                         i_matching_threshold=i_matching_threshold,
                                                         i_num_peptide_present=i_num_peptide_present,
                                                         str_turnstile_threshold=str_turnstile_threshold,
                                                         str_degen_diseasepep_againstTCR_filename=str_degen_diseasepep_againstTCR_filename)


    #######################################################
    # Degeneracy of SELF against all TCR repertoire
    # WARNING: using same function as above
    #######################################################
    str_degen_self_againstTCR_filename = "degen_SELF_against_TCRrepertoire_" + str_output_filename_survival + ".csv"
    degeneracy_pathogen_against_survivetcellpool_pathogen(list_generate_pathogen_peptide_pool=list_giant_ALL_COMBINED_peptides,
                                                         list_tcr_escape_thymus=list_all_tcell_sequence_escape,
                                                         str_misc_output_file_name=str_misc_output_file,
                                                         i_matching_threshold_within_peptide=i_matching_threshold_within_peptide,
                                                         i_matching_threshold=i_matching_threshold,
                                                         i_num_peptide_present=i_num_peptide_present,
                                                         str_turnstile_threshold=str_turnstile_threshold,
                                                         str_degen_diseasepep_againstTCR_filename=str_degen_self_againstTCR_filename)

    ####################################################################
    # Print parameters of run (again)
    ####################################################################

    print("\nParameters of run")
    print("*****************************************************")
    print("Number of peptides presented on each mTEC: ", i_num_peptide_present)
    print("Auto T-cell: ", b_auto_flag, "\n", 
          "Negative selection Threshold: ", i_matching_threshold_within_peptide, "\n",
          "TCR generated from proteome freq.?: ", b_generate_tcr_random_aa, "\n",
          "mTEC generated from proteome freq.?:", b_generate_mtec_random_aa, "\n",
          "Number of independent T-cells generated:", i_num_independent_tcells, "\n",
          "Number of interactions a single T-cell has with mTECs:", i_num_interactions, "\n",
          "Negative selection energy threshold", i_matching_threshold_within_peptide, "\n",
          "Number of immunological synapses", i_num_peptide_present, "\n",
          "Number of peptides that must match on a single interaction (TCR-pMHC) to negatively select:", i_matching_threshold, "\n"
          "Number of thymic epithelial cells", i_num_mtecs, "\n"
          )
    print("Number of Monte-Carlo runs:", i_num_montecarlo)          
    print("*****************************************************\n")


    #####################################################################
    # PRINTING MODULE
    #####################################################################

    # save these parameters to misc output file (NOTE: append mode here)
    # CAUTION: assumes that file already opened
    with open(str_misc_output_file, 'a', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)

        item = "Parameters of run"
        wr.writerow(item)

        item = "*****************************************************"
        wr.writerow(item)

        item = "Number of immunological synapses on each mTEC: "
        wr.writerow(item)

        item = [i_num_peptide_present, i_num_peptide_present]
        wr.writerow(item)

        item = "Auto T-cell: "
        wr.writerow(item)

        item = [b_auto_flag, b_auto_flag]
        wr.writerow(item)

        item = "Negative selection Threshold: "
        wr.writerow(item)

        item = [i_matching_threshold_within_peptide, i_matching_threshold_within_peptide]
        wr.writerow(item)

        item = "Number of thymic epithelial cells"
        wr.writerow(item)

        item = [i_num_mtecs, i_num_mtecs]
        wr.writerow(item)

        item = "TCR generated from proteome freq.?: "
        wr.writerow(item)

        item = [b_generate_tcr_random_aa, b_generate_tcr_random_aa]
        wr.writerow(item)

        item = "mTEC generated from proteome freq.?:"
        wr.writerow(item)

        item = [b_generate_mtec_random_aa, b_generate_mtec_random_aa]
        wr.writerow(item)

        item = "Number of independent T-cells generated:"
        wr.writerow(item)

        item = [i_num_independent_tcells, i_num_independent_tcells]
        wr.writerow(item)

        item = "Number of interactions a single T-cell has with mTECs:"
        wr.writerow(item)

        item = [i_num_interactions, i_num_interactions]
        wr.writerow(item)
        
        item = "Number of peptides that must match on a single interaction (TCR-pMHC) to negatively select:"
        wr.writerow(item)
        
        item = [i_matching_threshold, i_matching_threshold]
        wr.writerow(item)

        item = "Number of Monte-Carlo runs:"
        wr.writerow(item)

        item = [i_num_montecarlo, i_num_montecarlo]
        wr.writerow(item)

        item = "*****************************************************"
        wr.writerow(item)


    outfile.close()


    #####################################################################
    # VISUALIZATION MODULE
    #####################################################################
    
    #####################################################################
    # Plots for aa frequency histograms
    #####################################################################
    plt.figure()
    plt.hist([int(a) for a in aa_freq[:, 1]], bins=20)
    plt.title("Amino acid frequency distribution of TCR that survive negative selection in thymus")
    plt.xlabel("Amino acid labels")
    plt.ylabel("Frequency")
    plt.savefig("aa_freq_" + str_output_filename_survival + ".png")
    plt.close()
    
    #####################################################################
    # Downstream analysis for aa frequency of escaped TCR
    #####################################################################
    # TODO: pass in aa_freq filename
    # "aa_freq_tcrescape_" + str_output_filename_survival + ".csv"
    #os.system("cp results/analysis_aa_frequency.py .")
    #os.system("python3 analysis_aa_frequency.py " + " 'aa_freq_tcrescape_" + str_output_filename_survival + ".csv' ")
