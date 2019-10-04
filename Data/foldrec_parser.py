#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 11:56:25 2019

@author: c.papadopoulos
"""


#######################
#     FUNCTIONS       #
#######################


def parse_foldrec(file):
    '''
    Author: Chris Papadopoulos
    Last modified : 04 / 10 / 2019
    This function parses a foldrec file and keeps the important information 
    into a dictionary with different keys. 
    You just need to give as input the path of the foldrec file.
    
    PS: It is possible that you need to modify this function based on your 
        project's special needs. 
    '''
    # We read the foldrec file
    with open (file,"r") as inp_foldrec:
        dico = {}
        rank = 0
        # Line-by-line
        for line in inp_foldrec:
            # We keep the basic hit infos
            if line.startswith('Alignment :'):
                # Counter of hits
                rank = rank +1
                dico[rank]= {}
                name  = line.split('vs')[0].split(':')[1].split(',')[0].strip()
                Fam   = line.split('vs')[1].split(':')[0].strip()
                SCOP  = line.split('vs')[1].split(':')[1].strip()
                dico[rank]['name'] = name
                dico[rank]['SCOP'] = SCOP
                dico[rank]['Fam']  = Fam
                #print(IGORF,PU)
            # We keep the alignment infos
            if line.startswith('Score :'):
                # Here we clean the empty strings from the list
                clean_list = filter(None,line.split(' '))
                clean_list_to_use = []
                for elem in clean_list:
                    clean_list_to_use.append(elem)
                # We asign every element to a variable
                Score = clean_list_to_use[2]
                Norm_score = clean_list_to_use[7]
                Query_coverage = clean_list_to_use[12][0:-1]
                Identity = clean_list_to_use[16][0:-1]
                Alignment_length = clean_list_to_use[30]
                dico[rank]['Score'] = float(Score)
                dico[rank]['Norm_score'] = float(Norm_score)
                dico[rank]['Query_coverage'] = float(Query_coverage)
                dico[rank]['Identity'] = float(Identity)
                dico[rank]['Alignment_length'] = int(Alignment_length)
                
            # We keep the Query infos
            if line.startswith('Query'):
                clean_list = filter(None,line.split(' '))
                clean_list_to_use = []
                for elem in clean_list:
                    clean_list_to_use.append(elem)
                    
                query_start = clean_list_to_use[1].rstrip()
                query_stop  = clean_list_to_use[3].rstrip()
                dico[rank]['query_start'] = int(query_start)
                dico[rank]['query_stop']  = int(query_stop)
                query_seq   = clean_list_to_use[2]
                
                if 'W' in sorted(set(list(query_seq))) or \
                   'P' in sorted(set(list(query_seq))) or \
                   'N' in sorted(set(list(query_seq))) or \
                   'R' in sorted(set(list(query_seq))) or \
                   'Y' in sorted(set(list(query_seq))):
                       query_AA = query_seq
                       dico[rank]['query_AA'] = query_AA
                
                #elif query_seq[0] in ['0','1','2','3','4','5','6','7','8','9']:
                elif '9' in sorted(set(list(query_seq))) or \
                     '8' in sorted(set(list(query_seq))) or \
                     '7' in sorted(set(list(query_seq))) or \
                     '6' in sorted(set(list(query_seq))) or \
                     '5' in sorted(set(list(query_seq))) or \
                     '4' in sorted(set(list(query_seq))) or \
                     '3' in sorted(set(list(query_seq))) or \
                     '2' in sorted(set(list(query_seq))) or \
                     '1' in sorted(set(list(query_seq))):
                         confidence = query_seq
                         dico[rank]['query_confidence'] = confidence
                else:
                    psipred_seq = query_seq
                    dico[rank]['query_PSIPRED'] = psipred_seq
                    #print('Here is the AA:',PU,query_AA)
            
            # We keep the template infos
            if line.startswith('Template'):
                clean_list = filter(None,line.split(' '))
                clean_list_to_use = []
                for elem in clean_list:
                    clean_list_to_use.append(elem)
                    
                template_start = clean_list_to_use[1].rstrip()
                template_stop  = clean_list_to_use[3].rstrip()
                dico[rank]['template_start'] = int(template_start)
                dico[rank]['template_stop']  = int(template_stop)
                template_seq   = clean_list_to_use[2]
                
                # I try to fund some aminoacids that which are not C,E,H,G
                # If there is problem add some more to ensure the output
                if 'W' in sorted(set(list(template_seq))) or \
                   'P' in sorted(set(list(template_seq))) or \
                   'N' in sorted(set(list(template_seq))) or \
                   'R' in sorted(set(list(template_seq))) or \
                   'Y' in sorted(set(list(template_seq))):
                       template_AA = template_seq
                       dico[rank]['template_AA'] = template_AA
                
                #elif template_seq[0] in ['0','1','2','3','4','5','6','7','8','9']:
                elif '9' in sorted(set(list(template_seq))) or \
                     '8' in sorted(set(list(template_seq))) or \
                     '7' in sorted(set(list(template_seq))) or \
                     '6' in sorted(set(list(template_seq))) or \
                     '6' in sorted(set(list(template_seq))) or \
                     '5' in sorted(set(list(template_seq))) or \
                     '4' in sorted(set(list(template_seq))) or \
                     '3' in sorted(set(list(template_seq))) or \
                     '2' in sorted(set(list(template_seq))) or \
                     '1' in sorted(set(list(template_seq))):
                
                         confidence = template_seq
                         dico[rank]['template_confidence'] = confidence
                else:
                    psipred_seq = template_seq
                    dico[rank]['template_PSIPRED'] = psipred_seq
                    #print('Here is the AA:',PU,query_AA)

    return(dico)
