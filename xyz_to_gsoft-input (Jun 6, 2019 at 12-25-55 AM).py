#!/usr/bin/env python
# coding: utf-8

# Code to transform an .xyz snap file to Gaussian input exploring several Basis Sets and QM methods 

#%%
# Importing modules ######################################################################################

import sys, os 
import numpy as np 
import pandas as pd 
import math 
import fnmatch    			# To read files 
from django.utils.text import slugify   # To convert special characters to valid file_path names 
import itertools 
from unidecode import unidecode
import re
import argparse  

#%%

# Functions   ############################################################################################

def u_slugify(txt):
        """A custom version of slugify that retains non-ascii characters. The purpose of this
        function in the application is to make folder names readable by command line tools."""
        txt = txt.strip() # remove trailing whitespace
        txt = re.sub('\*', '-pol', txt, re.UNICODE) # replace * for text pol=polarization funcitons
        txt = re.sub('\+', '-dif', txt, re.UNICODE) # replace + for text dif=diffuse funcitons
        txt = re.sub('\(', '-', txt, re.UNICODE) # replace left parenthesis with dash 
        txt = re.sub(r'[?,:!@#~`=$%^&\)\[\]{}<>]','',txt, re.UNICODE) # remove some characters altogether
        return txt

def get_bse_local(loc_basis, element):
        """A function to extract element basis sets data from local EMSL database located 
        in /data/clark/ernesto/bss/emsl_basis_set_library/gbs """
        addrs = '/data/clark/ernesto/bss/emsl_basis_set_library/gbs/'
        with open(addrs + str(loc_basis) + '.gbs', 'r') as fbasis:
            lines = fbasis.readlines()
            bs_list = []
            for i in range(len(lines)):
                if lines[i].split()[0] == str('-') + element:
                    oidx = i
                    oline = lines[oidx]
                    while not oline.split()[0] == '****':
                        bs_list.append(oline)
                        oidx += 1
                        oline = lines[oidx]
                else:
                    continue
        return bs_list

#%%

# Main Function -- Program  #########################################################################
        
def main(fname, charge, spin, key_run_01):
    
    key_method_01 = 'giao'
    
    #Create dictionary with QM methods to be used in Gaussian calculation
    functionals = {'hf': ['hf'],
                   'mp2': ['mp2'],
                   'gga': ['pbe','blyp']
                  }
    #Create dictionary with basis sets to explore with each QM method
    basis_sets = {'core_valence': [],
		  'dunning': ['cc-pVDZ','cc-pVTZ','cc-pVQZ']
                 }
    
    core_valence_basis = list(itertools.product(['cc-pCV'], ['DZ','TZ','QZ']))
    add_val = []
    for k in core_valence_basis:
        add_val.append(str(k[0]+k[1]))
        add_val.append('aug-'+ str(k[0]+k[1]))

    for key,val in basis_sets.items():
        if key == 'core_valence':
            for p in add_val:
                if p not in val:
                    basis_sets[key].append(p)
                else:
                    continue
        else:
            continue
    
    cwd = os.getcwd()  # get Current Working Directory
    
    # Initializing input file reading and outputting file
    with open(fname, 'r') as f:
        cnt = 0
        a = []
        for line in f:
            if cnt == 0:
                N = line.split()[0]
                cnt += 1 
            elif cnt == 1:
                cnt += 1
            else:
                a.append(line)
                cnt += 1 
                
    for fkey,fval in functionals.items():
        for i in fval:
            try:
                os.makedirs(str(cwd)+'/'+i)
            except FileExistsError:
                pass
            
            for bkey,bval in basis_sets.items():
                try:
                    os.makedirs(str(cwd)+ '/' + i + '/' + bkey)
                except FileExistsError:
                    pass
                if  bkey == 'core_valence':
                    for j in bval:
                        try:
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        except FileExistsError:
                            pass
                        
                        al_bse = get_bse_local(j, 'Al')
                        
                        for p in [basis_sets['dunning'][indx] for indx,t in enumerate(basis_sets['dunning']) if t[-2]!='5']:
                            out_fname = str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '_' + u_slugify(p) + '.com'
                            os.chdir(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            o_bse = get_bse_local(p, 'O')
                            h_bse = get_bse_local(p, 'H')
                            with open(out_fname, 'w') as output:
                                output.write('%chk=' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) +'_' + u_slugify(p) + '\n'
                                            + f"#P {key_run_01}={key_method_01} {i}/Gen"                              
                                            + '\n\n' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) +'_' + u_slugify(p) + '\n\n' 
                                            + str(charge) + ' ' + str(spin) + '\n'                             
                                            + "".join(map(str, a))
                                            + '\n'
                                            + "".join(map(str,al_bse))
                                            + '****\n'
                                            + "".join(map(str,o_bse))
                                            + '****\n'
                                            + "".join(map(str,h_bse))
                                            + '****\n'
                                            + '\n'
                                        )
                            
                            sub_fname = 'sub_g09'+'_'+ u_slugify(p)+'.slurm'
                            with open(sub_fname, 'w') as sub_out:
                                sub_out.write('#!/bin/bash' +'\n'                                   
                                              +'#SBATCH --job-name=' + i + '_' + slugify(j) + '_' + u_slugify(p) + '   ###Job Name' + '\n'
                                              +'#SBATCH --partition=clark,kamiak,cas  ###Partition on which to run' + '\n'
                                              +'#SBATCH --nodes=1           ###Number of nodes to use' +'\n'
                                              +'#SBATCH --ntasks-per-node=20 ###Number of tasks per node (aka MPI processes)' +'\n'
                                              +'#SBATCH --cpus-per-task=1   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                              +'#SBATCH --time=7-00:00:00' + '\n'
                                              +'module load gaussian' + '\n'
                                              +'\n'
                                              +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) +'_' + u_slugify(p) + '"' + '\n'
                                              +'fend=' + '"' + '.com' + '"' + '\n'
                                              +'foutend=' + '"'+ '.out' + '"' +'\n'
                                              +'\n'
                                              +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                              +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                              +'\n'
                                            )
                        os.chdir(cwd)
                else:
                    continue
    os.chdir(cwd)
    
#%%    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Enter the xyz filename, i.e "aloh4.xyz" ')
    parser.add_argument('charge', help='Enter the system total charge, i.e "-1" or "1" ', type=int)
    parser.add_argument('spin', help='Enter spin state of the whole system i.e "1" ', type=int)
    parser.add_argument("-m", "--method", help='Enter job directive, i.e "nmr";' \
                        +'if no arg passed defaults to single point energy calc', type=str, choices=['nmr'])
    args = parser.parse_args()
    
    if args.method:
        print(f"Method selected is: {args.method}")
    else: 
        print("Default method. Generating input with single point energy calculation")
    main(args.filename, args.charge, args.spin, args.method)

