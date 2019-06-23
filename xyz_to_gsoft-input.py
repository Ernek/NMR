#!/usr/bin/env python
# coding: utf-8

# Code to transform an xyz snapshot file to Gaussian input exploring several Basis Sets and QM methods 

#%%
# Importing modules 


import sys, os 
import numpy as np 
import pandas as pd 
import math 
import fnmatch    # To read files 
from django.utils.text import slugify    # To convert special characters to valid file_path names 
import itertools 
from unidecode import unidecode
import re
import argparse  
import basis_set_exchange as bse



#%%

# Functions

def u_slugify(txt):
        """A custom version of slugify that retains non-ascii characters. The purpose of this
        function in the application is to make folder names readable by command line tools."""
        txt = txt.strip() # remove trailing whitespace
        #txt = re.sub('\s*-\s*','-', txt, re.UNICODE) # remove spaces before and after dashes
        #txt = re.sub('[\s/]', '_', txt, re.UNICODE) # replace remaining spaces with underscores
        #txt = re.sub('(\d):(\d)', r'\1-\2', txt, re.UNICODE) # replace colons between numbers with dashes
        txt = re.sub('\*', '-pol', txt, re.UNICODE) # replace * for text pol=polarization funcitons
        txt = re.sub('\+', '-dif', txt, re.UNICODE) # replace + for text dif=diffuse funcitons
        txt = re.sub('\(', '-', txt, re.UNICODE) # replace left parenthesis with dash 
        txt = re.sub(r'[?,:!@#~`=$%^&\)\[\]{}<>]','',txt, re.UNICODE) # remove some characters altogether
        return txt

def get_bse_local(loc_basis, element):
        """A function to extract element basis sets data from local EMSL database located 
        in /Users/ernesto/Main/Programming/Git/emsl_basis_set_library/gbs """
        addrs = '/Users/ernesto/Main/Programming/Git/emsl_basis_set_library/gbs/'
        with open(addrs + str(loc_basis) + '.gbs', 'r') as fbasis:
            lines = fbasis.readlines()
            #print(fbasis.readlines())
            #tlines = fbasis.readlines()
            bs_list = []
            #tlines = fbasis.readlines()
            #print(tlines)
            for i in range(len(lines)):
                #print(index, line.split()[0])
                if lines[i].split()[0] == str('-') + element:
                    #print('ok', lines[i])
                    #print(i)
                    oidx = i
                    #print(oidx)
                    #print(lines[oidx])
                    oline = lines[oidx]
                    #print(oline)
                    while not oline.split()[0] == '****':
                        bs_list.append(oline)
                        oidx += 1
                        oline = lines[oidx]
                else:
                    continue
                        
        return bs_list
    

    

#%%
        
def main(fname, charge, spin, key_run_01):
    
    
    
    
    # Calculation parameters
    
    
    
    #charge = -1
    #spin = 1
    #key_run_01 = 'nmr'
    key_method_01 = 'giao'
    
    
    
    
    #Create dictionary with QM methods to be used in Gaussian calculation
    
      
        
        
    #print(plus)
    
    functionals = {'hf': ['hf'],
                   'mp2': ['mp2'],
                   'hybrid': ['b3lyp', 'pbe0'],
                   'gga': ['pbe','blyp']
                  }
    #Create dictionary with basis sets to explore with each QM method
    basis_sets = {'pople': ['6-31G', '6-31G(d)','6-31+G', '6-31++G'],
                  'dunning':[],
                  'jensen': [],
                  'core_valence': []
                 }
    
    pople_basis_no_d = list(itertools.product(['6-31'],['G'],['d','2d','3d','3df'],['p','2p','3p','3pd']))
    add_val = []
    for k in pople_basis_no_d:
        add_val.append(str(k[0])+k[1]+'('+k[2]+','+k[3]+')')
        add_val.append(str(k[0])+'+'+k[1]+'('+k[2]+','+k[3]+')')
        add_val.append(str(k[0])+'++'+k[1]+'('+k[2]+','+k[3]+')')  
        print(add_val)
    
    for key,val in basis_sets.items():
        if key == 'pople':
            for p in add_val:
                if p not in val:
                    basis_sets[key].append(p)
                else:
                    continue
    dunning_basis_no_d = list(itertools.product(['cc-'],['pV'],['DZ','TZ','QZ','5Z']))
    add_val = []
    for k in dunning_basis_no_d:
        add_val.append(str(k[0])+k[1]+k[2])
        add_val.append('aug-' + str(k[0])+k[1]+k[2])
#        print(add_val)
        
    for key,val in basis_sets.items():
        if key == 'dunning':
            for p in add_val:
                if p not in val:
                    basis_sets[key].append(p)
                else:
                    continue
        else:
            continue
    
    jensen_basis = list(itertools.product(['pc'],['-' + str(i) for i in range(5)]))
    add_val = []
    for k in jensen_basis:
        add_val.append(str(k[0]+k[1]))
        add_val.append(str(k[0]+'s'+k[1]))
        add_val.append(str(k[0]+'sseg'+k[1]))
    
    for l in range(len(add_val)):
        add_val.append('aug-' + add_val[l])
        
    for key,val in basis_sets.items():
        if key == 'jensen':
            for p in add_val:
                if p not in val:
                    basis_sets[key].append(p)
                else:
                    continue
        else:
            continue
    
    core_valence_basis = list(itertools.product(['cc-pCV'], ['DZ','TZ','QZ','5Z']))
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
            
    #print(u_slugify('6-311++G(d,p)'))
    #print(basis_sets)    
    ###print(functionals.items())
    
    
    #%%
    
    
    #for key, val in basis_sets.items():
    #    for h in val:
    #        print(u_slugify(h))
    
    
    #%%
    
    
    cwd = os.getcwd()  # get Current Working Directory
    #directory = os.path.dirname(fpath)
    #fname = 'xyz_file.xyz'
    #print(cwd)
    
    
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
                #print(str(cwd)+'/'+i)
                os.makedirs(str(cwd)+'/'+i)
            except FileExistsError:
                #directory already exists
                pass
            
            for bkey,bval in basis_sets.items():
                try:
                    #print(str(cwd)+ '/' + i + '/' + bkey)
                    os.makedirs(str(cwd)+ '/' + i + '/' + bkey)
                except FileExistsError:
                    #directory already exists
                    pass
                if bkey == 'jensen':  #If the basis sets is Jensen type, read it from the database and put it explicitely at the end of the input
                    for j in bval:
                        try:
                            #print(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        except FileExistsError:
                            #directory already exists
                            pass
                        out_fname = str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '.com'
                        os.chdir(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        with open(out_fname, 'w') as output:
                            output.write('%chk=' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '\n'
                                        + f"#P {key_run_01}={key_method_01} {i}/Gen"                              
                                        + '\n\n' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '\n\n' 
                                        + str(charge) + ' ' + str(spin) + '\n'                             
                                        + "".join(map(str, a)) + '\n'
                                        + bse.get_basis(j, fmt='gaussian94', elements=['H','O','Al'],header=False)
                                        + '\n'
                                    )
                        
                        sub_fname = 'sub_g09.slurm'
                        with open(sub_fname, 'w') as sub_out:
                            sub_out.write('#!/bin/bash' +'\n'                                   
                                          +'#SBATCH --job-name=' + i + '_' + slugify(j) + '   ###Job Name' + '\n'
                                          +'#SBATCH --partition=clark,kamiak,cas  ###Partition on which to run' + '\n'
                                          +'#SBATCH --nodes=1           ###Number of nodes to use' +'\n'
                                          +'#SBATCH --ntasks-per-node=1 ###Number of tasks per node (aka MPI processes)' +'\n'
                                          +'#SBATCH --cpus-per-task=20   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                          +'#SBATCH --time=7-00:00:00' + '\n'
                                          +'##SBATCH --array=1-500  ##Array indexes ' + '\n'
                                          +'module load gaussian' + '\n'
                                          +'# Name of your com file' +'\n'
                                          +'#JobFile=$1' +'\n'
                                          +'\n'
                                          +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '"' + '\n'
                                          +'fend=' + '"' + '.com' + '"' + '\n'
                                          +'foutend=' + '"'+ '.out' + '"' +'\n'
                                          +'#index=$1' +'\n'
                                          +'\n'
                                          +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                          +'#g09 < ${1}.com > ${1}.out' +'\n'     
                                          +'#g09 < ${finit}${index}${fend} > ${finit}${index}${foutend}' + '\n'      
                                          +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                          +'\n'
                                        )
                        os.chdir(cwd)
                elif  bkey == 'dunning':
                    for j in bval:
                        try:
                            #print(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        except FileExistsError:
                            #directory already exists
                            pass
                        
                        al_bse = get_bse_local(j, 'Al')
                        
                        dplus = str(j[:-2]+'_'+j[-2]+'+d_'+j[-1])
                        al_bse_dplus = get_bse_local(dplus,'Al')
                        
                        try:
                            #print(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + dplus)
                        except FileExistsError:
                            #directory already exists
                            pass
                        
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
                                              +'#SBATCH --job-name=' + i + '_' + u_slugify(j) + '_' + u_slugify(p) + '   ###Job Name' + '\n'
                                              +'#SBATCH --partition=clark,kamiak,cas  ###Partition on which to run' + '\n'
                                              +'#SBATCH --nodes=1           ###Number of nodes to use' +'\n'
                                              +'#SBATCH --ntasks-per-node=1 ###Number of tasks per node (aka MPI processes)' +'\n'
                                              +'#SBATCH --cpus-per-task=20   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                              +'#SBATCH --time=7-00:00:00' + '\n'
                                              +'##SBATCH --array=1-500  ##Array indexes ' + '\n'
                                              +'module load gaussian' + '\n'
                                              +'# Name of your com file' +'\n'
                                              +'#JobFile=$1' +'\n'
                                              +'\n'
                                              +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) +'_' + u_slugify(p) + '"' + '\n'
                                              +'fend=' + '"' + '.com' + '"' + '\n'
                                              +'foutend=' + '"'+ '.out' + '"' +'\n'
                                              +'#index=$1' +'\n'
                                              +'\n'
                                              +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                              +'#g09 < ${1}.com > ${1}.out' +'\n'     
                                              +'#g09 < ${finit}${index}${fend} > ${finit}${index}${foutend}' + '\n'      
                                              +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                              +'\n'
                                            )
                        os.chdir(cwd)
                        
                        for p in [basis_sets['dunning'][indx] for indx,t in enumerate(basis_sets['dunning']) if t[-2]!='5']:
                            out_fname = str(fname).rpartition('.')[0] + '_' + i + '_' + dplus + '_' + u_slugify(p) + '.com'
                            os.chdir(str(cwd)+ '/' + i + '/' + bkey + '/' + dplus)
                            o_bse = get_bse_local(p, 'O')
                            h_bse = get_bse_local(p, 'H')
                            with open(out_fname, 'w') as output:
                                output.write('%chk=' + str(fname).rpartition('.')[0] + '_' + i + '_' + dplus +'_' + u_slugify(p) + '\n'
                                            + f"#P {key_run_01}={key_method_01} {i}/Gen"                              
                                            + '\n\n' + str(fname).rpartition('.')[0] + '_' + i + '_' + dplus +'_' + u_slugify(p) + '\n\n' 
                                            + str(charge) + ' ' + str(spin) + '\n'                             
                                            + "".join(map(str, a))
                                            + '\n'
                                            + "".join(map(str,al_bse_dplus))
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
                                              +'#SBATCH --job-name=' + i + '_' + dplus + '_' + u_slugify(p) + '   ###Job Name' + '\n'
                                              +'#SBATCH --partition=clark,kamiak,cas  ###Partition on which to run' + '\n'
                                              +'#SBATCH --nodes=1           ###Number of nodes to use' +'\n'
                                              +'#SBATCH --ntasks-per-node=1 ###Number of tasks per node (aka MPI processes)' +'\n'
                                              +'#SBATCH --cpus-per-task=20   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                              +'#SBATCH --time=7-00:00:00' + '\n'
                                              +'##SBATCH --array=1-500  ##Array indexes ' + '\n'
                                              +'module load gaussian' + '\n'
                                              +'# Name of your com file' +'\n'
                                              +'#JobFile=$1' +'\n'
                                              +'\n'
                                              +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + dplus +'_' + u_slugify(p) + '"' + '\n'
                                              +'fend=' + '"' + '.com' + '"' + '\n'
                                              +'foutend=' + '"'+ '.out' + '"' +'\n'
                                              +'#index=$1' +'\n'
                                              +'\n'
                                              +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                              +'#g09 < ${1}.com > ${1}.out' +'\n'     
                                              +'#g09 < ${finit}${index}${fend} > ${finit}${index}${foutend}' + '\n'      
                                              +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                              +'\n'
                                            )
                        os.chdir(cwd)
                elif  bkey == 'core_valence':
                    for j in bval:
                        try:
                            #print(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        except FileExistsError:
                            #directory already exists
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
                                              +'#SBATCH --ntasks-per-node=1 ###Number of tasks per node (aka MPI processes)' +'\n'
                                              +'#SBATCH --cpus-per-task=20   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                              +'#SBATCH --time=7-00:00:00' + '\n'
                                              +'##SBATCH --array=1-500  ##Array indexes ' + '\n'
                                              +'module load gaussian' + '\n'
                                              +'# Name of your com file' +'\n'
                                              +'#JobFile=$1' +'\n'
                                              +'\n'
                                              +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) +'_' + u_slugify(p) + '"' + '\n'
                                              +'fend=' + '"' + '.com' + '"' + '\n'
                                              +'foutend=' + '"'+ '.out' + '"' +'\n'
                                              +'#index=$1' +'\n'
                                              +'\n'
                                              +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                              +'#g09 < ${1}.com > ${1}.out' +'\n'     
                                              +'#g09 < ${finit}${index}${fend} > ${finit}${index}${foutend}' + '\n'      
                                              +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                              +'\n'
                                            )
                        os.chdir(cwd)
                else:
                    for j in bval:
                        try:
                            #print(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                            os.makedirs(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        except FileExistsError:
                            #directory already exists
                            pass
                        out_fname = str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '.com'
                        os.chdir(str(cwd)+ '/' + i + '/' + bkey + '/' + u_slugify(j))
                        with open(out_fname, 'w') as output:
                            output.write('%chk=' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '\n'
                                        + f"#P {key_run_01}={key_method_01} {i}/{j}"                              
                                        + '\n\n' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '\n\n' 
                                        + str(charge) + ' ' + str(spin) + '\n'                             
                                        + "".join(map(str, a)) + '\n'
                                    )
                            
                        sub_fname = 'sub_g09.slurm'
                        with open(sub_fname, 'w') as sub_out:
                            sub_out.write('#!/bin/bash' +'\n'                                   
                                          +'#SBATCH --job-name=' + i + '_' + slugify(j) + '   ###Job Name' + '\n'
                                          +'#SBATCH --partition=clark,kamiak,cas  ###Partition on which to run' + '\n'
                                          +'#SBATCH --nodes=1           ###Number of nodes to use' +'\n'
                                          +'#SBATCH --ntasks-per-node=1 ###Number of tasks per node (aka MPI processes)' +'\n'
                                          +'#SBATCH --cpus-per-task=20   ###Number of cpus per task (aka OpenMP threads)' +'\n'
                                          +'#SBATCH --time=7-00:00:00' + '\n'
                                          +'##SBATCH --array=1-500  ##Array indexes ' + '\n'
                                          +'module load gaussian' + '\n'
                                          +'# Name of your com file' +'\n'
                                          +'#JobFile=$1' +'\n'
                                          +'\n'
                                          +'finit=' + '"' + str(fname).rpartition('.')[0] + '_' + i + '_' + u_slugify(j) + '"' + '\n'
                                          +'fend=' + '"' + '.com' + '"' + '\n'
                                          +'foutend=' + '"'+ '.out' + '"' +'\n'
                                          +'#index=$1' +'\n'
                                          +'\n'
                                          +'export GAUSS_SCRDIR="$(mkworkspace -q)" ' +'\n'
                                          +'#g09 < ${1}.com > ${1}.out' +'\n'     
                                          +'#g09 < ${finit}${index}${fend} > ${finit}${index}${foutend}' + '\n'      
                                          +'g09 < ${finit}${fend} > ${finit}${foutend}' + '\n'                     
                                          +'\n'
                                        )
                        os.chdir(cwd)
    os.chdir(cwd)
    
    
#%%    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Enter the xyz filename, i.e "aloh4.xyz" ')
    parser.add_argument('charge', help='Enter the system total charge, i.e "-1" or "1" ', type=int)
    parser.add_argument('spin', help='Enter spin state of the whole system i.e "1" ', type=int)
    parser.add_argument("-m", "--method", help='Enter job directive, i.e "nmr";' \
                        +'if no arg passed defaults to single point energy calc', type=str, choices=['nmr'])
    #parser.add_argument("method", help='Enter job directive, i.e "nmr";' \
    #                    +'if no arg passed defaults to single point energy calc', type=str)
    args = parser.parse_args()
    
    if args.method:
        print(f"Method selected is: {args.method}")
    else: 
        print("Default method. Generating input with single point energy calculation")
    #print(args)
    #print(args.charge)
    main(args.filename, args.charge, args.spin, args.method)


