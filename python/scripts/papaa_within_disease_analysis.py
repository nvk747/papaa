#!/usr/bin/env python3

"""
Pancancer_Aberrant_Pathway_Activity_Analysis
scripts/within_tissue_analysis.py

Usage: Run in command line

        python within_tissue_analysis.py

with the following required flags:

        --genes         comma separated string of HUGO gene symbols

and the following optional flags:

        --diseases      comma separated string of disease types to include
        --l1_ratios     comma separated string of l1 parameters to test
        --folder        string indicating the location to save results
        --remove_hyper  if present, remove hypermutated tumors

Output:
Results of single tissue classifier run through pancancer_classifier.py
"""

import os
import sys
import subprocess
import argparse
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'papaa'))
from tcga_util import add_version_argument

parser = argparse.ArgumentParser()
add_version_argument(parser)
parser.add_argument('-g', '--genes',
                    help='string of the genes to extract or genelist file')
parser.add_argument('-d', '--diseases', default='Auto',
                    help='comma seperated diseases list in a file')
parser.add_argument('-a', '--alphas', default='0.1,0.15,0.2,0.5,0.8,1',
                    help='the alphas for parameter sweep')
parser.add_argument('-l', '--l1_ratios', default='0,0.1,0.15,0.18,0.2,0.3',
                    help='the l1 ratios for parameter sweep')
parser.add_argument('-n', '--num_features', default=8000, type=int,
                        help='Number of MAD genes to include in classifier')
parser.add_argument( '-y','--seed', default=None, type=int,
                        help='option to set seed')
parser.add_argument('-v', '--remove_hyper', action='store_true',
                    help='Remove hypermutated samples')
parser.add_argument('-f', '--classifier_results', default='Auto',
                    help='location to save')
parser.add_argument('-x', '--x_matrix', default=None,
                    help='Filename of features to use in model')
parser.add_argument( '--filename_mut', default=None,
                    help='Filename of sample/gene mutations to use in model')
parser.add_argument( '--filename_mut_burden', default=None,
                    help='Filename of sample mutation burden to use in model')
parser.add_argument( '--filename_sample', default=None,
                    help='Filename of patient/samples to use in model')
parser.add_argument( '--filename_copy_loss', default=None,
                    help='Filename of copy number loss')
parser.add_argument( '--filename_copy_gain', default=None,
                    help='Filename of copy number gain')
parser.add_argument( '--filename_cancer_gene_classification', default=None,
                    help='Filename of cancer gene classification table')

args = parser.parse_args()

if args.seed is not None:
    se = int(args.seed)
    np.random.seed(se)

# make it a little easier to pass forward filename args
args_dict = vars(args)
filename_arg_list = [ 'x_matrix' ]
for k in args_dict.keys():
    if k.startswith('filename_'):
        filename_arg_list.append(k)

# Load command arguments
# if list of the genes provided by file or comma seperated values:
try:
    genes = args.genes
    genes_df = pd.read_table(genes)
    genes = ','.join(genes_df['genes'].tolist())
except:
    genes = args.genes

# if list of the diseases provided by file or comma seperated values:
try:
    diseases = args.diseases
    diseases_df = pd.read_table(diseases)
    diseases = diseases_df['diseases'].tolist()
except:
    diseases = args.diseases.split(',')

folder = args.classifier_results
alphas = args.alphas
l1_ratios = args.l1_ratios
remove_hyper = args.remove_hyper
num_features_kept = args.num_features
expr_file = args.x_matrix
sample_freeze_file = args.filename_sample 

#base_folder = os.path.join('classifiers', 'within_disease',
#                           genes.replace(',', '_'))
base_folder = os.path.join('classifiers', 'within_disease',
                           genes[0])
if diseases == 'Auto':
    sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
    disease_types = sample_freeze['DISEASE'].unique().tolist()
else:
    disease_types = diseases

# Loop over disease types
for acronym in disease_types:
    print(acronym)
    if folder == 'Auto':
        result_folder = os.path.join(base_folder, acronym)
    else:
        result_folder = os.path.join(folder, 'within_disease', acronym)
   
    command = ['papaa_pancancer_classifier.py',
               '--genes', genes, '--diseases', acronym, '--drop',
               '--seed', str(se) , '--num_features', str(num_features_kept),
               '--copy_number', '--alphas', alphas, '--l1_ratios', l1_ratios,
               '--classifier_results', result_folder, '--shuffled', '--keep_intermediate']
    if remove_hyper:
        command += ['--remove_hyper']

    # Only set filename if it has been set
    for fname_arg in filename_arg_list:
        val = args_dict.get(fname_arg, None)   
        if val:
            command += ['--%s' % (fname_arg), val]
    subprocess.call(command)
