#!/usr/bin/env python                                                                                                                                                                                                               

__author__ = 'Patrick Campbell'
__email__ = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'

#Simple MONET utility to command line pair model vs. observations                                                                                                                                                                   

import os
from glob import glob
import sys

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet
from monet.util.tools import long_to_wide
import pandas as pd


if __name__ == '__main__':
    parser = ArgumentParser(
        description='combines paired data files for multiple model runs',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-p',
        '--paired_files',
        help='string input paired data files (>=2)',
        nargs='+',
        type=str,
        required=True)
    parser.add_argument(
        '-s',
        '--model_species',
        help='string/list input for obs species-variables to pair',
        type=str,
        nargs='+',
        required=False,
        default=['O3', 'PM25_TOT'])
    parser.add_argument(
        '-mdf',
        '--mergedf',
        help='boolean operator to merge entire dataframes (slow)',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-o',
        '--output',
        help='string output path/filename for combined paired dataframe',
        type=str,
        required=False,
        default='./AIRNOW_CMAQ_merged_pair')
    parser.add_argument(
        '-v',
        '--verbose',
        help='print debugging information',
        action='store_true',
        required=False)
    args = parser.parse_args()

    paired_files = args.paired_files
    species = args.model_species
    output = args.output
    verbose = args.verbose
    mergedf    = args.mergedf

    print('combining paired files...')
    print(paired_files)
    count=0
    for i in paired_files:
        df=pd.read_hdf(i) 
        if count == 0:
                df_merge=df
        else: 
                if mergedf is False:
                	df_species=df[species]
                	df_species = df_species.add_suffix('_'+str(count+1))
                	print('joining model data columns...')
                	df_merge=df_merge.join(df_species)
                else:
                	df[species] = df[species].add_suffix('_'+str(count+1))
                	print('merging entire dataframes...slow...')
                	df_merge=df_merge.merge(df,on='time')
        count=count+1
    print('final merged data frame...')
    print(df_merge.keys())

    df_merge.to_hdf(output+'.hdf',
            'df_merge',
            format='table',
            mode='w')
    df_merge.to_csv(output+'.csv')
    sys.exit(0)
