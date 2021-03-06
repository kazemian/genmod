#!/usr/bin/env python
# encoding: utf-8
"""
analyze_variants.py

Analyze the the variants in a vcf, the following will be printed:
    
    - How many variants found
    - How many mendelian violations
    - How many variants where not covered in all individuals. (Default depth 7)
    - How many variants did not satisfy the base call quality treshold. (Default 10)
    - How many variants followed each model:
        - AR_hom
        - AR_comp
        - AR_hom_dn
        - AR_comp_dn
        - AD
        - AD_dn
        - XD
        - XD_dn
        - XR
        - XR_dn
    - How many variants in genetic regions
    - How many rare variants (Default maf < 0.02)
    - How many high scored cadd. (Default cadd > 10)
    - How many rare + high score cadd
    - How many follow a genetic model + rare + high cadd

Created by Måns Magnusson on 2014-09-08.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import click

from codecs import open
from datetime import datetime
from pprint import pprint as pp

import pkg_resources

from vcf_parser import parser as vcf_parser

import genmod
from genmod import warning

###           This is for analyzing the variants       ###

@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
# @click.option('-c', '--config_file',
#                     type=click.Path(exists=True),
#                     help="""Specify the path to a config file."""
# )
@click.option('--frequency_treshold', '-freq',
                    default=0.05,
                    nargs=1,
                    help='Specify the ferquency treshold for variants to be considered. Default 0.05'
)
@click.option('--cadd_treshold', '-cadd',
                    default=10.0,
                    nargs=1,
                    help='Specify the cadd score treshold for variants to be considered. Default 10'
)
@click.option('--gq_treshold', '-gq',
                    default=10.0,
                    nargs=1,
                    help='Specify the genotype quality treshold for variants to be considered. Default 10'
)
@click.option('--read_depth_treshold', '-depth',
                    default=10.0,
                    nargs=1,
                    help='Specify the genotype quality treshold for variants to be considered. Default 10'
)
# @click.option('-p', '--patterns',
#                     type=click.Choice(['AR', 'AD', 'X']),
#                     multiple=True,
#                     help='Specify the inheritance patterns. Default is all patterns'
# )
# @click.option('-o', '--outfile',
#                     type=click.Path(exists=False),
#                     help='Specify the path to a file where results should be stored.'
# )
# @click.option('-v', '--verbose',
#                 is_flag=True,
#                 help='Increase output verbosity.'
# )
def summarize_variants(variant_file, frequency_treshold, cadd_treshold, gq_treshold, read_depth_treshold):
    """Analyze the annotated variants in a VCF file."""    
        
    freq_keyword = '1000GMAF'
    inheritance_keyword = 'GeneticModels'
    
    inheritance_dict = {'AR_hom':0, 'AR_hom_dn': 0, 'AR_comp':0, 'AR_comp_dn': 0 , 'AD':0, 'AD_dn':0, 
                            'XD':0, 'XD_dn':0, 'XR':0, 'XR_dn':0}
    number_of_variants = 0
    rare_variants = 0
    high_cadd_scores = 0
    no_cadd_score = 0
    high_cadd_and_rare = 0
    high_gq = 0
    covered_in_all = 0
    indels = 0
    indel_no_cadd = 0
    true_de_novos = 0
    
    analysis_start = datetime.now()
    
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock = sys.stdin)
    else:
        variant_parser = vcf_parser.VCFParser(infile = variant_file)
    
    for variant in variant_parser:
        models_found = variant['info_dict'].get(inheritance_keyword, None)
        maf = min([float(frequency) for frequency in variant['info_dict'].get(freq_keyword, '0').split(',')])
        cadd_score = max([float(cscore) for cscore in variant['info_dict'].get('CADD', '0').split(',')])
        reference = variant['REF']
        alternative = variant['ALT']

        number_of_variants += 1
        genotypes = variant.get('genotypes', {})
        
        correct_genotype = True
        adequate_depth = True
                
        for individual in genotypes:
            if genotypes[individual].genotype_quality < gq_treshold:
                    correct_genotype = False
            
            #If any individual has depth below "depth" we do not consider the variant
            if genotypes[individual].quality_depth < read_depth_treshold:
                adequate_depth = False
        
        # Check what variant models that are followed:
        if models_found:
            for model in models_found.split(','):
                if '_dn' in model:
                    if correct_genotype:
                        true_de_novos += 1
                        pp(variant)
                inheritance_dict[model] += 1
        
        # Check the frequency of the variants:
        
        if maf < frequency_treshold:
            rare_variants += 1
            if cadd_score > cadd_treshold:
                high_cadd_and_rare += 1
        
        # Check the cadd score:
        
        if cadd_score > cadd_treshold:
            high_cadd_scores += 1
        
        elif cadd_score == 0:
            no_cadd_score += 1
        
        # Check if indel:
        
        if len(reference) > 1 or len(alternative) > 1:
            indels += 1
            if cadd_score == 0:
                indel_no_cadd += 1

    pp(inheritance_dict)    
    print('Number of variants: %s' % number_of_variants)
    print('Number of rare: %s. Frequency of all: %.2f' % (rare_variants, rare_variants/number_of_variants))
    print('Number of high cadd scores: %s. Frequency of all: %.2f' % (high_cadd_scores, high_cadd_scores/number_of_variants))
    print('Number of high cadd scores and rare: %s. Frequency of all: %.2f' 
                % (high_cadd_and_rare, high_cadd_and_rare/number_of_variants))
    print('Number of no cadd scores: %s. Frequency of all: %.2f \n' % (no_cadd_score, no_cadd_score/number_of_variants))
    print('Number of indels: %s. Frequency of all: %.2f' % (indels, indels/number_of_variants))
    print('Number of indels and no cadd score: %s. Frequency of all: %.2f \n' 
                % (indel_no_cadd, indel_no_cadd/number_of_variants))
    print('"True" de novos: %s \n' % true_de_novos)
    print('Time for analysis: %s' % str(datetime.now()-analysis_start))

if __name__ == '__main__':
    summarize_variants()