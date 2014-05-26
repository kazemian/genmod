#!/usr/bin/env python
# encoding: utf-8
"""
strict_genetic_models.py

Genetic models take a family object with individuals and variants and annotates for each variant which models they follow in this family.

The following models are checked:

- Autosomal Dominant(AD)
- Autosomal Dominant De Novo(AD_DN)
- Autosomal Recessive(AR_hom)
- Autosomal Recessive De Novo(AR_DN)
- Autosomal Recesive Compound(AR_comp).

In this model a variant must imply affected status, otherwise it can not be dominant. All sick has to be ay least heterozygote for the variant and all healthy can not have it.

We will assume that each individual that we have information about is present among the individual in self.family.individuals.

Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import os
import sys
from datetime import datetime
from pprint import pprint as pp

from genmod.variants import genotype
from genmod.utils import pair_generator

def check_genetic_models(variant_batch, family, verbose = False, phased = False, proc_name = None):
    #A variant batch is a dictionary on the form {gene_id: {variant_id:variant_dict}}
    # Start by getting the genotypes for each variant:
    individuals = family.individuals.values()
    intervals = variant_batch.pop('haploblocks', {})
    for gene in variant_batch:
        
        compound_candidates = []
        compound_pairs = []
        for variant_id in variant_batch[gene]:
            genotypes = {}
            variant = variant_batch[gene][variant_id]
            for individual in family.individuals:
                try:
                    # print(variant_batch[gene][variant_id])
                    gt_info = variant[individual].split(':')[0]
                except KeyError:# If individual not in variant file
                    if verbose:
                        print('Warning! Individual %s is not in variant file!' % individual)
                    gt_info = './.'
                
                individual_genotype = genotype.Genotype(GT=gt_info)
                genotypes[individual] = individual_genotype
            variant['Genotypes'] = genotypes
            variant['Compounds'] = {}
            # Add information of models followed:
            variant['Inheritance_model'] = {'XR' : False, 'XR_dn' : False, 'XD' : False, 
                                            'XD_dn' : False, 'AD' : False, 'AD_dn' : False, 
                                            'AR_hom' : False, 'AR_hom_dn' : False, 'AR_comp' : False, 
                                            'AR_comp_dn' : False
                                            }
            if gene != '-':
                if check_compound_candidates(variant, family):
                    compound_candidates.append(variant_id)
            # Only check X-linked for the variants in the X-chromosome:
            # For X-linked we do not need to check the other models
            if variant['CHROM'] == 'X':
                
                check_X_recessive(variant, family)
                check_X_dominant(variant, family)
            else:
                # Check the dominant model:
                check_dominant(variant, family)
                # Check the recessive model:
                check_recessive(variant, family)
            
    
        # Now check the compound models:
        if len(compound_candidates) > 1:
            
            compound_pairs = pair_generator.Pair_Generator(compound_candidates)
            
            for pair in compound_pairs.generate_pairs():
                variant_1 = variant_batch[gene][pair[0]]
                variant_2 = variant_batch[gene][pair[1]]
                # We know from check_compound_candidates that all variants are present in all affected
                if not phased:
                    variant_1['Compounds'][pair[1]] = 0
                    variant_2['Compounds'][pair[0]] = 0
                    variant_1['Inheritance_model']['AR_comp'] = True
                    variant_2['Inheritance_model']['AR_comp'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('compound', individual, family, variant_1, variant_2)
                
                elif check_compounds(variant_1, variant_2, family, intervals):
                    variant_1['Compounds'][pair[1]] = 0
                    variant_2['Compounds'][pair[0]] = 0
                    variant_1['Inheritance_model']['AR_comp'] = True
                    variant_2['Inheritance_model']['AR_comp'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('compound', individual, family, variant_1, variant_2)
                    
    return

def check_compound_candidates(variant, family):
    """Sort out the variants that are potential compound candidates. 
        This function is used to reduce the number of potential candidates for the future analysis.
        It will go through all variants in a batch and filter out those variants that not fit the model.
        Creates a new dictionary with variants that is a subset of the batch
        
        Cases:
            Affected:
                - If individual is affected it needs to be variant otherwise it can not be a compound
                - It is ok to be 0/1(het.) or 1/1(hom. alt.) for affected individuals
                
            Healthy:
                - Can not be hom. alt for any variant in a potential compound pair.
                
        Args:
            variant(dict): A dictionay with information about a variant.
            family: A family object with information about the family for this analysis
            
        Returns:
            bool: depending on if the variant is a potential compound candidate according to the
            rules stated above
        
    """
    # This is the case when the variant is located in an uninteresting region:
    
    if not variant.get('comp_candidate',True):
        return False
    
    for individual in family.individuals:
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        if family.individuals[individual].affected:
        # Affected can not be hom. ref. for compounds
            if not individual_genotype.heterozygous:
                return False
        elif not family.individuals[individual].affected:
            if individual_genotype.homo_alt:
                return False
        
    return True

def check_compounds(variant_1, variant_2, family, intervals):
    """Check if two variants of a pair follow the compound heterozygous model. 
        At this stage we know that none of the individuals are homozygote alternative for the variants.
        
    """
    # Check in all individuals what genotypes that are in the trio based of the individual picked.
    for individual in family.individuals:
        genotype_1 = variant_1['Genotypes'].get(individual, genotype.Genotype())
        genotype_2 = variant_2['Genotypes'].get(individual, genotype.Genotype())
        #check if variants are in the same phased interval:
        variant_1_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
        variant_2_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
        
        if family.individuals[individual].healthy:
        # If the family is phased we need to check if a healthy individual have both variants on same allele
        # This individual intervals can not overlap
        # If the variants are not in the same phased interval we can not say that the model is not followed
            if variant_1_interval == variant_2_interval:
                # If the variants are on different alleles it can not be a compound pair:
                if genotype_1.has_variant or genotype_2.has_variant:
                    if genotype_1.allele_1 != '0':
                        if genotype_2.allele_2 != '0':
                            return False
                        if genotype_1.allele_2 != '0':
                            if genotype_2.allele_1 != '0':
                                return False        
            # The case where the individual is affected
        elif family.individuals[individual].affected:
            #If the individual is sick and phased it has to have one variant on each allele
            if variant_1_interval == variant_2_interval:
                if (genotype_1.allele_1 == genotype_2.allele_1) or (genotype_1.allele_2 == genotype_2.allele_2):
                    return False
    
    return True


def check_dominant(variant, family):
    """Check if the variant follows the dominant pattern in this family."""
    for individual in family.individuals: 
        # Check in all individuals what genotypes that are in the trio based of the individual picked.
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        if family.individuals[individual].healthy:# The case where the individual is healthy
            if individual_genotype.has_variant:
                # If the individual is healthy and have a variation on one or both alleles it can not be dominant.
                variant['Inheritance_model']['AD'] = False
                variant['Inheritance_model']['AD_dn'] = False
                return
        elif family.individuals[individual].affected:
            # The case when the individual is sick
            if individual_genotype.genotyped:
                if not individual_genotype.heterozygote:
                # Individual has to be heterozygote i AD can be true
                    variant['Inheritance_model']['AD'] = False
                    variant['Inheritance_model']['AD_dn'] = False
                    return
            # Now the ind is sick and have a variant ≠ ref, check parents for de novo
            if family.individuals[individual].has_parents:
                check_parents('dominant', individual, family, variant)
    return

def check_recessive(variant, family):
    """Check if the variant follows the autosomal recessive pattern in this family."""
    for individual in family.individuals:
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy:
        if family.individuals[individual].healthy:
        # If the individual is healthy and homozygote alt the model is broken.
            if individual_genotype.homo_alt:
                variant['Inheritance_model']['AR_hom'] = False
                variant['Inheritance_model']['AR_hom_dn'] = False
                return
        # The case when the individual is sick:
        elif family.individuals[individual].affected:
        # In the case of a sick individual it must be homozygote alternative for Autosomal recessive to be true.
        # Also, we can not exclude the model if no call.
            if not individual_genotype.homo_alt:
                variant['Inheritance_model']['AR_hom'] = False
                variant['Inheritance_model']['AR_hom_dn'] = False
                return
            #Models are followed but we need to check the parents to see if de novo is followed or not.
            elif family.individuals[individual].has_parents:
                check_parents('recessive', individual, family, variant)
    return

def check_X_recessive(variant, family):
    """Check if the variant follows the x linked heterozygous pattern of inheritance in this family."""
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        
        # The case where the individual is healthy
        if family.individuals[individual].healthy:
            # If individual is healthy and homozygote alternative the variant can not be deleterious:
            if individual_genotype.homo_alt:
                variant['Inheritance_model']['XR'] = False
                variant['Inheritance_model']['XR_dn'] = False
                return
        #The case where the individual is a male
            if family.individuals[individual].sex == 1:
                if individual_genotype.has_variant:
        # If the individual is healthy, male and have a variation it can not be x-linked-recessive.
                    variant['Inheritance_model']['XR'] = False
                    variant['Inheritance_model']['XR_dn'] = False
                    return
        
        # The case when the individual is sick
        elif family.individuals[individual].affected:
        #If the individual is sick it has to be homozygote alt
            if not individual_genotype.homo_alt:
                variant['Inheritance_model']['XR'] = False
                variant['Inheritance_model']['XR_dn'] = False
                return
            
            if family.individuals[individual].has_parents:
                check_parents('X_recessive', individual, family, variant)
    return

def check_X_dominant(variant, family):
    """Check if the variant follows the x linked dominant pattern of inheritance in this family."""
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy
        if not family.individuals[individual].affected():
        # Healthy womans can be carriers but not homozygote:
            if family.individuals[individual].sex == 2:
                if individual_genotype.homo_alt:
                    variant['Inheritance_model']['XD'] = False
                    variant['Inheritance_model']['XD_dn'] = False
                    return
        # Males can not carry the variant:
            elif family.individuals[individual].sex == 1:
                if individual_genotype.has_variant:
                    variant['Inheritance_model']['XD'] = False
                    variant['Inheritance_model']['XD_dn'] = False
                    return
        # The case when the individual is sick
        elif family.individuals[individual].affected:
        #If the individual is sick and homozygote ref it can not be x-linked-dominant
            if individual_genotype.homo_ref:
                variant['Inheritance_model']['XD'] = False
                variant['Inheritance_model']['XD_dn'] = False
                return
            elif individual_genotype.has_variant:
                if family.individuals[individual].has_parents:
                    check_parents('X_dominant', individual, family, variant)
    return

def check_parents(model, individual, family, variant, variant_2={}):
    """Check if information in the parents can tell us if model is de novo or not. Model in ['recessive', 'compound', 'dominant']."""
    sex = family.individuals[individual].sex
    individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())

    mother_id = family.individuals[individual].mother
    mother_genotype = variant['Genotypes'].get(mother_id, genotype.Genotype())
    mother_phenotype = family.get_phenotype(mother_id)

    father_id = family.individuals[individual].father
    father_genotype = variant['Genotypes'].get(father_id, genotype.Genotype())
    father_phenotype = family.get_phenotype(father_id)


    if model == 'recessive':
        # If a parnent is homozygote or if both parents are heterozygote the variant is not denovo
        if ((mother_genotype.homo_alt or father_genotype.homo_alt) or
                (mother_genotype.has_variant and father_genotype.has_variant)):
            variant['Inheritance_model']['AR_hom_dn'] = False
        # If both parents are called but none of the above is fullfilled it is denovo
        elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['AR_hom'] = False
                    
    elif model == 'dominant':
        # If one or both parents are affected it is de novo if none of them have a variant
        if mother_genotype.has_variant or father_genotype.has_variant:
            variant['Inheritance_model']['AD_dn'] = False
        # If both parents are called but none of them carry the variant it is denovo
        elif mother_genotype.genotyped and father_genotype.genotyped:
            variant['Inheritance_model']['AD'] = False
            
    elif model == 'X_recessive':
        #If the individual is a male we only need if the mother carry the variant:
        if sex == 1:
            if mother_genotype.has_variant:
                variant['Inheritance_model']['XR_dn'] = False
            elif mother_genotype.genotyped:
                variant['Inheritance_model']['XR'] = False
        #If female, both parents must have the variant otherwise denovo is true
        elif sex == 2:
            if (mother_genotype.has_variant and father_genotype.has_variant):
                variant['Inheritance_model']['XR_dn'] = False
        #If both parents are genotyped but they both are not carriers XR is not true
            elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['XR'] = False
    
    elif model == 'X_dominant':
        #If the individual is a male we only need to look at the mother:
        if sex == 1:
            if mother_genotype.has_variant:
                variant['Inheritance_model']['XD_dn'] = False
            elif mother_genotype.genotyped:
                variant['Inheritance_model']['XD'] = False
        #If female, one of the parents must have the variant otherwise denovo is true
        elif sex == 2:
            if (mother_genotype.has_variant or father_genotype.has_variant):
                variant['Inheritance_model']['XD_dn'] = False
            elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['XD'] = False
    
    elif model == 'compound':
        individual_genotype_2 = variant_2['Genotypes'].get(individual, genotype.Genotype())
        mother_genotype_2 = variant_2['Genotypes'].get(mother_id, genotype.Genotype())
        father_genotype_2 = variant_2['Genotypes'].get(father_id, genotype.Genotype())
        # One of the variants must come from father and one from mother
        if not (mother_genotype.has_variant or mother_genotype_2.has_variant):
            variant['Inheritance_model']['AR_comp_dn'] = True
            variant_2['Inheritance_model']['AR_comp_dn'] = True
            if mother_genotype.genotyped and mother_genotype_2.genotyped:
                variant['Inheritance_model']['AR_comp'] = False
                variant_2['Inheritance_model']['AR_comp'] = False
        if not (father_genotype.has_variant or father_genotype_2.has_variant):
            variant['Inheritance_model']['AR_comp_dn'] = True
            variant_2['Inheritance_model']['AR_comp_dn'] = True
            if father_genotype.genotyped and father_genotype_2.genotyped:
                variant['Inheritance_model']['AR_comp'] = False
                variant_2['Inheritance_model']['AR_comp'] = False
            
            
        

def main():
    from ped_parser import family, individual
    from interval_tree import interval_tree
    
    duo_family = family.Family(family_id = '1')
    sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
    healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
    healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
    duo_family.add_individual(sick_son)
    duo_family.add_individual(healthy_mother)
    duo_family.add_individual(healthy_father)
    
    pp(duo_family.individuals)
    intervals = {ind_id:interval_tree.IntervalTree([[1,100, '1']], 1, 100) for ind_id in duo_family.individuals}

    pp(intervals)
    
    #Setup two variants with autosomal recessive compound pattern
    recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                                             '1':'0|1', '2':'1|0', '3':'0|0'}
    
    recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                            '1':'1|0', '2':'1|0', '3':'0|0'}
    
    
    batch = {'ABC':{'1_5_A_C':recessive_comp_simple_1, '1_10_C_T':recessive_comp_simple_2}}
    batch['haploblocks'] = intervals
    
    check_genetic_models(batch, duo_family, phased=True)
    for gene in batch:
        pp(batch[gene])
    

if __name__ == '__main__':
    main()

