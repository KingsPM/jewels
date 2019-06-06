#!/usr/bin/env python
# coding: utf-8

import os
import string
import sys
import requests
from pysqvd import SQVD
import argparse

if '-h' in sys.argv[1]:
 parser = argparse.ArgumentParser(
  description='Creates new sqvd studies and add corresponding vcf files from batch of vcf files in cuurent directory. Requires pysqvd to be installed in python environment')
 parser.add_argument('required input:', help='batch_pysqvd.py --usename --password')
 #args = parser.parse_args()
 parser.print_help()
else:

####extract names from *.vcf files in current folder####

 input_file_name=[f for f in os.listdir(".")
  if f.endswith(".vcf")]

 outputbase=open("basename.txt", "w+")
# hold_file_name=[]
 hold_base_name=[]

 print(input_file_name) 

####extract base name from vcf file names####

 for line in input_file_name:
#  hold_file_name.append(line)
#  base=line.spilt(".")
   base=(line.split(".")[0])	
   hold_base_name.append(base)
   outputbase.write(base)
   outputbase.write("\n")	
 outputbase.close()

 print(hold_base_name)

###open connection to sqvd (currently local)###
 sqvd = SQVD(username=sys.argv[1], password=sys.argv[2], host='127.0.0.1:3000')

####### define pysqvd submission and upload details with sqvd.create.Study ###
 with sqvd:
  for basename in hold_base_name:
    new_study = {
              'study_name': basename,
              'sample_id': basename,
              'panel_id': 'MYELOID',
              'panel_version': 1,
              'subpanels': ['FULL'],
   	      'workflow': 'dna_somatic',
              'group': 'precmed'
    }
  
    sqvd.createStudy(new_study)
  print(new_study)  

#####upload VCF files corresponding to studies created above 

 with sqvd:
  for id in hold_base_name:   
    study_id = id
   #for file_name in input_file_id:   
    vcfFile = id + (".vcf")
   
    sqvd.upload([vcfFile],study_id)

