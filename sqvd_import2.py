#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import os
import string
import sys
import requests
from pysqvd import SQVD
import argparse


# possible extensions:
# - check if study already exists before creating it
# - check if file has been uplaoded before
# - Tell what print statements show

def main(args):
    # extract names from *.vcf files in current folder
    input_file_name=[f for f in os.listdir(".") if f.endswith(".vcf")]
    print(input_file_name, file=sys.stderr) 

    # extract base name from vcf file names and write to basename.txt
    # what do we need the basename.txt file for? Would it make sense to use it as a log, eg when file has been uploaded successfully.
    hold_base_name = []
    with open("basename.txt", "w+") as outputbase:
        for line in input_file_name:
            base = line.split(".")[0]  # not safe if filename contains .
            hold_base_name.append(base)
            outputbase.write(base+"\n")

    print(hold_base_name, file=sys.stderr)

    # define pysqvd submission and upload details with sqvd.create.Study
    with SQVD(username=args.username, password=args.password, host=args.host) as sqvd:
        for basename in hold_base_name:
            new_study = {
              'study_name': basename,
              'sample_id': basename,
              'panel_id': args.panel,
              'panel_version': args.panelversion,
              'subpanels': ['FULL'],
              'workflow': args.workflow,
              'group': args.group
            }
            sqvd.createStudy(new_study)
            print(new_study, file=sys.stderr)  

    # upload VCF files corresponding to studies created above    
    with sqvd:
        for file_basename in hold_base_name:   
            study_id = file_basename
            vcfFile = file_basename + (".vcf")
            sqvd.upload([vcfFile],study_id)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Creates new sqvd studies and adds corresponding vcf files in batch.')
    parser.add_argument('-dir', nargs='?' , default=os.getcwd())
    parser.add_argument('-u', dest='username', help="Username")
    parser.add_argument('-p', dest='password', help="Password")
    parser.add_argument('-ho', dest='host', help="Host")
    parser.add_argument('--panel', help="Panel")
    parser.add_argument('--panelversion', help="Panel Version")
    parser.add_argument('--workflow', default='dna_somatic', help="Workflow")
    parser.add_argument('--group', default='precmend', help="Owning group")

    args = parser.parse_args()
    main(args)
