#!/usr/bin/python
import sys
import subprocess
import shlex
from multiprocessing import Pool
import glob
import os
import re
import argparse
import pandas as pd
import numpy as np
from process import process


#python run_ribosomal_tRNA_matR.py --useR2 ./data/metadata.txt sample.type~cell.line SW480-Par tRNA_LvM2-v-Par.txt
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--printMode", help="Only print commands", dest='runMode', action='store_false')
  parser.add_argument("-r", "--runMode", help="Run commands (default)", dest='runMode', action='store_true')
  parser.add_argument("-l", "--read", help="R1, R2, both; Use both R1 and R2 (miSeq) as opposed to R2 (10X). Default R2", type=str)
  parser.add_argument("-a", "--aligner", help="Choose aligner package (default BWA)", type=str)
  parser.add_argument("--paired", help="Paired end (default)", dest='paired', action='store_true')
  parser.add_argument("--single", help="Single end", dest='paired', action='store_false')
  parser.add_argument("--hasUMI", help="UMI version (default)", dest='umi', action='store_true')
  parser.add_argument("--noUMI", help="Without UMI", dest='umi', action='store_false')
  parser.add_argument("--ref", help="The sample that is the reference in univariate analysis.", type=str)
  parser.add_argument("metadata", help="Specifications of the samples", type=str)
  parser.add_argument("formula", help="Design formula: e.g. type~cell", type=str)
  parser.add_argument("outfile", help="Output file format.", type=str)
  parser.set_defaults(species="human", read="R2", aligner="BWA", paired=True, umi=True, runMode=True)
  args = parser.parse_args()

  matRdir='./'
  if ('matRdir' in os.environ):
    matRdir = os.environ['matRdir']
  else:
    print('matRdir is not set. Please uese `export matRdir=<dir>` to set the directory path for matRdir.')
    exit(1)

  tRNAref = "{}/tRNAs/human/human_mature_tRNA_ref_reduced.fa".format(matRdir)
  if (args.species=="mouse"):
    tRNAref = "{}/tRNAs/mouse/mouse_mature_tRNA_ref_reduced.fa"
  meta = pd.read_csv(args.metadata, sep="\t", header=0, index_col=0)
  print(meta)
  met = os.path.abspath(args.metadata)
  indir = os.path.dirname(met)
  os.chdir(indir)

  log= open("pipeline.txt", "w+")
  runner = process(metadata=meta, matRdir=matRdir, reffile=tRNAref, log=log, isPaired=args.paired, useRead=args.read, hasUMI=args.umi, aligner=args.aligner, runMode=args.runMode)
  print(matRdir, tRNAref, args.paired, args.read, args.umi, args.aligner, args.runMode)
  runner.umi_extract()
  runner.trim()
  runner.merge()
  runner.align()
  runner.dedup()
  runner.count()
  runner.make_sam()
  runner.make_mut_profile()
  runner.count_file()


  analysis_type = "univariate"
  dummy = args.formula.split('~')
  if (len(dummy[0])>0 and len(dummy[1])>0):
    analysis_type = "logit"
  print("Based on the formula, the analysis type is identified as {}".format(analysis_type))

  if (analysis_type == "univariate"):
    runner.univariate(os.path.basename(args.metadata),args.formula,args.ref,args.outfile)
  elif (analysis_type == "logit"):
    runner.logit(os.path.basename(args.metadata),args.formula,args.outfile)

  if(args.formula.split('~')[0]==''):
    covariate = args.formula.split('~')[1]
  else:
    covariate = args.formula.split('~')[0]
  runner.mlogit(os.path.basename(args.metadata),covariate,re.sub('.txt', '.mprofile.txt',args.outfile))
  log.close()
