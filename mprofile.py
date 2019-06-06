import subprocess
import re
import os
import pandas as pd
from collections import defaultdict

class mprofile:
  def __init__(self, metadata, matRdir, reffile, log, hasUMI=True, runMode=True):
    self.metadata = metadata
    self.matRdir = matRdir
    self.reffile = reffile
    self.umi = hasUMI
    self.log = log
    self.counts=defaultdict(dict)
    self.rm = runMode #print mode if set to false
  
  def pileup(self):
    meta = self.metadata
    log = self.log
    if (self.umi):
      cmd = 'bcftools mpileup --skip-indels --annotate FORMAT/AD,FORMAT/DP,INFO/AD -f {} {} > pileup.txt'.format(self.reffile, " ".join([sample+".dd.bam" for sample in meta.index.tolist()]))
    else:
      cmd = 'bcftools mpileup --skip-indels --annotate FORMAT/AD,FORMAT/DP,INFO/AD -f {} {} > pileup.txt'.format(self.reffile, " ".join([sample+".srt.bam" for sample in meta.index.tolist()]))
    print(cmd)
    log.write("%s\n" % (cmd))
    if self.rm: subprocess.call(cmd,shell=True)
    log.write("\n\n")
  
  def mcount(self):
    meta = self.metadata
    log = self.log
    with open('pileup.txt', 'rt') as vcf:
      for l in vcf:
        if (l.startswith('##')): continue
        
        if (l.startswith('#')):
          l = re.sub('\s+$', '', l)
          header = l.split('\t')
          samples = [re.sub('\.\S+\.bam$', '', x) for x in header[-meta.shape[0]:]]
          print(samples)
          continue
        
        l = re.sub('\s+$', '', l)
        a = l.split('\t')
        if (a[4]=='<*>'): continue
        
        ref = a[0]
        pos = a[1]
        n = '{}.p{}'.format(ref, pos)
        AD = re.sub('\S+;AD=(\S+);I16\S+$', '\\1', a[7])
        b = [int(x) for x in AD.split(',')]
        idx = b.index(max(b))
        for i,d in enumerate(samples):
          col = a[-(i+1)].split(':')[-1]
          s = samples[-(i+1)]
          b = [int(x) for x in col.split(',')]
          r = b.pop(idx) #get ref allele
          m = sum([int(k) for k in b]) #combine all other alleles
          #print(n, s, r, m)
          self.counts[s+'.rDA'][n] = r
          self.counts[s+'.mDA'][n] = m
    
    log.write("Parsing pileup file and writing to pileup.parsed.txt\n\n")
    if self.rm: 
      count_df = pd.DataFrame.from_dict(self.counts)
      count_df.to_csv('pileup.parsed.txt',sep="\t", header=True, index=True)

  def logit(self, metadata, covariate, outfile):
      log = self.log
      log.write("####Performing logistic regression for mutational analysis#####\n")
      cmd = 'Rscript {}/logit_mut_analysis.R {} {} {} {}'.format(self.matRdir, 'pileup.parsed.txt', metadata, covariate, outfile)
      print(cmd)
      if self.rm: subprocess.call(cmd,shell=True)
      log.write("\n\n")
