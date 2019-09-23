import subprocess
import re
import os
from ext_mut_profile import Mut_Profile
from Lcount import merger

#overall structure
#NNNN (umi1) ACTGGATAC TGGN (tRNA) GTATCCAGT NNNN (umi2)
#all reads start with NNNNACTGGATACTGGN
#reverse read starts with NNNNACTGGATAC
#R1 will be 26nt and will read UMI1
#R2 will be 98nt and will read UMI2
class process:
    def __init__(self, metadata, matRdir, reffile, log, isPaired=True, useRead="R2", hasUMI=True, aligner="BWA", runMode=True):
        self.metadata = metadata
        self.matRdir = matRdir
        self.reffile = reffile
        self.log = log
        self.paired = isPaired
        self.read = useRead
        self.umi = hasUMI
        self.aligner = aligner
        self.rm = runMode #print mode if set to false

    def umi_extract (self):
        if (not self.umi): return True
        
        meta = self.metadata
        log = self.log
        log.write("####Extracting UMIs#####\n")
        
        if (self.paired):
            for sample in meta.index:
                print(sample)
                r1 = meta.loc[sample, 'R1']
                r2 = meta.loc[sample, 'R2']
                o1 = re.sub(".fastq.gz", ".c.fastq.gz", r1)
                o2 = re.sub(".fastq.gz", ".c.fastq.gz", r2)
                meta.loc[sample, 'R1'] = o1
                meta.loc[sample, 'R2'] = o2
                cmd = 'umi_tools extract --stdin={} --read2-in={} --bc-pattern=NNNN --bc-pattern2=NNNN -L log.out --stdout={} --read2-out={}'.format(r1, r2, o1, o2)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        else:
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                o1 = re.sub(".fastq.gz", ".c.fastq.gz", r1)
                meta.loc[sample, 'R1'] = o1
                cmd = 'umi_tools extract --stdin={} --bc-pattern=NNNN -L log.out --stdout={}'.format(r1, o1)      
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)      
        log.write("\n\n")

    def trim(self):
        meta = self.metadata
        log = self.log
        log.write("####Removing adaptors#####\n")
        
        if (self.paired and self.read=="both"):
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                r2 = meta.loc[sample, 'R2']
                o1 = re.sub(".fastq.gz", ".trim.fastq.gz", r1)
                o2 = re.sub(".fastq.gz", ".trim.fastq.gz", r2)
                cmd = 'cutadapt -e 0.12 -m 20 -q 15 -a ^ACTGGATACTGGN...GTATCCAGT -A ^ACTGGATAC...NCCAGTATCCAGT -o {} -p {} {} {}'.format(o1, o2, r1, r2)
                meta.loc[sample, 'R1'] = o1
                meta.loc[sample, 'R2'] = o2
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        elif (self.paired and self.read=="R2"):
            for sample in meta.index:
                r2 = meta.loc[sample, 'R2']
                o2 = re.sub(".fastq.gz", ".trim.fastq.gz", r2)
                cmd = 'cutadapt -j 8 --trimmed-only  -e 0.12 -m 55 -q 15 -a ^ACTGGATAC...NCCAGTATCCAGT -o {} {}'.format(o2,r2)
                meta.loc[sample, 'R2'] = o2
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        else:
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                o1 = re.sub(".fastq.gz", ".trim.fastq.gz", r1)
                cmd = 'cutadapt -j 8 --trimmed-only  -e 0.12 -m 55 -q 15 -a ^ACTGGATACTGGN...GTATCCAGT -o {} {}'.format(o1,r1)
                meta.loc[sample, 'R1'] = o1
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)           
        log.write("\n\n")       

    def merge(self):
        meta = self.metadata
        log = self.log
        log.write("####Merging reads if required#####\n")
        
        if (self.paired and self.read=="both"):
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                r2 = meta.loc[sample, 'R2']
                cmd = 'pear -n 30 -f {} -r {} -o {} '.format(r1, r2, sample)
                cmd = cmd + ' ; rm *.discarded.fastq; rm *.unassembled.*; '
                cmd = cmd + ' ; mv {}.assembled.fastq {}.post.fastq; gzip {}.post.fastq'.format(sample, sample, sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        elif (self.paired and self.read=="R2"):
            for sample in meta.index:
                r2 = meta.loc[sample, 'R2']
                cmd = 'zcat {} | fastx_reverse_complement -z -i - -o {}.post.fastq.gz'.format(r2,sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        else:
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                cmd = 'mv {} {}.post.fastq.gz'.format(r1, sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)     
        log.write("\n\n")       
            
    def align(self):
        meta = self.metadata
        log = self.log
        for sample in meta.index:
            if (self.aligner=="BWA"):
                cmd = 'bwa mem -t 8 {} {}.post.fastq.gz > {}.sam'.format(self.reffile,sample,sample)
                cmd = cmd + '\nsamtools view -F 0x14 -Sb {}.sam > {}.bam; rm {}.sam; samtools sort -@ 12 -o {}.srt.bam {}.bam; samtools index {}.srt.bam'.format(sample,sample,sample,sample,sample,sample)
            print(cmd)
            log.write("%s\n" % (cmd))
            if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

    def dedup(self):
        if (not self.umi): return True
        
        meta = self.metadata
        log = self.log
        for sample in meta.index:
            cmd = 'umi_tools dedup -I {}.srt.bam --output-stats=deduplicated -S {}.dd.bam'.format(sample,sample)   # sample : SW480_LVM2_ribo_r2.srt.bam
            print(cmd)
            log.write("%s\n" % (cmd))
            if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

    def count(self):
        meta = self.metadata
        log = self.log
        for sample in meta.index:
            if (self.umi):
                cmd = 'samtools view %s.dd.bam | cut -f3 | sort | uniq -c | awk \'{ print $2 "\t" $1}\' > %s.cnt' % (sample,sample)
            else:
                cmd = 'samtools view %s.srt.bam | cut -f3 | sort | uniq -c | awk \'{ print $2 "\t" $1}\' > %s.cnt' % (sample,sample)
            print(cmd)
            log.write("%s\n" % (cmd))
            if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

    def make_sam(self):
        meta = self.metadata
        log = self.log
        for sample in meta.index:
            if (self.umi):
                cmd = 'samtools  view -h -o %s.sam %s.dd.bam' %(sample,sample)
            else:
                cmd = 'samtools  view -h -o %s.sam %s.srt.bam' %(sample,sample)
            print(cmd)
            log.write("%s\n" % (cmd))
            if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

    def make_mut_profile(self):
        meta = self.metadata
        log = self.log
        log.write("####Making mutation profiles from sam files#####\n")
        #cmd = 'for f in *.sam\n' + 'do\n' + 'out=${f/.sam/.tRNAcnt}\n'+ 'python3 ext_mut_profile.py < $f > $out\n' + 'done'
        #if self.rm: subprocess.call(cmd,shell=True)
        for sample in meta.index:
            MutP = Mut_Profile()
            MutP.read_from_file(open('%s.sam'%(sample), "rt"))
            MutP.export(open('%s.tRNAcnt'%(sample), "wt"))
        log.write("\n\n")

    def count_file(self):
        meta = self.metadata
        log = self.log
        log.write("####generating the count matrix using Lcount.py#####\n")
        my_merger = merger()
        cntfiles = [str(x)+".tRNAcnt" for x in meta.index]
        for f in cntfiles:
            my_merger.make_mutation_tables(f)
        my_merger.make_H_file()
        my_merger.export_mutations_table()
        log.write("\n\n")

    def mlogit(self, metadata, covariate, outfile):
      log = self.log
      log.write("####Performing logistic regression for mutational analysis#####\n")
      cmd = 'Rscript {}/logit_mut_analysis.R {} {} {} {}'.format(self.matRdir, 'countfile.txt', metadata, covariate, outfile)
      print(cmd)
      if self.rm: subprocess.call(cmd,shell=True)
      log.write("\n\n")

    def univariate(self, metadata, design, ref, outfile):
        log = self.log
        log.write("####Performing univariate analysis using DESeq2#####\n")
        cmd = 'Rscript {}/univariate_matR_analysis.R {} {} {} {}'.format(self.matRdir, metadata, design, ref, outfile)
        print(cmd)
        if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

    def logit(self, metadata, design, outfile):
        log = self.log
        log.write("####Performing logistic regression for logit analysis#####\n")
        cmd = 'Rscript {}/logit_matR_analysis.R {} {} {}'.format(self.matRdir, metadata, design, outfile)
        print(cmd)
        if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")
