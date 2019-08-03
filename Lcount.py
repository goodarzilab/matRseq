import os
import pandas as pd
from ast import literal_eval
import numpy as np
from collections import defaultdict
class merger():

    def __init__(self):
        self.genes = {}                ## for tRNA count in all samples
        self.muts = {}                # for tRNA mutations
        self.counts=defaultdict(dict)
    def add_gene(self, gene, pos):

        if not (gene in self.genes.keys() ):
            self.genes[gene] = []
        if not pos in self.genes[gene]:
            self.genes[gene].append(pos)

    def add_muts(self, gene, sample, pos, cnt, totalcnt):
        if not (sample in self.muts.keys() ):
            self.muts[sample] = {}
        if not (gene in self.muts[sample]):
            self.muts[sample][gene]={"Count": 0, "Pos":{}}
        ttcnt = totalcnt.strip()
        pos = pos.strip()
        self.muts[sample][gene]["Count"] = ttcnt
        self.muts[sample][gene]["Pos"][pos] = cnt

    def make_mutation_tables(self,file):
        with open (file, "r") as myfile:
            dat = myfile.readlines()[1:]
        current_gene = ""
        cnt = 0
        for i in range(len(dat)):
            data = dat[i].strip().split('\t')
            pos = data[0]
            dicti = data[1]

            if(data[0][0] == 't'):
                current_gene = data[0].strip()
                totalcnt = data[1]
            else:
                cnt = literal_eval(dicti[1:])["count"]
                self.add_gene(current_gene,pos)
                self.add_muts(current_gene, file, pos, cnt, totalcnt)

            
    def make_H_file(self):
        positions = []
        for tRNA , poses in self.genes.items():
            for p in sorted(poses):
                n = tRNA.strip()+".p"+p
                positions.append(n.rstrip())

        for i,j in self.muts.items():
            samp = i[:-8]
            for g in sorted(positions):
                t = g.split(".")[0]
                p = g.split(".")[1].strip()[1:]
                if t in list(j.keys()):
                    for tR, v in j.items():
                        if tR == t :
                            if p in list(v["Pos"].keys()):
                                self.counts[samp+".rDA"][g] = str(int(v["Count"]) - int(v["Pos"][p]))
                                self.counts[samp+".mDA"][g] = (v["Pos"][p])
                            else:
                                self.counts[samp+'.rDA'][g] = "0"
                                self.counts[samp+'.mDA'][g] = "0"
                else:
                    self.counts[samp+'.rDA'][g] = "0"
                    self.counts[samp+'.mDA'][g] = "0"


    def export_mutations_table(self):

        count_df = pd.DataFrame.from_dict(self.counts)
        count_df.to_csv('countfile.txt',sep="\t", header=True, index=True)


def main():
    my_merger = merger()
    files = os.listdir()
    cntfiles = sorted(filter(lambda x: x[-8:] == ".tRNAcnt", files))
    for f in cntfiles:
        my_merger.make_mutation_tables(f)
    my_merger.make_H_file()
    my_merger.export_mutations_table()

if __name__ == "__main__":
    main()
