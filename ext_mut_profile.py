import re
import os
import subprocess
import sys
class Mut_Profile():


    def __init__(self):
        self.genes = {}
        
    def add_gene(self, name):
        self.genes[name] = {"miss_matched": {}, "count": 0 }

    def mut(self, gene_name, read, start_index, mutation):
        pattern = re.compile('\d+|[\^CGTA]+')
        tokens = pattern.findall(mutation)
        index = start_index
        missed = self.genes[gene_name]["miss_matched"]

        self.genes[gene_name]["count"] += 1
        for token in tokens:
            delete_flag = False
            if '0' <= token[0] <= '9':
                index += int(token)
                continue

            if token[0] == '^':
                delete_flag = True
                token = token[1:]
            
            for nt in token:
                
                if index not in missed:
                    missed[index] = {"count": 0}
                try:
                    read_nt = read[index - start_index]
                except:
                    pass
                if delete_flag:
                    read_nt = "*"
                diff = (nt, read_nt)
                if diff not in missed[index]:
                    missed[index][diff] = [0, 0, 0]
                missed[index][diff][0] += 1
                # missed[index]["count"] += 1
                if not delete_flag:
                    missed[index]["count"] += 1
                index += 1
    def set_ratio(self):
        for gene in self.genes:
            missed = self.genes[gene]["miss_matched"]
            for index in missed:
                for diff in missed[index]:
                    if diff == "count":
                        continue
                    missed[index][diff][1] = '{0:.10f}'.format((missed[index][diff][0] / missed[index]["count"]))
                    missed[index][diff][2] = '{0:.10f}'.format(missed[index][diff][0] / self.genes[gene]["count"])

    def read_from_file(self, file):
        file = file
        for line in file:
            data = line.strip().split('\t')
            if len(data) == 3:
                self.add_gene(data[1].split(':')[1])
            elif '@' not in data[0] and data[2] != '*':
                gene_name = data[2]
                read = data[9]
                start_index = int(data[3])
                mutation = data[12].split(':')[2]
                self.mut(gene_name, read, start_index, mutation)

    def export(self, outfile=False):
        if (not outfile):
            out=sys.stdout
        for gene in self.genes:
            #print(gene, "\t", self.genes[gene]["count"])
            out.write(gene, "\t", self.genes[gene]["count"], "\n")
            for index in self.genes[gene]["miss_matched"]:
                if self.genes[gene]["miss_matched"][index]["count"] ==0:
                    continue
                #print(index,"\t", self.genes[gene]["miss_matched"][index])  #{('G', '-'): [1, 0.5, 0.1111111111111111], 'count': 2, ('G', 'C'): [1, 0.5, 0.1111111111111111]}
                out.write(index,"\t", self.genes[gene]["miss_matched"][index], "\n")  #{('G', '-'): [1, 0.5, 0.1111111111111111], 'count': 2, ('G', 'C'): [1, 0.5, 0.1111111111111111]}

def main():
    MutP = Mut_Profile()
    a = sys.stdin
    MutP.read_from_file(a)
    MutP.export()

if __name__ == "__main__":
    main()
