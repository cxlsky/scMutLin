'''
   split10x: split 10x bam file into several sub-bam files based their CB tag
       input:
            -i <inputfile>, the bam file from 10x experiment
            -b <barcodefile>,the cell barcode file that is used to filter output [Optional]
       output:
            -o <outputpath>, sub-bam files the named as the CB tag
   author: Xiaolong from Wei Li Lab
   email:xcheng@childrensnational.org
'''
import os
from sys import argv
import pysam
import getopt
import pandas as pd

# parse parameters
dir_out = './'
data_bam = ''
data_bc = ''
flag_limit = False
opts, args = getopt.getopt(argv[1:],"i:o:b:")
for opt, arg in opts:
    if opt == "-i":
        data_bam = arg
    elif opt == "-o":
        dir_out = arg
    elif opt == "-b":
        data_bc = arg

# cell barcode used in scRNA-seq analysis
cbc = []
if data_bc:
    df = pd.read_csv(data_bc,  header=None, names=['cbc'], index_col=None)
    cbc = [i[0].split('-')[0] for i in df.values]
    cbc = list(set(cbc))

# open bam file
fin = pysam.AlignmentFile(data_bam, "rb")

# open the number of bam files and map the out file handler to the CB id, write to a bam with wb
fouts_dict = {}
for read in fin:
    tags = read.tags
    cell_barcode = ''
    for x in tags:
        if x[0] == 'CB':
            cell_barcode = x[1].split('-')[0]
            if len(cbc) > 0:
                if cell_barcode not in cbc:
                    break
    if len(cell_barcode) == 0:
        continue
    if cell_barcode in fouts_dict.keys():
        fouts_dict[cell_barcode].write(read)
    elif cell_barcode in cbc:
        fouts_dict[cell_barcode] = pysam.AlignmentFile(os.path.join(dir_out, cell_barcode+".bam"), "wb", template = fin)
        fouts_dict[cell_barcode].write(read)

fin.close()
for fout in fouts_dict.values():
    fout.close()

