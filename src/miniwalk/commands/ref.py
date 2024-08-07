import pandas as pd
import sys
import gzip
import re
import os, glob
import argparse
import re
import warnings
from Bio.Seq import Seq
from Bio import Align
import math

def main(options):

    def print_vcf(vcf,n,node=None):
        """Print SV as a new vcf line"""
        if node:
            print(node)
        with open(options.out,'a') as o:
            print(vcf.loc[[n]].to_string(index=False,header=False).replace(' ','\t'),file=o)

    warnings.simplefilter(action='ignore', category=FutureWarning)

    aligner = Align.PairwiseAligner()
    aligner.mode='local'
    aligner.open_gap_score=-2
    aligner.extend_gap_score=-0.5

    vcf = pd.read_csv(options.vcf, sep='\t', header=None,comment='#')
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Ref', 'Sam']
    vcf.columns = header[:len(vcf.columns)]

    with open(options.out,'w') as o:
        print("##fileformat=VCFv4.2\n##ALT=<ID=CNV,Description=\"description\">\n##contig=<ID=NC_000962.3,length=4411532>\n##FORMAT=<ID=GT0,Number=1,Type=String,Description=\"Original genotype\">\n##ALT=<ID=X:1,Description=\"Allele 1\">\n##ALT=<ID=X:2,Description=\"Allele 2\">\n##ALT=<ID=X:3,Description=\"Allele 3\">\n##ALT=<ID=X:4,Description=\"Allele 4\">\n##ALT=<ID=X:5,Description=\"Allele 5\">\n##ALT=<ID=X:6,Description=\"Allele 6\">\n##ALT=<ID=X:7,Description=\"Allele 7\">\n##ALT=<ID=X:8,Description=\"Allele 8\">\n##ALT=<ID=X:9,Description=\"Allele 9\">\n##ALT=<ID=X:10,Description=\"Allele 10\">\n##ALT=<ID=X:11,Description=\"Allele 11\">\n##ALT=<ID=X:12,Description=\"Allele 12\">\n##ALT=<ID=X:13,Description=\"Allele 13\">\n##ALT=<ID=X:14,Description=\"Allele 14\">\n##ALT=<ID=X:15,Description=\"Allele 15\">\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Number of alleles\">\n##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count\">\n##INFO=<ID=ALEN,Number=.,Type=Integer,Description=\"Length of each allele\">\n##INFO=<ID=ANNO,Number=1,Type=String,Description=\"Annotation\">\n##INFO=<ID=VS,Number=1,Type=String,Description=\"Start vertex\">\n##INFO=<ID=VE,Number=1,Type=String,Description=\"End vertex\">\n##INFO=<ID=AWALK,Number=.,Type=String,Description=\"Walk of each allele\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=CSTRAND,Number=1,Type=String,Description=\"Contig strand\">\n##FORMAT=<ID=CTG,Number=1,Type=String,Description=\"Contig name\">\n##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Contig start, BED-like\">\n##FORMAT=<ID=CE,Number=1,Type=String,Description=\"Contig end, BED-like\">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NC_000962.3 %s" %(options.vcf),file=o)

    vcf['Prefix'] = vcf['ID'].str.split('.').str[0]
    # Assign sort order based on 'Prefix' ('INS' comes before 'DEL')
    vcf['Sort_Order'] = vcf['Prefix'].apply(lambda x: 0 if x == 'DEL' else 1)

    # Sort DataFrame by 'POS' in ascending order and 'Sort_Order' in ascending order
    vcf = vcf.sort_values(by=['POS', 'Sort_Order']).reset_index().drop('index',axis=1)

    # Drop the temporary columns 'Prefix' and 'Sort_Order'
    vcf = vcf.drop(columns=['Prefix', 'Sort_Order'])

    #vcf = vcf.sort_values(by=['ID'],ascending=['False'],key=lambda x: x.str.lower()).reset_index().drop('index',axis=1)
    n_cluster=[]

    for n in vcf.index:
        if isinstance(vcf.iloc[n,2],float) and math.isnan(vcf.iloc[n,2]):
            if n+1 == len(vcf.index):
                sv_end = vcf.iloc[n, 1]
            elif isinstance(vcf.iloc[n+1,2],float) and math.isnan(vcf.iloc[n+1,2]):
                sv_end = vcf.iloc[n+1, 1]
            else:
                sv_end = vcf.iloc[n+1, 1] + int(vcf.iloc[n+1, 2].split('.')[1]) + 150 
            print_vcf(vcf,n)
            continue
        if vcf.iloc[n,2].split('.')[0] == "INS":
            ins=1
        if vcf.iloc[n,1] in n_cluster:
            if n == len(vcf.index)-1:
                n = n - len(n_cluster)
                if len(vcf.iloc[n,3]) > len(vcf.iloc[n,4]):
                    if vcf.iloc[n,4] == '*':
                        vcf.iloc[n, 2] = "DEL.%i" % (len(vcf.iloc[n, 3]))
                    else:
                        vcf.iloc[n,2] = "DEL.%i" %(len(vcf.iloc[n,3])-len(vcf.iloc[n,4]))
                        if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                            vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                else:
                    if vcf.iloc[n, 3] == '*':
                        vcf.iloc[n, 2] = "INS.%i" % (len(vcf.iloc[n, 4]))
                    else:
                        vcf.iloc[n,2] = "INS.%i" %(len(vcf.iloc[n,4])-len(vcf.iloc[n,3]))
                            #print(vcf.iloc[n,1],aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0],aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0])
                        if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                            vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                if int(vcf.iloc[n,2].split('.')[1]) > 49:
                    print_vcf(vcf, n)
            continue
        if vcf.iloc[n,2].split('.')[0] == 'INV':
            sv_end = vcf.iloc[n+1, 1] + int(vcf.iloc[n+1, 2].split('.')[1]) + 150 
            print_vcf(vcf,n)
            continue
        if n < 15:
            sv_end = vcf.iloc[n,1] + int(vcf.iloc[n,2].split('.')[1]) + 150 
        if n == len(vcf.index)-1:
            print_vcf(vcf,n)
        n_cluster=[]
        for x in vcf.iloc[n+1:,].index:
            if isinstance(vcf.iloc[x,2],float) and math.isnan(vcf.iloc[x,2]):
                if n_cluster != []:
                    if len(vcf.iloc[n, 3]) > len(vcf.iloc[n, 4]):
                        if vcf.iloc[n, 4] == '*':
                            vcf.iloc[n, 2] = "DEL.%i" % (len(vcf.iloc[n, 3]))
                        else:
                            vcf.iloc[n, 2] = "DEL.%i" % (
                                len(vcf.iloc[n, 3])-len(vcf.iloc[n, 4]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))    
                    else:
                        if vcf.iloc[n, 3] == '*':
                            vcf.iloc[n, 2] = "INS.%i" % (len(vcf.iloc[n, 4]))
                        else:
                            vcf.iloc[n, 2] = "INS.%i" % (
                                len(vcf.iloc[n, 4])-len(vcf.iloc[n, 3]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                    if n+1 == len(vcf.index)-1:
                        sv_end = vcf.iloc[n, 1]
                    elif isinstance(vcf.iloc[n+1,2],float) and math.isnan(vcf.iloc[n+1,2]):
                        sv_end = vcf.iloc[n+1, 1]
                    else:
                        sv_end = vcf.iloc[n+1, 1] + int(vcf.iloc[n+1, 2].split('.')[1]) + 150 
                else:
                    if n+1 == len(vcf.index):
                        sv_end = vcf.iloc[n, 1]
                    elif isinstance(vcf.iloc[n+1,2],float) and math.isnan(vcf.iloc[n+1,2]):
                        sv_end = vcf.iloc[n+1, 1]
                    else:
                        sv_end = vcf.iloc[n+1, 1] + int(vcf.iloc[n+1, 2].split('.')[1]) + 150 
                if int(vcf.iloc[n, 2].split('.')[1]) > 49:
                    print_vcf(vcf, n)
                break
            if vcf.iloc[x, 2].split('.')[0] == 'INV':
                if n_cluster != []:
                    if len(vcf.iloc[n, 3]) > len(vcf.iloc[n, 4]):
                        if vcf.iloc[n, 4] == '*':
                            vcf.iloc[n, 2] = "DEL.%i" % (len(vcf.iloc[n, 3]))
                        else:
                            vcf.iloc[n, 2] = "DEL.%i" % (
                                len(vcf.iloc[n, 3])-len(vcf.iloc[n, 4]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                    else:
                        if vcf.iloc[n, 3] == '*':
                            vcf.iloc[n, 2] = "INS.%i" % (len(vcf.iloc[n, 4]))
                        else:
                            vcf.iloc[n, 2] = "INS.%i" % (
                                len(vcf.iloc[n, 4])-len(vcf.iloc[n, 3]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                    sv_end = vcf.iloc[x, 1] + \
                        int(vcf.iloc[x, 2].split('.')[1]) + 150 
                else:
                    sv_end = vcf.iloc[n, 1] + \
                        int(vcf.iloc[n, 2].split('.')[1]) + 150 
                if int(vcf.iloc[n, 2].split('.')[1]) > 49:
                    print_vcf(vcf, n)
                break
            if vcf.iloc[x,1] > sv_end:
                if n_cluster != []:
                    if len(vcf.iloc[n,3]) > len(vcf.iloc[n,4]):
                        if vcf.iloc[n,4] == '*':
                            vcf.iloc[n, 2] = "DEL.%i" % (len(vcf.iloc[n, 3]))
                        else:
                            vcf.iloc[n,2] = "DEL.%i" %(len(vcf.iloc[n,3])-len(vcf.iloc[n,4]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                    else:
                        if vcf.iloc[n, 3] == '*':
                            vcf.iloc[n, 2] = "INS.%i" % (len(vcf.iloc[n, 4]))
                        else:
                            vcf.iloc[n,2] = "INS.%i" %(len(vcf.iloc[n,4])-len(vcf.iloc[n,3]))
                            if int(aligner.score(vcf.iloc[n, 3],vcf.iloc[n, 4])) >= min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))*0.85 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[0][0][0] < 20 and aligner.align(vcf.iloc[n, 3],vcf.iloc[n, 4])[0].aligned[1][0][0] < 15:
                                vcf.iloc[n,1]+=min(len(vcf.iloc[n, 3]),len(vcf.iloc[n, 4]))
                    if x >= len(vcf)-2:
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif isinstance(vcf.iloc[x+2,2],float) and math.isnan(vcf.iloc[x+2,2]):
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif isinstance(vcf.iloc[x+1,2],float) and math.isnan(vcf.iloc[x+1,2]):
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif vcf.iloc[x+2,2].split('.')[0] == "INS" and vcf.iloc[x+1,2].split('.')[0] == "INS":
                        sv_end = vcf.iloc[x,1] + 250
                    else:
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                else:
                    if x >= len(vcf)-2:
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif isinstance(vcf.iloc[x+2,2],float) and math.isnan(vcf.iloc[x+2,2]):
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif isinstance(vcf.iloc[x+1,2],float) and math.isnan(vcf.iloc[x+1,2]):
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                    elif vcf.iloc[x+2,2].split('.')[0] == "INS" and vcf.iloc[x+1,2].split('.')[0] == "INS":
                        sv_end = vcf.iloc[x,1] + 250
                    else:
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                if int(vcf.iloc[n,2].split('.')[1]) > 49:
                    print_vcf(vcf, n)
                break
            else:
                n_cluster.append(vcf.iloc[x,1])
                #if vcf.iloc[n,2] == vcf.iloc[x,2]:
                #    continue
                if vcf.iloc[x,2].split('.')[0] == "DEL":
                    if vcf.iloc[n,3] == '*':
                        vcf.iloc[n, 3] = vcf.iloc[x, 3]
                    else:
                        vcf.iloc[n,3] += vcf.iloc[x,3]
                else:
                    if vcf.iloc[n, 4] == '*':
                        vcf.iloc[n, 4] = vcf.iloc[x, 4]
                    else:
                        vcf.iloc[n, 4] += vcf.iloc[x, 4]
                if x >= len(vcf)-1:
                    continue
                elif isinstance(vcf.iloc[x+1,2],float) and math.isnan(vcf.iloc[x+1,2]):
                    if (sv_end < (vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]))):
                        sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
                elif vcf.iloc[n,2].split('.')[0]=="INS" and vcf.iloc[x+1,2].split('.')[0]=="INS":
                    sv_end = vcf.iloc[x,1] + 250
                elif (sv_end < (vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]))):
                    sv_end = vcf.iloc[x,1] + int(vcf.iloc[x,2].split('.')[1]) + 150 
    
