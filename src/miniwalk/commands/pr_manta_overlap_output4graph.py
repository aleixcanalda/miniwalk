import pandas as pd
import sys
import gzip
import re
import os, glob
import argparse
import re
import warnings
import numpy as np
from Bio.Seq import Seq
from Bio import Align
import math

def clean_sequence(seq):
    valid_chars = set("ACGTacgt")  # Assuming DNA sequences
    return ''.join([c for c in seq if c in valid_chars])

def main(options):

    warnings.simplefilter(action='ignore', category=FutureWarning)

    aligner = Align.PairwiseAligner()

    call = pd.read_csv(options.call, sep='\t', header=None,comment='#')
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Sample']
    call.columns = header[:len(call.columns)]

    vcf = pd.read_csv(options.vcf, sep='\t', header=None,comment='#')
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Ref','Sam']
    vcf.columns = header[:len(vcf.columns)]

    rep = pd.read_csv(options.repeat, sep='\t', header=None)
    rep = rep.dropna(axis=1,how='all')

    fp = 0
    fn = 0
    tp = 0
    h37rv_seq = ""
    with open(options.ref,'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                h37rv_seq += line
    vcf_tp = []
    call_fp = []
    call_fn = []
    call_span = []

    for c in call.index:
        if call.iloc[c,7].split(';')[1].split('=')[0] != 'SVTYPE' or call.iloc[c,7].split(';')[2].split('=')[0] != 'SVLEN':
            continue
        elif int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-','')) < 50:
            continue 
        if call.iloc[c,1] in call_span:
            continue
        if c !=1:
            if call.iloc[c,1] == call.iloc[c-1,1]:
                continue
        sv = False
        repe = False
        for r in rep.index:
            if call.iloc[c,1] >= rep.iloc[r,0]-50 and call.iloc[c,1] <= rep.iloc[r,1] + rep.iloc[r,2] + 50:
                repe=True
                rep_pos_start = rep.iloc[r,0] - 50
                rep_pos_end = rep.iloc[r,1] + rep.iloc[r,2] + 50
                break
        for v in vcf.index:
            if isinstance(vcf.iloc[v,2],float) and math.isnan(vcf.iloc[v,2]):
                continue
            P1 = call.iloc[c,1] + int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))
            F2 = call.iloc[c,1]
            F1 = vcf.iloc[v,1] + int(vcf.iloc[v,2].split('.')[1])
            P2 = vcf.iloc[v,1]
            ref_len = int(vcf.iloc[v,2].split('.')[1])*0.25
            sam_len = int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))*0.25
            same=False
            same2=False
            if call.iloc[c,7].split(';')[1].split('=')[1] == 'DUP' and vcf.iloc[v,2].split('.')[0] == 'INS':
                same=True
            if v+1 != len(vcf.index):
                if P1 > vcf.iloc[v+1,1] + int(vcf.iloc[v+1,2].split('.')[1]):
                    P1 = vcf.iloc[v+1,1] + int(vcf.iloc[v+1,2].split('.')[1])
                if call.iloc[c,7].split(';')[1].split('=')[1] == 'DUP' and vcf.iloc[v+1,2].split('.')[0] == 'INS':
                    same2=True
                if P2 > F2:
                    F2 = P2
                if F1 > F2 and P1 >= vcf.iloc[v+1,1] and F1 - F2 >= ref_len and P1 - vcf.iloc[v+1,1] >= int(vcf.iloc[v+1,2].split('.')[1])*0.25 and (F1 - F2) + (P1 - vcf.iloc[v+1,1]) >= sam_len and ((call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0] and call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v+1,2].split('.')[0]) or (same == True and same2==True) or (call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0] and same2==True) or (call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v+1,2].split('.')[0] and same == True)):
                    tp += 1
                    vcf_tp.append(vcf.iloc[v,1])
                    vcf_tp.append(vcf.iloc[v+1,1])
                    sv = True
                    if call.iloc[c,7].split(';')[1].split('=')[1] == "DEL":
                        print("%s,-%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                    else:
                        print("%s,%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                    break
            if c+1 != len(call.index) and call.iloc[c+1,7].split(';')[1].split('=')[0] == 'SVTYPE' and call.iloc[c+1,7].split(';')[2].split('=')[0] == 'SVLEN':
                if int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-','')) > 50 or call.iloc[c+1,7].split(';')[0].split('=')[1] != 'BND':
                    if call.iloc[c+1,7].split(';')[1].split('=')[1] == 'DUP' and vcf.iloc[v,2].split('.')[0] == 'INS':
                        same2=True
                    if F1 > call.iloc[c+1,1] + int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-','')):
                        F1 = call.iloc[c+1,1] + int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-',''))
                    if F2 > P2:
                        P2 = F2
                    if P1 > P2 and F1 >= call.iloc[c+1,1] and P1 - P2 >= sam_len and F1 - call.iloc[c+1,1] >= int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-',''))*0.25 and (P1 - P2) + (F1 - call.iloc[c+1,1]) >= ref_len and ((vcf.iloc[v,2].split('.')[0] == call.iloc[c,7].split(';')[1].split('=')[1] and call.iloc[c+1,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0]) or (same == True and same2 ==True) or (vcf.iloc[v,2].split('.')[0] == call.iloc[c,7].split(';')[1].split('=')[1] and same2 ==True) or (same == True and call.iloc[c+1,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0])):
                        tp += 2
                        vcf_tp.append(vcf.iloc[v,1])
                        sv = True
                        call_span.append(call.iloc[c+1,1])
                        if call.iloc[c,7].split(';')[1].split('=')[1] == "DEL":
                            print("%s,-%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                        else:
                            print("%s,%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                        if call.iloc[c+1,7].split(';')[1].split('=')[1] == "DEL":
                            print("%s,-%i,TP"%(call.iloc[c+1,7].split(';')[1].split('=')[1],int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-',''))))
                        else:
                            print("%s,%i,TP"%(call.iloc[c+1,7].split(';')[1].split('=')[1],int(call.iloc[c+1,7].split(';')[2].split('=')[1].replace('-',''))))
                        break
            if vcf.iloc[v,1] >= call.iloc[c,1] and P1 > F1:
                P1 = F1
            elif call.iloc[c,1] >= vcf.iloc[v,1] and F1 > P1:
                F1 = P1
            if ((vcf.iloc[v,1] >= call.iloc[c,1] and P1 - P2 >= sam_len and P1 - P2 >= ref_len) or (call.iloc[c,1] >= vcf.iloc[v,1] and F1 - F2 >= sam_len and F1 - F2 >= ref_len)) and (call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0] or same ==True):
                tp += 1
                vcf_tp.append(vcf.iloc[v,1])
                sv = True
                if call.iloc[c,7].split(';')[1].split('=')[1] == "DEL":
                    print("%s,-%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                else:
                    print("%s,%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                break
            if repe==True: #if our call and std SVs arent overlapping BUT they are in a TR region, this may still be the same SV, we will compare their adjacent sequences to make sure
                if (call.iloc[c,7].split(';')[1].split('=')[1] == vcf.iloc[v,2].split('.')[0] or same == True) and vcf.iloc[v,1] >= rep_pos_start and vcf.iloc[v,1] <= rep_pos_end and int(aligner.score(clean_sequence(h37rv_seq[call.iloc[c,1]-50:call.iloc[c,1]+int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))+50]),clean_sequence(h37rv_seq[vcf.iloc[v,1]-50:vcf.iloc[v,1]+50+int(vcf.iloc[v,2].split('.')[1])]))) >= int((min(int(vcf.iloc[v,2].split('.')[1]),int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-','')))+100)*0.5):
                    tp += 1
                    vcf_tp.append(vcf.iloc[v,1])
                    sv = True
                    if call.iloc[c,7].split(';')[1].split('=')[1] == "DEL":
                        print("%s,-%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                    else:
                        print("%s,%i,TP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
                    break
        if sv == False:
            fp += 1
            call_fp.append(call.iloc[c,1])
            if call.iloc[c,7].split(';')[1].split('=')[1] == "DEL":
                print("%s,-%i,FP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
            else:
                print("%s,%i,FP"%(call.iloc[c,7].split(';')[1].split('=')[1],int(call.iloc[c,7].split(';')[2].split('=')[1].replace('-',''))))
    for v in vcf.index:
        if isinstance(vcf.iloc[v,2],float) and math.isnan(vcf.iloc[v,2]):
            continue
        if vcf.iloc[v,1] in vcf_tp:
            continue
        elif int(vcf.iloc[v,2].split('.')[1]) < 50:
            continue
        else:
            fn += 1
            call_fn.append(vcf.iloc[v,1])
            if vcf.iloc[v,2].split('.')[0] == "DEL":
                print("%s,-%i,FN"%(vcf.iloc[v,2].split('.')[0],int(vcf.iloc[v,2].split('.')[1])))
            else:
                print("%s,%i,FN"%(vcf.iloc[v,2].split('.')[0],int(vcf.iloc[v,2].split('.')[1])))

