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
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Ref', 'Sam']
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
        if isinstance(call.iloc[c,2],float) and math.isnan(call.iloc[c,2]):
            continue
        if call.iloc[c,1] in call_span:
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
            P1 = call.iloc[c,1] + int(call.iloc[c,2].split('.')[1])
            F2 = call.iloc[c,1]
            F1 = vcf.iloc[v,1] + int(vcf.iloc[v,2].split('.')[1])
            P2 = vcf.iloc[v,1]
            ref_len = int(vcf.iloc[v,2].split('.')[1])*0.25
            sam_len = int(call.iloc[c,2].split('.')[1])*0.25
            if v+1 != len(vcf.index):
                if P1 > vcf.iloc[v+1,1] + int(vcf.iloc[v+1,2].split('.')[1]):
                    P1 = vcf.iloc[v+1,1] + int(vcf.iloc[v+1,2].split('.')[1])
                if P2 > F2:
                    F2 = P2
                if F1 > F2 and P1 >= vcf.iloc[v+1,1] and F1 - F2 >= ref_len and P1 - vcf.iloc[v+1,1] >= int(vcf.iloc[v+1,2].split('.')[1])*0.25 and (F1 - F2) + (P1 - vcf.iloc[v+1,1]) >= sam_len and call.iloc[c,2].split('.')[0] == vcf.iloc[v,2].split('.')[0] and call.iloc[c,2].split('.')[0] == vcf.iloc[v+1,2].split('.')[0]:
                    tp += 1
                    vcf_tp.append(vcf.iloc[v,1])
                    vcf_tp.append(vcf.iloc[v+1,1])
                    sv = True
                    break
            if c+1 != len(call.index):
                if F1 > call.iloc[c+1,1] + int(call.iloc[c+1,2].split('.')[1]):
                    F1 = call.iloc[c+1,1] + int(call.iloc[c+1,2].split('.')[1])
                if F2 > P2:
                    P2 = F2
                if P1 > P2 and F1 >= call.iloc[c+1,1] and P1 - P2 >= sam_len and F1 - call.iloc[c+1,1] >= int(call.iloc[c+1,2].split('.')[1])*0.25 and (P1 - P2) + (F1 - call.iloc[c+1,1]) >= ref_len and vcf.iloc[v,2].split('.')[0] == call.iloc[c,2].split('.')[0] and call.iloc[c+1,2].split('.')[0] == vcf.iloc[v,2].split('.')[0]:
                    tp += 2
                    vcf_tp.append(vcf.iloc[v,1])
                    sv = True
                    call_span.append(call.iloc[c+1,1])
                    break
            if vcf.iloc[v,1] >= call.iloc[c,1] and P1 > F1:
                P1 = F1
            elif call.iloc[c,1] >= vcf.iloc[v,1] and F1 > P1:
                F1 = P1
            if ((vcf.iloc[v,1] >= call.iloc[c,1] and P1 - P2 >= sam_len and P1 - P2 >= ref_len) or (call.iloc[c,1] >= vcf.iloc[v,1] and F1 - F2 >= sam_len and F1 - F2 >= ref_len)) and call.iloc[c,2].split('.')[0] == vcf.iloc[v,2].split('.')[0]:
                tp += 1
                vcf_tp.append(vcf.iloc[v,1])
                sv = True
                break
            if repe==True: #if our call and std SVs arent overlapping BUT they are in a TR region, this may still be the same SV, we will compare their adjacent sequences to make sure
                if call.iloc[c,2].split('.')[0] == vcf.iloc[v,2].split('.')[0] and vcf.iloc[v,1] >= rep_pos_start and vcf.iloc[v,1] <= rep_pos_end and int(aligner.score(clean_sequence(h37rv_seq[call.iloc[c,1]-50:call.iloc[c,1]+int(call.iloc[c,2].split('.')[1])+50]),clean_sequence(h37rv_seq[vcf.iloc[v,1]-50:vcf.iloc[v,1]+50+int(vcf.iloc[v,2].split('.')[1])]))) >= int((min(int(vcf.iloc[v,2].split('.')[1]),int(call.iloc[c,2].split('.')[1]))+100)*0.5):
                    tp += 1
                    vcf_tp.append(vcf.iloc[v,1])
                    sv = True
                    break
        if sv == False:
            fp += 1
            call_fp.append(call.iloc[c,1])
    for v in vcf.index:
        if isinstance(vcf.iloc[v,2],float) and math.isnan(vcf.iloc[v,2]):
            continue
        if vcf.iloc[v,1] in vcf_tp:
            continue
        else:
            fn += 1
            call_fn.append(vcf.iloc[v,1])

    precision = tp/(tp+fp)
    recall = len(vcf_tp)/(len(vcf_tp)+fn)

    print("%f\tMTB-Graph\tPrecision\n%f\tMTB-Graph\tRecall\n%f\tMTB-Graph\tF1-score" %(precision,recall,2*precision*recall/(precision+recall)))