import pandas as pd
import os, glob
import argparse
import warnings
import numpy as np
import subprocess
import math

def main(options):
    warnings.simplefilter(action='ignore', category=FutureWarning)
    call = pd.read_csv(options.call, sep='\t', header=None,comment='#')
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Ref', 'Sam']
    call.columns = header[:len(call.columns)]
    
    h37rv_seq = ""
    with open(options.ref,'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                h37rv_seq += line.rstrip()
    
    for c in call.index:
        if isinstance(call.iloc[c,2],float) and math.isnan(call.iloc[c,2]):
            continue
        if call.iloc[c,2].split('.')[0] == 'INS':
            with open(f"{options.call}.fa",'w') as f:
                print(f">Seq\n{h37rv_seq[call.iloc[c,1]-len(call.iloc[c,4])*2:call.iloc[c,1]]+call.iloc[c,4]+h37rv_seq[call.iloc[c,1]:call.iloc[c,1]+len(call.iloc[c,4])*2]}",file=f)
            args2=['repeat-match','-t',f'{options.call}.fa']
            p=subprocess.Popen(args2,stdout=subprocess.PIPE)
            stdout=p.stdout.read()
            try: #We only investigate further those INS instances that have TRs according to repeat-match
                lengths=[]
                stdout.decode('utf-8')[57]
                for n in range(78,len(stdout.decode('utf-8')),31):
                    try:
                        lengths.append(int(stdout.decode('utf-8')[n]+stdout.decode('utf-8')[n+1]+stdout.decode('utf-8')[n+2]))
                    except:
                        break
                for x in lengths: #Only those TRs that are minimum half the size of the INS will be considered as real DUP
                    if x >= int(call.iloc[c,2].split('.')[1])*0.5:
                        subprocess.call(f"sed -i s/'{call.iloc[c,1]}\tINS'/'{call.iloc[c,1]}\tDUP'/g {options.call}",shell=True)
                        break            
            except IndexError:
                continue
        
