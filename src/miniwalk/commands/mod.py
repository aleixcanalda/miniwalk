import pandas as pd
import sys
import gzip
import re
import os, glob
import argparse
import re
from Bio.Seq import Seq
import warnings
from Bio import Align

def main(options):

    def print_vcf(vcf,n,node=None):
        """Print SV as a new vcf line"""
        if node:
            print(node)
        with open(options.out,'a') as o:
            print(vcf.loc[[n]].to_string(index=False,header=False).replace(' ','\t'),file=o)

    aligner = Align.PairwiseAligner()
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    bed = pd.read_csv(options.bed, sep='\t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'startNode', 'endNode', 'path']
    bed.columns = header[:len(bed.columns)]
    
    vcf = pd.read_csv(options.vcf, sep='\t', header=None,comment='#')
    header = ['chrom', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Ref', 'Sam']
    vcf.columns = header[:len(vcf.columns)]
    
    gfa = pd.read_csv(options.gfa, sep='\t', header=None,comment='L',index_col=1)
    header = ['SL', 'sequence', 'length','contig', 'contig_position', 'path']
    gfa.columns = header[:len(gfa.columns)]
    
    with open(options.out,'w') as o:
        print("##fileformat=VCFv4.2\n##contig=<ID=NC_000962.3,length=4411532>\n##ALT=<ID=CNV,Description=\"description\">\n##FORMAT=<ID=GT0,Number=1,Type=String,Description=\"Original genotype\">\n##ALT=<ID=X:1,Description=\"Allele 1\">\n##ALT=<ID=X:2,Description=\"Allele 2\">\n##ALT=<ID=X:3,Description=\"Allele 3\">\n##ALT=<ID=X:4,Description=\"Allele 4\">\n##ALT=<ID=X:5,Description=\"Allele 5\">\n##ALT=<ID=X:6,Description=\"Allele 6\">\n##ALT=<ID=X:7,Description=\"Allele 7\">\n##ALT=<ID=X:8,Description=\"Allele 8\">\n##ALT=<ID=X:9,Description=\"Allele 9\">\n##ALT=<ID=X:10,Description=\"Allele 10\">\n##ALT=<ID=X:11,Description=\"Allele 11\">\n##ALT=<ID=X:12,Description=\"Allele 12\">\n##ALT=<ID=X:13,Description=\"Allele 13\">\n##ALT=<ID=X:14,Description=\"Allele 14\">\n##ALT=<ID=X:15,Description=\"Allele 15\">\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Number of alleles\">\n##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count\">\n##INFO=<ID=ALEN,Number=.,Type=Integer,Description=\"Length of each allele\">\n##INFO=<ID=ANNO,Number=1,Type=String,Description=\"Annotation\">\n##INFO=<ID=VS,Number=1,Type=String,Description=\"Start vertex\">\n##INFO=<ID=VE,Number=1,Type=String,Description=\"End vertex\">\n##INFO=<ID=AWALK,Number=.,Type=String,Description=\"Walk of each allele\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=CSTRAND,Number=1,Type=String,Description=\"Contig strand\">\n##FORMAT=<ID=CTG,Number=1,Type=String,Description=\"Contig name\">\n##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Contig start, BED-like\">\n##FORMAT=<ID=CE,Number=1,Type=String,Description=\"Contig end, BED-like\">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NC_000962.3 %s" %(options.vcf),file=o)
    
    #We get all the bubbles where the pangenome graph was unable to resolve the paths due to its complexity
    nas = [bed.iloc[m,1] for m in bed.index if bed.iloc[m,5] == '.']
     
    
    for n in vcf.index:
        #if the row we're looking at has an NA in the minigraph bed output, we cannot determine the SV, mark it as NA
        if options.na:
            if vcf.iloc[n,1] in nas:
                vcf.iloc[n,2] = "NA"
                vcf.iloc[n,3] = "NA"
                vcf.iloc[n,4] = "NA"
                print_vcf(vcf,n)
                continue
            #If there is no SV, skip line
        if vcf.iloc[n,4] == '.':
            continue#print(vcf.loc[[n]].to_string(index=False,header=False).replace(' ','\t'))
        #otherwise, determine which SV we're looking at
        else:
            #Firstly, we get the ref path and the sample path so we can compare them (they are separated by a comma), we also separate each node, however maintaining those that are reversed "<"
            walk = vcf.iloc[n,7].split(';')[-1].replace('AWALK=','')
            nodes_ref = walk.split(',')[0].split('>')
            nodes_sam = walk.split(',')[1].split('>')
            if nodes_ref[0] == '' :
                nodes_ref = nodes_ref[1:]
            if nodes_sam[0] == '':
                nodes_sam = nodes_sam[1:]
            #Biallelic DEL; if the sample node is empty
            if nodes_sam[0] == '*':
                seq=""
                for node in nodes_ref:
                    #the deletion sequence and length depends on the nodes that make it up, therefore adding them up is necessary, as well as taking into acount those nodes with "<"
                    if re.search('<',node) is not None:
                        comps = node.split('<')
                        for comp in comps:
                            if comp == '':
                                continue
                            if re.search('<',comp) is not None:
                                seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                            else:
                                seq += gfa.loc[comp,'sequence']
                    else:
                        seq += gfa.loc[node,'sequence']
                vcf.iloc[n,2] = "DEL.%i" %(len(seq))
                vcf.iloc[n,3] = seq
                vcf.iloc[n,4] = "*"
                if len(seq) > 49:
                    print_vcf(vcf,n)
            # Biallelic INS
            elif nodes_ref[0] == '*':
                seq=""
                #Now we're going to get the INS sequence. If one of the nodes is on the - strand (has a "<"), we'll get its reverse complement sequence.
                for node in nodes_sam:
                    if re.search('<',node) is not None:
                        comps = node.split('<')
                        for comp in comps:
                            if comp == '':
                                continue
                            if re.search('<',comp) is not None:
                                seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                            else:
                                seq += gfa.loc[comp,'sequence']
                    else:
                        seq += gfa.loc[node,'sequence']                        
                vcf.iloc[n,2] = "INS.%i" %(len(seq))
                vcf.iloc[n,4] = seq
                vcf.iloc[n,3] = "*"
                if len(seq) > 49:
                    print_vcf(vcf,n)
            else:
                #we get all the nodes that are unique to the reference or the sample
                diff_nodes = set(nodes_ref).symmetric_difference(set(nodes_sam))
                del_index = {}
                ins_index = {}
                #we get the nodes in the ref and in the sam, separately, as well as their index
                for diff in diff_nodes:
                    try:
                        del_index[diff] = nodes_ref.index(diff)
                    except ValueError:
                        ins_index[diff] = nodes_sam.index(diff)
                sort_ins_index = sorted(ins_index.items(), key=lambda x:x[1])
                sort_del_index = sorted(del_index.items(), key=lambda x:x[1])
                del_svs = []
                ins_svs = []
                for sort in sort_ins_index:
                    if len(ins_svs) == 0:
                        ins_svs.append([sort[0]])
                        i = sort[1]
                    else:
                        if sort[1] == i + 1:
                            ins_svs[-1].append(sort[0])
                            i = sort[1]
                        else:
                            ins_svs.append([sort[0]])
                            i = sort[1]
                #We sort the unique nodes by their index; if nodes are contiguous to each other, they form the whole SV, otherwise they are different SVs
                for sort in sort_del_index:
                    if len(del_svs) == 0:
                        del_svs.append([sort[0]])
                        i = sort[1]
                    else:
                        if sort[1] == i + 1:
                            del_svs[-1].append(sort[0])
                            i = sort[1]
                        else:
                            del_svs.append([sort[0]])
                            i = sort[1]
                #print(ins_svs,del_svs)
                pos = vcf.iloc[n,1]
                sv_inv = []
                inv_bub=False
                #now we go through each SV unique to the sample (making it an INS, INV or candidate DUP)
                for sv in ins_svs:
                    type_sv="INS"
                    seq=""
                    vcf.iloc[n,1] = pos
                    inv_size=0
                    inv_ref=''
                    inv_alt=''
                    #We check whether the whole bubble is an INV
                    if re.search('<',sv[0]) is not None and len(del_svs) !=0:
                        if sorted(del_svs[0]) == sorted(sv[0].split('<')[1:]) and len(sv) == 1:
                            for saminv in sv[0].split('<')[1:]:
                                inv_size+=len(str(Seq(gfa.loc[saminv,'sequence']).reverse_complement()))
                                inv_alt+=str(Seq(gfa.loc[saminv,'sequence']).reverse_complement())
                            for refinv in del_svs[0]:
                                inv_ref+=gfa.loc[refinv,'sequence']
                            if nodes_sam.index(sv[0]) == 0:
                                ran = 0
                            else:
                                ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                            for x in range(ran):
                                if re.search('<',nodes_ref[x]) is not None:
                                    comps = nodes_ref[x].split('<')
                                    for comp in comps:
                                        if comp =='':
                                            continue
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                else:
                                    vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                            
                            vcf.iloc[n,2] = "INV.%i" %(inv_size)
                            vcf.iloc[n,4] = inv_alt
                            vcf.iloc[n,3] = inv_ref
                            if len(gfa.loc[sv[0].split('<')[1],'sequence']) > 49:
                                print_vcf(vcf,n)
                            inv_bub=True
                            break
                    #we're going to look through each SV's nodes and sum them up to define the SV
                    if inv_bub == True:
                        continue
                    for node in sv:
                        current_p = vcf.iloc[n,1]
                        if re.search('<',node) is not None:
                            #for the first round looking through each node, we will check whether it is an INV (otherwise it will be confused as an INS)
                            for nodes in node.split('<'):
                                if nodes == '' or len(gfa.loc[nodes,'sequence']) < 50:
                                    continue
                                for s in del_svs:
                                    for d in s:
                                        if re.search('<',d) is not None:
                                            for comp in d.split('<'):
                                                #A sample node will only be compared with a reference node if it is larger than 50bp and +-5% of the reference node length
                                                #Possibility I am not taking into account; a group of nodes constitute an INV; I have yet to see this happen, though
                                                if comp == '' or len(gfa.loc[comp,'sequence']) < 50 or len(gfa.loc[comp,'sequence']) < len(gfa.loc[nodes,'sequence']) * 0.95 or len(gfa.loc[comp,'sequence']) > len(gfa.loc[nodes,'sequence']) * 1.05 or nodes == comp:
                                                    continue
                                                if d.split('<').index(comp) == 0 and node.split('<').index(nodes) == 0:
                                                    if d.split('<')[0] != '' and node.split('<')[0] != '':
                                                        score = int(aligner.score(gfa.loc[nodes,'sequence'],str(Seq(gfa.loc[comp,'sequence']).reverse_complement())))
                                                elif d.split('<').index(comp) != 0 and node.split('<').index(nodes) == 0:
                                                    if node.split('<')[0] != '':
                                                        score = int(aligner.score(gfa.loc[nodes,'sequence'],gfa.loc[comp,'sequence']))
                                                elif d.split('<').index(comp) == 0 and node.split('<').index(nodes) != 0:
                                                    if d.split('<')[0] != '':
                                                        score = int(aligner.score(gfa.loc[nodes,'sequence'],gfa.loc[comp,'sequence']))
                                                elif d.split('<').index(comp) != 0 and node.split('<').index(nodes) != 0:
                                                    score = int(aligner.score(gfa.loc[nodes,'sequence'],str(Seq(gfa.loc[comp,'sequence']).reverse_complement())))
                                                #if the sample node reverse complement aligns at a 95% accuracy or more, we have found an INV
                                                #print(nodes,comp)
                                                if score > int(len(gfa.loc[nodes,'sequence']) * 0.95):
                                                    #print(nodes,comp,node)
                                                    #print(nodes_ref,nodes_sam)
                                                    #print(del_svs,ins_svs)
                                                    #if we find an INV, we print it!
                                                    if sv.index(node) == 0:
                                                        if nodes_sam.index(sv[0]) == 0:
                                                            ran = 0
                                                        else:
                                                            if re.search('<',nodes_sam[nodes_sam.index(sv[0])]) is not None and nodes_sam[nodes_sam.index(sv[0])].split('<')[0] in nodes_ref:
                                                                ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])].split('<')[0])+1
                                                            else:
                                                                ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                                        for x in range(ran):
                                                            if re.search('<',nodes_ref[x]) is not None:
                                                                comps = nodes_ref[x].split('<')
                                                                for comp in comps:
                                                                    if comp =='':
                                                                        continue
                                                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                            else:
                                                                vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                                    pos_other = vcf.iloc[n,1]
                                                    for x in range(sv.index(node)):
                                                        if re.search('<',sv[x]) is not None:
                                                            comps = sv[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            vcf.iloc[n,1] += len(gfa.loc[sv[x],'sequence'])
                                                    for x in range(s.index(d)):
                                                        if re.search('<',s[x]) is not None:
                                                            comps = s[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                pos_other += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            pos_other += len(gfa.loc[s[x],'sequence'])
                                                    if vcf.iloc[n,1] + 1000 > pos_other and vcf.iloc[n,1] - 1000 < pos_other:
                                                        vcf.iloc[n,2] = "INV.%i" %(len(gfa.loc[nodes,'sequence']))
                                                        vcf.iloc[n,3] = str(Seq(gfa.loc[nodes,'sequence']).reverse_complement())
                                                        vcf.iloc[n,4] = gfa.loc[nodes,'sequence']
                                                        if len(gfa.loc[nodes,'sequence']) > 49:
                                                            print_vcf(vcf,n)
                                                        sv_inv.append(nodes)
                                                        sv_inv.append(comp)
                                                        break
                                                else:
                                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        #if the node we are comparing is small or not the approximate same length as our candidate inversion, skip it
                                        else:
                                            if len(gfa.loc[d,'sequence']) < 50 or len(gfa.loc[d,'sequence']) < len(gfa.loc[nodes,'sequence']) * 0.95 or len(gfa.loc[d,'sequence']) > len(gfa.loc[nodes,'sequence']) * 1.05 or nodes == d:
                                                continue
                                            if node.split('<').index(nodes) == 0:
                                                if node.split('<')[0] == '':
                                                    score = int(aligner.score(gfa.loc[nodes,'sequence'],str(Seq(gfa.loc[nodes,'sequence']).reverse_complement())))
                                                else:
                                                    if d == node:
                                                        continue
                                                    score = int(aligner.score(gfa.loc[nodes,'sequence'],gfa.loc[d,'sequence']))
                                            else:
                                                score = int(aligner.score(gfa.loc[nodes,'sequence'],gfa.loc[d,'sequence']))
                                            #print(nodes,d,score)
                                            #if the alignment is more than 95% accurate, we have found an inv!
                                            if score > int(len(gfa.loc[nodes,'sequence']) * 0.95):
                                                #if we find an INV, calculate its position on the reference genome, print it!!
                                                if sv.index(node) == 0:
                                                    if nodes_sam.index(sv[0]) == 0:
                                                        ran = 0
                                                    else:
                                                        if re.search('<',nodes_sam[nodes_sam.index(sv[0])]) is not None and nodes_sam[nodes_sam.index(sv[0])].split('<')[0] in nodes_ref:
                                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])].split('<')[0])+1
                                                        else:
                                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                                    for x in range(ran):
                                                        if re.search('<',nodes_ref[x]) is not None:
                                                            comps = nodes_ref[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                                pos_other = vcf.iloc[n,1]
                                                for x in range(sv.index(node)):
                                                        if re.search('<',sv[x]) is not None:
                                                            comps = sv[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            vcf.iloc[n,1] += len(gfa.loc[sv[x],'sequence'])
                                                for x in range(s.index(d)):
                                                        if re.search('<',s[x]) is not None:
                                                            comps = s[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                pos_other += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            pos_other += len(gfa.loc[s[x],'sequence'])
                                                if vcf.iloc[n,1] + 1000 > pos_other and vcf.iloc[n,1] -1000 < pos_other:
                                                    vcf.iloc[n,2] = "INV.%i" %(len(gfa.loc[nodes,'sequence']))
                                                    vcf.iloc[n,3] = str(Seq(gfa.loc[nodes,'sequence']).reverse_complement())
                                                    vcf.iloc[n,4] = gfa.loc[nodes,'sequence']
                                                    if len(gfa.loc[nodes,'sequence']) > 49:
                                                        print_vcf(vcf,n)
                                                    sv_inv.append(nodes)
                                                    sv_inv.append(d)
                                                    break
                                            else:
                                                vcf.iloc[n,1] += len(gfa.loc[nodes,'sequence'])
                            vcf.iloc[n,1] = current_p
                            inv = 0
                            #if the node is an INV, split the current SV and do not take this node into account
                            for comp in node.split('<'):
                                if comp in sv_inv:
                                    if seq == "":
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        continue
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,3] = '*'
                                    vcf.iloc[n,4] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                    #we now change the position to the new SV
                                    vcf.iloc[n,1] += len(seq)
                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    #now we set up the new seq
                                    seq = ""
                                    inv = 1
                                    break
                            if inv == 1:
                                continue
                            extra = 0
                            first = node.split('<')[0]
                            node_badmaybe = []
                            total = 0
                            #Let's look at those nodes that are in fact part of the reference (due to "<")
                            for ref in del_svs:
                                for r in ref:
                                    if r == first:
                                        extra = 1
                                        node_badmaybe.append(r)
                            for comp in node.split('<'):
                                if comp != '':
                                    total += 1
                                for ref in del_svs:
                                        for r in ref:
                                            if re.search('<',r) is not None:
                                                for ir in r.split('<'):
                                                    if ir == comp:
                                                        extra = 1
                                                        node_badmaybe.append(comp)
                            #print(node_badmaybe,extra)
                            #if the whole node is part of the ref
                            if len(node_badmaybe) == total:
                                if sv.index(node) == 0:
                                    if nodes_sam.index(sv[0]) == 0:
                                        ran = 0
                                    else:
                                        ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                    for x in range(ran):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                if seq == "":
                                    for comp in node.split('<'):
                                        if comp == '':
                                            continue
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    continue
                                vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                vcf.iloc[n,3] = '*'
                                vcf.iloc[n,4] = seq
                                if len(seq) > 49:
                                    print_vcf(vcf,n)
                                #we now change the position to the new SV
                                vcf.iloc[n,1] += len(seq)
                                comps = node.split('<')
                                for comp in comps:
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                #now we set up the new seq
                                seq = ""
                                continue
                            #if the node we're looking at is actually part of the reference, we will not take it into account
                            if extra == 1:
                                #if the reference node is at the beginning of the SV, we will calculate the length
                                if sv.index(node) == 0:
                                    if nodes_sam.index(sv[0]) == 0:
                                        ran = 0
                                    else:
                                        if re.search('<',nodes_sam[nodes_sam.index(sv[0])]) is not None and nodes_sam[nodes_sam.index(sv[0])].split('<')[0] in nodes_ref:
                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])].split('<')[0])+1
                                        else:
                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                    for x in range(ran):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp == '':
                                            continue
                                        elif comp in node_badmaybe:
                                            if seq == "":
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                continue
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,3] = '*'
                                            vcf.iloc[n,4] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                            #we now change the position to the new SV
                                            vcf.iloc[n,1] += len(seq)
                                            vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            #now we set up the new seq
                                            seq = ""
                                            comps = node.split('<')
                                        #if the node is the very very last and not part of the reference, print it
                                        elif sv.index(node) == len(sv)-1 and comps.index(comp) == len(comps)-1:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,3] = '*'
                                            vcf.iloc[n,4] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                                break
                                        
                                        else:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                #if it is not at the beginning of the SV, we have to start a new SV, after printing the previous one
                                else:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp == '':
                                            continue
                                        elif comp in node_badmaybe:
                                            if seq == "":
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                continue
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,3] = '*'
                                            vcf.iloc[n,4] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                            #we now change the position to the new SV
                                            vcf.iloc[n,1] += len(seq)
                                            vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            #now we set up the new seq
                                            seq = ""
                                            comps = node.split('<')
                                        #if the node is the very very last and not part of the reference
                                        elif sv.index(node) == len(sv)-1 and comps.index(comp) == len(comps)-1:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,3] = '*'
                                            vcf.iloc[n,4] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                        
                                        else:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                        
                            #if the node is not part of the reference, we will use it for the SV. Now we must see whether the present node is at the beginning or end of the SV. Beginning: we update the genome position and add the sequence, End: we print the SV, Neither: just add the sequence.
                            else:
                                if sv.index(node) == 0:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp =='':
                                            continue
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                    
                                    if nodes_sam.index(sv[0]) == 0:
                                        ran = 0
                                    else:
                                        if re.search('<',nodes_sam[nodes_sam.index(sv[0])]) is not None and nodes_sam[nodes_sam.index(sv[0])].split('<')[0] in nodes_ref:
                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])].split('<')[0])+1
                                        else:
                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                    for x in range(ran):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    
                                    #if this is the last node of the SV, we can print it
                                    if sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,3] = '*'
                                        vcf.iloc[n,4] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                elif sv.index(node) == len(sv)-1:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp =='':
                                            continue
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,3] = '*'
                                    vcf.iloc[n,4] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                else:
                                    comps = node.split('<')
                                    for comp in comps:
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                        #if the current node is not inversed, we must still check whether the node is present in an inversion of the reference list of nodes. For the rest, repeat code.
                        else:
                            #print(node)
                            current_p = vcf.iloc[n,1]
                            #We're going to determine whether one of our nodes is an INV, before anything else
                            for s in del_svs:
                                if len(gfa.loc[node,'sequence']) < 50:
                                    continue    
                                for d in s:
                                    if re.search('<',d) is not None:
                                        for comp in d.split('<'):
                                            if d.split('<').index(comp) == 0:
                                                if comp == '' or len(gfa.loc[comp,'sequence']) < 50 or len(gfa.loc[comp,'sequence']) < len(gfa.loc[node,'sequence']) * 0.95 or len(gfa.loc[comp,'sequence']) > len(gfa.loc[node,'sequence']) * 1.05:
                                                    continue
                                                elif d.split('<')[0] != '':
                                                    score = int(aligner.score(gfa.loc[node,'sequence'],str(Seq(gfa.loc[comp,'sequence']).reverse_complement())))
                                            else:
                                                if comp == node or len(gfa.loc[comp,'sequence']) < 50 or len(gfa.loc[comp,'sequence']) < len(gfa.loc[node,'sequence']) * 0.95 or len(gfa.loc[comp,'sequence']) > len(gfa.loc[node,'sequence']) * 1.05:
                                                    continue
                                                score = int(aligner.score(gfa.loc[node,'sequence'],gfa.loc[comp,'sequence']))
                                            #print(node,comp,score)
                                            if score > int(len(gfa.loc[node,'sequence']) * 0.95):
                                                #if we find an INV, we print it
                                                if sv.index(node) == 0:
                                                    if nodes_sam.index(sv[0]) == 0:
                                                        ran = 0
                                                    else:
                                                        ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                                    for x in range(ran):
                                                        if re.search('<',nodes_ref[x]) is not None:
                                                            comps = nodes_ref[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                                pos_other = vcf.iloc[n,1]
                                                for x in range(sv.index(node)):
                                                        if re.search('<',sv[x]) is not None:
                                                            comps = sv[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            vcf.iloc[n,1] += len(gfa.loc[sv[x],'sequence'])
                                                for x in range(s.index(d)):
                                                        if re.search('<',s[x]) is not None:
                                                            comps = s[x].split('<')
                                                            for comp in comps:
                                                                if comp =='':
                                                                    continue
                                                                pos_other += len(gfa.loc[comp,'sequence'])
                                                        else:
                                                            pos_other += len(gfa.loc[s[x],'sequence'])
                                                if vcf.iloc[n,1] + 1000 > pos_other and vcf.iloc[n,1] - 1000 < pos_other:
                                                    vcf.iloc[n,2] = "INV.%i" %(len(gfa.loc[node,'sequence']))
                                                    vcf.iloc[n,3] = str(Seq(gfa.loc[node,'sequence']).reverse_complement())
                                                    vcf.iloc[n,4] = gfa.loc[node,'sequence']
                                                    if len(gfa.loc[node,'sequence']) > 49:
                                                        print_vcf(vcf,n)
                                                    sv_inv.append(node)
                                                    sv_inv.append(comp)
                                                    break
                                            else:
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    else:
                                        if node == d or len(gfa.loc[d,'sequence']) < 50 or len(gfa.loc[d,'sequence']) < len(gfa.loc[node,'sequence']) * 0.95 or len(gfa.loc[d,'sequence']) > len(gfa.loc[node,'sequence']) * 1.05:
                                            continue
                                        #print(node,d)
                                        if int(aligner.score(gfa.loc[node,'sequence'],str(Seq(gfa.loc[d,'sequence']).reverse_complement()))) > int(len(gfa.loc[node,'sequence']) * 0.95):
                                            #if we find an INV, we print it!
                                            if sv.index(node) == 0:
                                                if nodes_sam.index(sv[0]) == 0:
                                                    ran = 0
                                                else:
                                                    ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                                for x in range(ran):
                                                    if re.search('<',nodes_ref[x]) is not None:
                                                        comps = nodes_ref[x].split('<')
                                                        for comp in comps:
                                                            if comp =='':
                                                                continue
                                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                    else:
                                                        vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                            pos_other = vcf.iloc[n,1]
                                            for x in range(sv.index(node)):
                                                if re.search('<',sv[x]) is not None:
                                                    comps = sv[x].split('<')
                                                    for comp in comps:
                                                        if comp =='':
                                                            continue
                                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                else:
                                                    vcf.iloc[n,1] += len(gfa.loc[sv[x],'sequence'])
                                            for x in range(s.index(d)):
                                                if re.search('<',s[x]) is not None:
                                                    comps = s[x].split('<')
                                                    for comp in comps:
                                                        if comp =='':
                                                            continue
                                                        pos_other += len(gfa.loc[comp,'sequence'])
                                                else:
                                                    pos_other += len(gfa.loc[s[x],'sequence'])
                                            if vcf.iloc[n,1] + 1000 > pos_other and vcf.iloc[n,1] - 1000 < pos_other:
                                                vcf.iloc[n,2] = "INV.%i" %(len(gfa.loc[node,'sequence']))
                                                vcf.iloc[n,3] = gfa.loc[d,'sequence']
                                                vcf.iloc[n,4] = gfa.loc[node,'sequence']
                                                if len(gfa.loc[node,'sequence']) > 49:
                                                    print_vcf(vcf,n)
                                                sv_inv.append(node)
                                                sv_inv.append(d)
                                                break
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                            vcf.iloc[n,1] = current_p
                            inv = 0
                            if node in sv_inv:
                                if seq == "":
                                    vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                    continue
                                vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                vcf.iloc[n,3] = '*'
                                vcf.iloc[n,4] = seq
                                if len(seq) > 49:
                                    print_vcf(vcf,n)
                                #we now change the position to the new SV
                                vcf.iloc[n,1] += len(seq)
                                vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                #now we set up the new seq
                                seq = ""
                                inv = 1
                                break
                            if inv == 1:
                                continue
                            extra = 0
                            for ref in del_svs:
                                for r in ref:    
                                    try:
                                        if node == r.split('<')[0]:
                                            extra = 1
                                            break
                                    except:
                                        continue
                            #if the node we're looking at is part of the reference, we must determine if it is in the middle of an SV or not.
                            if extra == 1:
                                if len(sv) == 1:
                                    continue
                                else:
                                    if sv.index(node) == 0:
                                        if nodes_sam.index(sv[0]) == 0:
                                            ran = 0
                                        else:
                                            ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                        for x in range(ran):
                                            if re.search('<',nodes_ref[x]) is not None:
                                                comps = nodes_ref[x].split('<')
                                                for comp in comps:
                                                    if comp =='':
                                                        continue
                                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            else:
                                                vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                        vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                        continue
                                    elif sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,3] = '*'
                                        vcf.iloc[n,4] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                    else:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,3] = '*'
                                        vcf.iloc[n,4] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                        #we now change the position to the new SV
                                        vcf.iloc[n,1] += len(seq)
                                        vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                        #now we set up the new seq
                                        seq = ""
                            else:
                                if sv.index(node) == 0:
                                    seq += gfa.loc[node,'sequence']
                                    if nodes_sam.index(sv[0]) == 0:
                                        ran = 0
                                    else:
                                        ran = nodes_ref.index(nodes_sam[nodes_sam.index(sv[0])-1])+1
                                    for x in range(ran):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    
                                    #if this is the last node of the SV, we can print it
                                    if sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,3] = '*'
                                        vcf.iloc[n,4] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                elif sv.index(node) == len(sv)-1:
                                    seq += gfa.loc[node,'sequence']
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,3] = '*'
                                    vcf.iloc[n,4] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                else:
                                    seq += gfa.loc[node,'sequence']
                #Now we will find the DEL
                vcf.iloc[n,1] = pos            
                for sv in del_svs:
                    type_sv="DEL"
                    seq=""
                    vcf.iloc[n,1] = pos
                    if inv_bub==True:
                        break
                    #we're going to look through each SV's nodes and sum them up to define the SV
                    for node in sv:
                        if re.search('<',node) is not None:                            
                            inv = 0
                            for comp in node.split('<'):
                                if comp in sv_inv:
                                    if seq == "":
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        continue
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,4] = '*'
                                    vcf.iloc[n,3] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                    #we now change the position to the new SV
                                    vcf.iloc[n,1] += len(seq)
                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    #now we set up the new seq
                                    seq = ""
                                    inv = 1
                                    break
                            if inv == 1:
                                continue
                            extra = 0
                            first = node.split('<')[0]
                            node_badmaybe = []
                            total = 0
                            for ref in ins_svs:
                                for r in ref:
                                    if r == first:
                                        extra = 1
                                        node_badmaybe.append(r)
                            for comp in node.split('<'):
                                if comp != '':
                                    total += 1
                                for ref in ins_svs:
                                        for r in ref:
                                            if re.search('<',r) is not None:
                                                for ir in r.split('<'):
                                                    if ir == comp:
                                                        extra = 1
                                                        node_badmaybe.append(comp)
                            #print(node_badmaybe,extra)
                            #if the whole node is part of the ref
                            if len(node_badmaybe) == total:
                                if sv.index(node) == 0:
                                    for x in range(nodes_ref.index(sv[0])):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                if seq == "":
                                    for comp in node.split('<'):
                                        if comp == '':
                                            continue
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    continue
                                vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                vcf.iloc[n,4] = '*'
                                vcf.iloc[n,3] = seq
                                if len(seq) > 49:
                                    print_vcf(vcf,n)
                                #we now change the position to the new SV
                                vcf.iloc[n,1] += len(seq)
                                comps = node.split('<')
                                for comp in comps:
                                        vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                #now we set up the new seq
                                seq = ""
                                continue
                            #if the node we're looking at is actually part of the reference, we will not take it into account
                            if extra == 1:
                                #if the reference node is at the beginning of the SV
                                if sv.index(node) == 0:
                                    for x in range(nodes_ref.index(sv[0])):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp == '':
                                            continue
                                        elif comp in node_badmaybe:
                                            if seq == "":
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                continue
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,4] = '*'
                                            vcf.iloc[n,3] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                            #we now change the position to the new SV
                                            vcf.iloc[n,1] += len(seq)
                                            vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            #now we set up the new seq
                                            seq = ""
                                            comps = node.split('<')
                                        #if the node is the very very last and not part of the reference
                                        elif sv.index(node) == len(sv)-1 and comps.index(comp) == len(comps)-1:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,4] = '*'
                                            vcf.iloc[n,3] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                                break
                                        
                                        else:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                #if it is not at the beginning of the SV, we have to start a new SV, after printing the previous one
                                else:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp == '':
                                            continue
                                        elif comp in node_badmaybe:
                                            if seq == "":
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                                continue
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,4] = '*'
                                            vcf.iloc[n,3] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                            #we now change the position to the new SV
                                            vcf.iloc[n,1] += len(seq)
                                            vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            #now we set up the new seq
                                            seq = ""
                                            comps = node.split('<')
                                        #if the node is the very very last and not part of the reference
                                        elif sv.index(node) == len(sv)-1 and comps.index(comp) == len(comps)-1:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                            vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                            vcf.iloc[n,4] = '*'
                                            vcf.iloc[n,3] = seq
                                            if len(seq) > 49:
                                                print_vcf(vcf,n)
                                        
                                        else:
                                            seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                        
                            #if the node is not part of the reference, we will use it for the SV. Now we must see whether the present node is at the beginning or end of the SV. Beginning: we update the genome position and add the sequence, End: we print the SV, Neither: just add the sequence.
                            else:
                                if sv.index(node) == 0:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp =='':
                                            continue
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                    
                                    for x in range(nodes_ref.index(sv[0])):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    
                                    #if this is the last node of the SV, we can print it
                                    if sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,4] = '*'
                                        vcf.iloc[n,3] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                elif sv.index(node) == len(sv)-1:
                                    comps = node.split('<')
                                    for comp in comps:
                                        if comp =='':
                                            continue
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,4] = '*'
                                    vcf.iloc[n,3] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                else:
                                    comps = node.split('<')
                                    for comp in comps:
                                        seq += str(Seq(gfa.loc[comp,'sequence']).reverse_complement())
                        #if the current node is not inversed, we must still check whether the node is present in an inversion of the reference list of nodes. For the rest, repeat code.
                        else:
                            inv = 0
                            if node in sv_inv:
                                #print("hi2",node,vcf.iloc[n,1],sv_inv)
                                if seq == "":
                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                    continue
                                vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                vcf.iloc[n,4] = '*'
                                vcf.iloc[n,3] = seq
                                if len(seq) > 49:
                                    print_vcf(vcf,n)
                                #we now change the position to the new SV
                                vcf.iloc[n,1] += len(seq)
                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                #now we set up the new seq
                                seq = ""
                                inv = 1
                                break
                            if inv == 1:
                                continue
                            extra = 0
                            for ref in ins_svs:
                                for r in ref:    
                                    try:
                                        if node == r.split('<')[0]:
                                            extra = 1
                                            break
                                    except:
                                        continue
                            #if the node we're looking at is part of the reference, we must determine if it is in the middle of an SV or not.
                            if extra == 1:
                                if len(sv) == 1:
                                    continue
                                else:
                                    if sv.index(node) == 0:
                                        for x in range(nodes_ref.index(node)):
                                            if re.search('<',nodes_ref[x]) is not None:
                                                comps = nodes_ref[x].split('<')
                                                for comp in comps:
                                                    if comp =='':
                                                        continue
                                                    vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                            else:
                                                vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                        vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                        continue
                                    elif sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,4] = '*'
                                        vcf.iloc[n,3] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                    else:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,4] = '*'
                                        vcf.iloc[n,3] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                        #we now change the position to the new SV
                                        vcf.iloc[n,1] += len(seq)
                                        vcf.iloc[n,1] += len(gfa.loc[node,'sequence'])
                                        #now we set up the new seq
                                        seq = ""
                            else:
                                if sv.index(node) == 0:
                                    seq += gfa.loc[node,'sequence']
                                    
                                    for x in range(nodes_ref.index(node)):
                                        if re.search('<',nodes_ref[x]) is not None:
                                            comps = nodes_ref[x].split('<')
                                            for comp in comps:
                                                if comp =='':
                                                    continue
                                                vcf.iloc[n,1] += len(gfa.loc[comp,'sequence'])
                                        else:
                                            vcf.iloc[n,1] += len(gfa.loc[nodes_ref[x],'sequence'])
                                    #print(vcf.iloc[n,1])
                                    #if this is the last node of the SV, we can print it
                                    if sv.index(node) == len(sv)-1:
                                        vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                        vcf.iloc[n,4] = '*'
                                        vcf.iloc[n,3] = seq
                                        if len(seq) > 49:
                                            print_vcf(vcf,n)
                                elif sv.index(node) == len(sv)-1:
                                    #print("hi")
                                    seq += gfa.loc[node,'sequence']
                                    vcf.iloc[n,2] = "%s.%i" %(type_sv,len(seq))
                                    vcf.iloc[n,4] = '*'
                                    vcf.iloc[n,3] = seq
                                    if len(seq) > 49:
                                        print_vcf(vcf,n)
                                else:
                                    seq += gfa.loc[node,'sequence']
