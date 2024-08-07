[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Static Badge](https://img.shields.io/badge/language-Python_3-blue)
<!---[![Build Status]()]()
[![github release version]()]()
[![DOI]()]()
--->

# Miniwalk

Miniwalk is a tool for genotyping SVs from a minigraph graph.

## Motivation

I created this tool as I was interested in capturing SVs from short-read data and there doesn't seem to be such an approach with minigraph graphs as they are incompatible with the vg toolkit.

It is only able to capture SVs from **assemblies**, not raw short-read data. Therefore, it is necessary to create an assembly (i.e. using Flye for long-read data or Shovill for short-read data).

## Installation

### `pip`
```
pip install miniwalk
```

### Install locally
```
git clone git@github.com:aleixcanalda/miniwalk.git

cd miniwalk

pip install .
```

## Pre-requisites

This tool complements the minigraph output, therefore it is necessary to obtain a vcf file using minigraph. Moreover, you will need to have ![k8](https://github.com/attractivechaos/k8) installed in your environment.

This can be achieved using the following commands:

```
#we map the assembly to the graph and get the paths it traverses
minigraph -cxasm --call graph.gfa assembly.fasta > sample.bed

#Run the previous command with the reference genome
minigraph -cxasm --call graph.gfa reference.fasta > reference.bed

#After that, here we are going to paste the reference paths next to our sample's paths, IMPORTANT to keep in this order
paste reference.bed sample.bed | k8 /path/to/minigraph/misc/mgutils.js merge - > merged.bed

#We add the names of the reference and our sample in a file - not too important
echo reference > names.txt
echo sample >> names.txt

#merge2vcf will output a vcf file that points to the paths that are different and we will use that file to determine our SVs
k8 /path/to/minigraph/misc/mgutils-es6.js merge2vcf -s names.txt merged.bed > pan.vcf

```

## Workflow

There are 4 main pipelines within miniwalk: mod, ref, ins2dup and bench.

```
usage: miniwalk [-h] {mod,ref,ins2dup,bench} ...

miniwalk - A tool for genotyping SVs from minigraph graphs.

positional arguments:
  {mod,ref,ins2dup,bench}
                        Choose one of the following pipelines (order for genotyping: mod --> ref --> ins2dup; for pipeline-specific flags run,
                        for example, "miniwalk mod -h")
    mod                 This script takes a bed file outputted from minigraph -xasm --call and the vcf from merge2vcf to modify the vcf file
                        to show the SV type and exact positions. Moreover, the gfa file of the pangenome will also be necessary to determine
                        the exact SV position and length.
    ref                 This script takes the output from mod and refines the vcf file by sorting and clustering SVs
    ins2dup             This script looks at the contiguous bases of each INS to determine if it is instead a DUP. MuMmer must be installed or
                        available locally.
    bench               This benchmarks called SVs from minigraph or manta to a standard

options:
  -h, --help            show this help message and exit
```

There is an inherent order in the way these pipelines work: mod --> ref --> ins2dup

**mod**
mod takes the merge2vcf vcf file and goes through those bubbles where the paths from the reference and the sample are different. Inside each bubble, mod will look at specific nodes that are different, those unique to the reference (DEL/INV) and those unique to the sample (INS/DUP/INV).

```
usage: miniwalk mod [-h] -v VCF -g GFA -o OUT -b BED [-na NA]

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     The vcf file output from merge2vcf.
  -g GFA, --gfa GFA     The gfa pangenome file.
  -o OUT, --output OUT  The new vcf file with the SVs.
  -b BED, --bed BED     The minigraph --call output bed file for finding NAs.
  -na NA                Whether we're interested in highlighting NA regions
```

**ref**
ref takes the mod output and refines it: it sorts the SVs by position and clusters those SVs that are in proximity or completely overlapping each other.

```
usage: miniwalk ref [-h] -v VCF -o OUT

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     The vcf file output from mod.
  -o OUT, --output OUT  The new vcf file with the SVs.
```

**ins2dup**
ins2dup is an extra pipeline I added as I was interested in finding out which insertion calls were actually duplications.

IMPORTANT: this script uses MUMmer's tool "repeat-match" to look for duplications, therefore nucmer must be available in your environment in order to run this tool.

```
usage: miniwalk ins2dup [-h] -c CALL -r REF

options:
  -h, --help            show this help message and exit
  -c CALL, --call CALL  The vcf file output from ref.
  -r REF, --reference REF
                        The reference genome used.
```

**bench**
I created this pipeline as I needed to benchmark the SV calls I got from this tool with calls using other, widely used tools.
It is quite a specific benchmark which may not be useful in many cases. Specifically, I compare my graph vcfs or ![manta](https://github.com/Illumina/manta) vcfs against SVIM-ASM standard vcfs, or a graph vcf.

```
usage: miniwalk bench [-h] -c CALL -v VCF -r REPEAT -e REF -t {gon,gos,gnn,gns,mon,mos,mnn,mns}

options:
  -h, --help            show this help message and exit
  -c CALL, --call CALL  The called vcf file.
  -v VCF, --vcf VCF     The standard vcf file.
  -r REPEAT, --repeat REPEAT
                        The tandem repeat regions found in the genome of MTB. Those regions will be treated differently.
  -e REF, --reference REF
                        The reference genome used.
  -t {gon,gos,gnn,gns,mon,mos,mnn,mns}
                        Choose an option: 
                        - gon: minigraph-called vcf vs minigraph long-read standard; SV-csv output.
                        - gos: minigraph-called vcf vs svim-asm long-read standard; SV-csv output.
                        - gnn: minigraph-called vcf vs minigraph long-read standard; Precision-Recall output.
                        - gns: minigraph-called vcf vs svim-asm long-read standard; Precision-Recall output.
                        - mon: manta-called vcf vs minigraph long-read standard; SV-csv output.
                        - mos: manta-called vcf vs svim-asm long-read standard; SV-csv output.
                        - mnn: manta-called vcf vs minigraph long-read standard; Precision-Recall output.
                        - mns: manta-called vcf vs svim-asm long-read standard; Precision-Recall output."""
```

## Output

The vcf file outputted by my tool uses the INFO column to keep the information of the graph.

The information of the SV is specified in the ID column such as {SV type}.{SV size}. The REF columnholds the nodes found in the reference path whereas the ALT column holds the nodes in the sample paths. In the case of an asterisk '*', there are extra nodes in one of the paths.

Example:
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NC_000962.3 /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/SRR12395048_mod.vcf
NC_000962.3	425339	DEL.60	CGGTGCCGGTATTGCCGATGCCGATGTTTCCGTTTCCGGTGTTGCCGAAGCCGATGTTGC	*	30	PASS	END=425399;AN=2;NS=2;NA=2;ALEN=60,0;AC=1;VS=>s412;VE=>s414;AWALK=>s413,*	GT:GT0	0:0	1:1
NC_000962.3	802477	DUP.162	CGGCCTAGCCCGGCGACGATGCAGAGCGAAGCGATGAGGAGGAGCAGGGC	CGGCCTAGCCCGGCGACGATGCAGAGCGAAGCGATGAGGAGGAGCAGGGCAATGCGGCCTAGCCCGGCGACGATGCAGAGCGAAGCGATGAGGAGGAGCAGGGCACGATGCAGAGCGAAGCGATGAGGAGGAGCAGGGCAATGCGGCCTAGCCCGGCGACGATGCAGAGCGAAGCGATGAGGAGGAGCAGGGCAATGCGGCCTAGCCCGGCG	30	PASS	END=802481;AN=2;NS=2;NA=2;ALEN=54,108;AC=1;VS=>s677;VE=>s681;AWALK=>s678>s679>s680,<s4831>s680	GT:GT0	0:0	1:1
NC_000962.3	917629	DEL.61	GCCCGCCTCCTGCTCATCGCGCTGCGCGCTCTGCATCGTCGCCGGGCTGGGTTGGATTGCC	*	30	PASS	END=917690;AN=2;NS=2;NA=2;ALEN=61,3;AC=1;VS=>s828;VE=>s832;AWALK=>s829>s6647,>s831	GT:GT0	0:0	1:1
NC_000962.3	986670	DEL.146	CGCGTGCGGCTCGCCGGCAACGTTGGCAACATCCCGATTCCCATTGATTGCACGTTGCGCGGCCTAACCCAATATTCCCGGACGAACAACGCCGAGGTCGTGCAGAGCGTCGAGACACACCACCGTCCCGCTAACTTTGATGCCCT	*	30	PASSEND=986816;AN=2;NS=2;NA=2;ALEN=146,0;AC=1;VS=>s951;VE=>s953;AWALK=>s952,*	GT:GT0	0:0	1:1
NC_000962.3	1481197	INS.1674	TCGTCAAGAGCGCCGTGCCAACACCCCAGAACATGAGGTGGCCCACGGCGACTCGCCAGGGGGCCAGGAAGGTGCGGCGCTCCTCCTCACGAGTCGGTTTCCGTCCTTCGATCGCCCAGCGCAGGGCTTGCACGGTCTGCCTGGTCAGTGCGTAGCTACCCAAAGCGAGGGCTAGCAGGACATAGCCCGGTACCACCCCGAACGTGAGCCACCGTGGCGTGTCGCGAACGATGCTCGGTTCGGGGATGGCGATCGTCACCAATAGCAGGGCAACCCCGATGCCGAGCAGGTTCGCGGTCACGACCAGCGCGGTCAGCATGACCTGGATCCGTACCCGTCGGCGCCGTTGGCTTTCCGAAACCCGCCCAAGCAGCCAGGAGCCGTACGCGGGAGTTTCTGGCAGCCGGCCGCTCTGCCGGGTCACCGTCTCCAGCACCCGACCCAAGCGTTGCGCCGTGCTCTTCTTGGCCGACATTGTGGCGTCAGACTAGTTTGTCGAAGAGTCGGGTGCGACCGGTTGGCGCGCTCGTGTTGTTTGCCCGGCTTAGGTGGGCACGGCCAGCCGAGTCGGCTGCTCATGTCCGCGCAGCGTCACCGTCTCGCCCAAAGACCAATGGGCACGTTCGGTTTCGCTGGCAGCGTGCAGTGTGTCCGAGGATGCTAGCAATCGCGCGGGGTGTGATTTGGCCAGTTCGCACAATCGGGCCGCCTGGTTGACCGGCTTGCCGACCACTGTGTATTCGAATCTTTGCTTGGCGCCGACATTGCCGGCGACGATCTGGCCTGCCGCCACCCCGATGCCGGCTTGGACCTCGGGCATCTCGTTGGCCAGCCGATCGGCTATGGCC	TCGTCAAGAGCGCCGTGCCAACACCCCAGAACATGAGGTGGCCCACGGCGACTCGCCAGGGGGCCAGGAAGGTGCGGCGCTCCTCCTCACGAGTCGGTTTCCGTCCTTCGATCGCCCAGCGCAGGGCTTGCACGGTCTGCCTGGTCAGTGCGTAGCTACCCAAAGCGAGGGCTAGCAGGACATAGCCCGGTACCACCCCGAACGTGAGCCACCGTGGCGTGTCGCGAACGATGCTCGGTTCGGGGATGGCGATCGTCACCAATAGCAGGGCAACCCCGATGCCGAGCAGGTTCGCGGTCACGACCAGCGCGGTCAGCATGACCTGGATCCGTACCCGTCGGCGCCGTTGGCTTTCCGAAACCCGCCCAAGCAGCCAGGAGCCGTACGCGGGAGTTTCTGGCAGCCGGCCGCTCTGCCGGGTCACCGTCTCCAGCACCCGACCCAAGCGTTGCGCCGTGCTCTTCTTGGCCGACATTGTGGCGTCAGACTAGTTTGTCGAAGAGTCGGGTGCGACCGGTTGGCGCGCTCGTGTTGTTTGCCCGGCTTAGGTGGGCACGGCCAGCCGAGTCGGCTGCTCATGTCCGCGCAGCGTCACGGTTTCGCCCAAAGACCAATGGGCACGTTCGGTTTCGCTGGCAGCGTGCAGTGTGTCCGAGGATGCTAGCAATCGCGCGGGGTGTGATTTGGCCAGTTCGCACAATCGGGCCGCCTGGTTGACCGGCTTGCCGACCACTGTGTATTCGAATCTTTGCTTGGCGCCGACATTGCCGGCGACGATCTGGCCTGCCGCCACCCCGATGCCGGCTTGGACCTCGGGCATCTCGTTGGCCAGCCGATCGGCTATGGCCCGGGCGGCGGCCAGCGCGGCGTCTTCGGGACGGTCGAGGCGGTTCGGGGCTCCGAAGATGGCCAGGGCGGCGTCGCCTGCGAACTTGTTGATCAGTCCGTGGTGACGGTCGACCTCGTTGACGACGATCGCGAAAAACCGGTTGAGGAGCTTGACCACGTGGGCGGCAGGTTGGTTGTCCACCAGCTGGGTGGAGCCGACGATGTCGACGAAGACGACGGCGGCGTGGCGGTCTTCGCCGCCTAGCTGTGGTCGTTCACGCTCGGCGGCGGCGGCGACTTCGCGTCCGACGTGGCGGCCGAAAAGGTCGCGCACGCGTTCGCGCTCGCGCAGGCCGTTGACCATCGCGTTGAAACCACGCTGCAGCTCACCGAGTTCGGTGCCGTCGAACACCACCAGATCCCCTCGCAGATCCCCCTGCTCGACACGCTTGAGCGCAGCGCGCACCACTCGCACCGGCGCCGCCGTCAGCCACGACATGATCCACATCAAAAGGAACCCGCCGATCAGCGTTGCCGACGACACGATCAGCACACCGTATGCGAACTGGATGTCGGTGAAGTTCTGCAGCACCATCTCAAGGATCCCTAGCAAGGCGATACCGAAGAGCGGCACGCCCGAGCTGAGCAGCCACACCACCATGGTCCGGCCGATGATTCCCAGTGCGAACCGTCCCGGCGGTAGCCCGGCCTCGAGTGCCTGGGCGGCCACCGGGCGCAGAGCGAACTCGGTGAGCATGTAGCAGGCGGTGGCGACCAGAACGCCGCAAAACAGTGTCGAGAACATGACCTGCGGGATGAATAGGCGGTTAGCCAACCCGTAGAGGGTGGTAAACAACGCCGCTCCGACTGCCCACAGGATGAGGGTCCCGACCGCCACCCGCAGCGGGGCCAGCAGGGTGTTGCGTTGGTCGGCCCGACTTGGCGTGTGTTCCTCGCTCGCCCACCGTAGCGCGGCCAAGGTCTGCACCTTGACGTAGTAGCTGCCCAACGCCATCCCGATGACCGTATAAGCGGGTCCGGCCACGAACGTGATCCACAAGGGTGCGTCGCTGACGATGCTAGGTGTCGGAAAAGCGACATTTTCTAACAGCACAGCCAGGCCGATCCCTACCAGGTTTGCTGTGATGAGGGCGATGGTCAGGATGATCTGAATCCGCACCCGGCGGCGGAGCGGGCGTTCTGACACCCACCCAAGCAGCCAGGAGCCGTACGCGGGCGGATCTCGCAGCCGGTCGTGCTGACGGGTCATCGTCTCCAGAACCCGGCGCAGGCGTTGCTCCCGCGTCTTGGACATGGTGGGTCAGCCTAATCTGTCGGTCCTGATTTCCTCGGCATGCTCTGCGGTGAGGTGGATGCTTCCCGGTGGCCGCCGATGCGGCGAGGCTAGTTGGGTGGGCTGGTGGTAACCGCGCAGCGTCACGGTTTCGCCCAAAGACCAATGGGCACGTTCGGTTTCGCTGGCAGCGTGCAGTGTGTCCGAGGATGCTAGCAATCGCGCGGGGTGTGATTTGGCCAGTTCGCACAATCGGGCCGCCTGGTTGACCGGCTTGCCGACCACTGTGTATTCGAATCTTTGCTTGGCGCCGACATTGCCGGCGACGATCTGGCCTGCCGCCACCCCGATGCCGGCTTGGACCTCGGGCATTTCGCTGGCTAGCCGGTCGGCGATGGCC	30	PASS	END=1482758;AN=2;NS=2;NA=2;ALEN=4295,5969;AC=1;VS=>s1399;VE=>s1413;AWALK=>s1400>s1401>s1402>s1403>s1404>s1405>s1406>s1407>s1408>s1409>s1410>s1411>s1412,>s1400>s1401<s8602<s6587<s6586<s6585<s7001<s7000<s6999<s6998<s6997<s6996<s6584>s1406>s1407>s1408>s1409>s1410>s1411>s1412	GT:GT0	0:0	1:1
```

## Limitations

* Only works with assemblies. Short-read assemblies may introduce more errors captured as false positives (even though this still remains the more precise approach compared to mapping raw short reads to a linear reference)
* Tested only on *Mycobacterium tuberculosis*, not tested on eukaryotes (where heterozygosity might be important).
* Inversions are only detected when they form a single node or an entire bubble, but not contiguous nodes in a bigger bubble; certain INVs could be missed.