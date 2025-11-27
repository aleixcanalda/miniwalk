import argparse
from miniwalk.commands import ref,mod,ins2dup,minisr,pr_graph_overlap_output4graph_svimasm,pr_graph_overlap_output4graph,pr_graph_overlap,pr_graph_overlap_svimasm,pr_manta_overlap_output4graph_svimasm,pr_manta_overlap_output4graph,pr_manta_overlap,pr_manta_overlap_svimasm

def main():
    parser = argparse.ArgumentParser(description="miniwalk - A tool for genotyping SVs from minigraph graphs.")
    subparsers = parser.add_subparsers(dest='command', help='Choose one of the following pipelines (order for genotyping: mod --> ref --> ins2dup; for pipeline-specific flags run, for example, "miniwalk mod -h")')

    # mod
    parser_command1 = subparsers.add_parser('mod', help='This script takes a bed file outputted from minigraph -xasm --call and the vcf from merge2vcf to modify the vcf file to show the SV type and exact positions. Moreover, the gfa file of the pangenome will also be necessary to determine the exact SV position and length.')
    
    parser_command1.add_argument('-v', '--vcf',
                        dest = "vcf",
                        action = "store",
                        required = True,
                        help = "The vcf file output from merge2vcf.")

    parser_command1.add_argument('-g', '--gfa',
                        dest = "gfa",
                        action = "store",
                        required = True,
                        help = "The gfa pangenome file.")

    parser_command1.add_argument('-o', '--output',
                        dest = "out",
                        action = "store",
                        required = True,
                        help = "The new vcf file with the SVs.")

    parser_command1.add_argument('-b', '--bed',
                        dest = "bed",
                        action = "store",
                        required = True,
                        help = "The minigraph --call output bed file for finding NAs.")

    parser_command1.add_argument('-na',
                        dest = "na",
                        action = "store_true",
                        required = False,
                        help = "Whether we're interested in highlighting NA regions")

    # ref
    parser_command2 = subparsers.add_parser('ref', help='This script takes the output from mod and refines the vcf file by sorting and clustering SVs')
    parser_command2.add_argument('-v', '--vcf',
                    dest = "vcf",
                    action = "store",
                    required = True,
                    help = "The vcf file output from mod.")

    parser_command2.add_argument('-o', '--output',
                        dest = "out",
                        action = "store",
                        required = True,
                        help = "The new vcf file with the SVs.")

    #ins2dup
    parser_command3 = subparsers.add_parser('ins2dup', help='This script looks at the contiguous bases of each INS to determine if it is instead a DUP. MuMmer must be installed or available locally.')    
    parser_command3.add_argument('-c', '--call',
                    dest = "call",
                    action = "store",
                    required = True,
                    help = "The vcf file output from ref.")
    parser_command3.add_argument('-r', '--reference',
                        dest = "ref",
                        action = "store",
                        required = True,
                        help = "The reference genome used.")

    #bench
    parser_command4 = subparsers.add_parser('bench', help='This benchmarks called SVs from minigraph or manta to a standard')    
    parser_command4.add_argument('-c', '--call',
                    dest = "call",
                    action = "store",
                    required = True,
                    help = "The called vcf file.")

    parser_command4.add_argument('-v', '--vcf',
                        dest = "vcf",
                        action = "store",
                        required = True,
                        help = "The standard vcf file.")

    parser_command4.add_argument('-r', '--repeat',
                        dest = "repeat",
                        action = "store",
                        required = True,
                        help = "The tandem repeat regions found in the genome of MTB. Those regions will be treated differently.")

    parser_command4.add_argument('-e', '--reference',
                        dest = "ref",
                        action = "store",
                        required = True,
                        help = "The reference genome used.")
    
    parser_command4.add_argument('-t',
                        dest = "type",
                        action = "store",
                        required = True,
                        choices = ["gon","gos","gnn","gns","mon","mos","mnn","mns"],
                        help = """Choose an option: \n
        - gon: minigraph-called vcf vs minigraph long-read standard; SV-csv output.\n
        - gos: minigraph-called vcf vs svim-asm long-read standard; SV-csv output.\n
        - gnn: minigraph-called vcf vs minigraph long-read standard; Precision-Recall output.\n
        - gns: minigraph-called vcf vs svim-asm long-read standard; Precision-Recall output.\n
        - mon: manta-called vcf vs minigraph long-read standard; SV-csv output.\n
        - mos: manta-called vcf vs svim-asm long-read standard; SV-csv output.\n
        - mnn: manta-called vcf vs minigraph long-read standard; Precision-Recall output.\n
        - mns: manta-called vcf vs svim-asm long-read standard; Precision-Recall output.""")
    
    #minisr
    parser_command5 = subparsers.add_parser('minisr', help='This script takes a GAF alignment file and looks at node mapping depth to determine the paths traversed through the graph, creating a BED file ready to be input to miniwalk mod.')
    parser_command5.add_argument('gaf_file', help='Input GAF file')
    parser_command5.add_argument('gfa_file', help='Input GFA file')
    parser_command5.add_argument('bubble_file', help='Input bubble file')
    parser_command5.add_argument('--sample', '-s',dest="sample", required=True, help='Sample name')
    parser_command5.add_argument('--ploidy', '-p',dest="ploidy", choices=['haploid', 'diploid'], 
                        default='diploid', help='Sample ploidy (default: diploid)')
    parser_command5.add_argument('--output', '-o',dest="output", default='output', 
                        help='Output prefix (default: output)')
    parser_command5.add_argument('--min-reads', '-m',dest="min_reads", type=int, default=5,
                        help='Minimum read count for a node to be considered (default: 5)')
    parser_command5.add_argument('--min-uniformity', '-u',dest="min_uniformity", type=int, default=0.3,
                        help='Minimum uniformity score using exponential decay of coefficient of variation across a node (default: 0.3)')
    parser_command5.add_argument('--min-coverage-fraction', '-c',dest="min_coverage_fraction", type=int, default=0.5,
                        help='Minimum coverage across a node to be considered as fully mapped (default: 0.5)')
    parser_command5.add_argument('--min-read-cov', '-r',dest="min_read_cov", type=int, default=0.9,
                        help='Minimum coverage of alignment on read (default: 0.9)')
    parser_command5.add_argument('--min-node-cov', '-n',dest="min_node_cov", type=int, default=0.9,
                        help='Minimum coverage of alignment on node (default: 0.9)')
    parser_command5.add_argument('--max-node-read-ratio', '-rr',dest="max_read_node_ratio", type=int, default=2,
                        help='Maximum number of times to consider a read mapped to a smaller node (default: reads 2 times larger than a node, solely mapped to that node)')

    args = parser.parse_args()

    if args.command == 'mod':
        mod.main(args)
    elif args.command == 'ref':
        ref.main(args)
    elif args.command == 'ins2dup':
        ins2dup.main(args)
    elif args.command == 'minisr':
        minisr.main(args)
    elif args.command == "bench":
        if args.type == "gon":
            pr_graph_overlap_output4graph.main(args)
        elif args.type == "gos":
            pr_graph_overlap_output4graph_svimasm.main(args)
        elif args.type == "gnn":
            pr_graph_overlap.main(args)
        elif args.type == "gns":
            pr_graph_overlap_svimasm.main(args)
        elif args.type == "mon":
            pr_manta_overlap_output4graph.main(args)
        elif args.type == "mos":
            pr_manta_overlap_output4graph_svimasm.main(args)
        elif args.type == "mnn":
            pr_manta_overlap.main(args)
        elif args.type == "mns":
            pr_manta_overlap_svimasm.main(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
