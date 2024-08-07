from miniwalk.main import main
import sys
import pytest

def test_main_command1(mocker):
    mocker.patch('sys.argv', ['miniwalk', 'bench', '-t', 'gns','-c', '~/Documents/PHD/SV/SV_analysis/bench/synthetic_longgraph.vcf','-v','~/Documents/PHD/SV/SV_analysis/bench/synthetic_svimasm.vcf','-r','~/Documents/PHD/SV/SV_analysis/bench/TR_full_100_h37rv.txt','-e','/home/student.unimelb.edu.au/acanaldabalt/Documents/PHD/SV/SV_analysis/bench/TB.REF.GENOME.H37RV.fa'])
    main()