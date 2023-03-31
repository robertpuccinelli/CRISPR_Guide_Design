import logging
from crRNA_generator import GuideRNAGeneratorBase as gb
from crRNA_generator import T7crRNAGenerator as tg

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s [%(name)s]', level=logging.INFO, datefmt='%H:%M:%S')

print('Basic gen')
gen = gb('GACCACCCCAAAAAUGAAGGGGACUAAAAC')
gen.seq_target = 'TTACAAACATTGGCCGCAAA'
gen.crRNAGenerate()

seq_list = ['GGACCCCAAAATCAGCGAAA',\
            'CGCATTACGTTTGGTGGACC',\
            'TTCAACTGGCAGTAACCAGA',\
            'TTACAAACATTGGCCGCAAA',\
            'CAATTTGCCCCCAGCGCTTC',\
            'TTCTTCGGAATGTCGCGCAT',\
            'GGGAGCCTTGAATACACCAA',\
            'CACATTGGCACCCGCAATCC',\
            'AATGCTGCAATCGTGCTACA',\
            'GGGGAACTTCTCCTGCTAGA',\
            'CTTTGCTGCTGCTTGACAGA',\
            'CAGCTTGAGAGCAAAATGTC',\
            'AGGUUCUUGACUACCGUAAU'] # <- last entry is nontargetting

print('T7 gen')
t7gen = tg('GACCACCCCAAAAAUGAAGGGGACUAAAAC')

for target in seq_list:
    t7gen.seq_target = target
    t7gen.crRNAGenerate()
 #   print(t7gen.seq_crRNA_DNA_template)
