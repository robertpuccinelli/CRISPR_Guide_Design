from crRNA_generator import GuideRNAGeneratorBase as gb
from crRNA_generator import T7crRNAGenerator as tg
gen = gb('GACCACCCCAAAAAUGAAGGGGACUAAAAC')
gen.seq_target = 'TTACAAACATTGGCCGCAAA'
gen.crRNAGenerate()
print(gen.seq_crRNA)
print(gen.seq_crRNA_DNA_template)

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
    t7gen.crRNAGenerateT7()
    print(t7gen.seq_crRNA_DNA_template)
