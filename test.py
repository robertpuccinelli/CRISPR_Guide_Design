from crRNA_generator import GuideRNAGeneratorBase as gb
gen = gb('GACCACCCCAAAAAUGAAGGGGACUAAAAC')
gen.seq_target = 'TTACAAACATTGGCCGCAAA'
gen.crRNAGenerate()
print(gen.seq_crRNA)
print(gen.seq_crRNA_DNA_template)