import os
import sys

for file in os.listdir('vcftrue'):
    if not file.endswith('.rf'): continue
    prefix = ''
    prefix = '.'.join(os.path.basename(file).split('.')[0:-1])
    vcfFile = 'vcftrue/' + prefix + '.vcf'
    if os.path.exists(vcfFile): continue
    fh = open('vcftrue/' + file)
    oh = open(vcfFile, 'w')
    oh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    for line in fh.readlines():
        pos,ref, alt = line.strip().split()
        oh.write('ref' + '\t' + pos + '\t' + '.' + '\t' + ref + '\t' + alt + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\t' + 'GT' + '\t' + '1\n')
        
    fh.close()
    oh.close()
