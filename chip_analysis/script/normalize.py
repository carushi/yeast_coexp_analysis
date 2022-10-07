#!/usr/bin/env python
import os
import sys
from optparse import OptionParser
import pandas as pd
# ```
# python normalize.py --process intersect [file_list.csv]
# python normalize.py --process scale [file_list.csv]
# python normalize.py --process compare [file_list.csv]
# python normalize.py --process merge 
#
# For TF dataset
# python normalize.py --tf --process intersect
# ```


DATA_DIR='./'
ANN_DIR='../../data/ann/'
parser = OptionParser()
parser.add_option("-p", "--process", dest="process",
                  help="write report to FILE", metavar="PROCESS")
parser.add_option("-t", "--tf", action='store_true')
parser.add_argument('input', nargs='+', help='input')

(options, args) = parser.parse_args()

job_flag = ['intersect', 'scale', 'compare', 'merge']
file_list = options.input[0]

if options.process not in job_flag:
    print('Option does not match.')
    sys.exit()

if options.tf:
    genotype_list = ['WT', 'ixr1d', 'swi4d', 'swi6d']
else: # CP dataset 
    genotype_list =  ['WT', 'Rad53K', 'Mrc1d']

if options.process != 'merge':
    df = pd.read_csv(file_list)
    for i, row in df.iterrows():
        case = row['case']
        control = row['control']
        rep = row['rep']
        index = row['index']
        TAIL='__trimmed_R1.sorted.bam'
        N_TAIL='__trimmed_R1.sorted_nvalid.bam'
        O_TAIL='__trimmed_R1.sorted_valid.bam' # origin
        valid_file = os.path.join(ANN_DIR, "valid_origin_pm_15000.bed")
        case_file = case+TAIL
        cont_file = cont+TAIL
        FILE = 'chip_case_cont_'+str(index)+'_'+str(rep)+'.bigWig'
        N_FILE = 'chip_case_cont_'+str(index)+'_'+str(rep)+'_non.bigWig'
        O_FILE = 'chip_case_cont_'+str(index)+'_'+str(rep)+'_ori.bigWig'

        if options.process == 'intersect':
            os.system('bedtools intersect -v -a '+case_file+' -b '+valid_file+' > '+case+N_TAIL)
            os.system('bedtools intersect    -a '+case_file+' -b '+valid_file+' > '+case+O_TAIL)
            os.system('bedtools intersect -v -a '+cont_file+' -b '+valid_file+' > '+cont+N_TAIL)
            os.system('bedtools intersect    -a '+cont_file+' -b '+valid_file+' > '+cont+O_TAIL)
        elif options.process == 'scale':
            os.system('samtools index '+cont+N_TAIL)
            os.system('samtools index '+case+N_TAIL)
            os.system('samtools index '+cont+O_TAIL)
            os.system('samtools index '+case+O_TAIL)
            os.system('samtools index '+case+TAIL)
            os.system('samtools index '+cont+TAIL)
            os.system('samtools index '+case+TAIL)
            os.system('samtools index '+cont+TAIL)
            os.system('samtools view -c -F 260 '+case+N_TAIL+' > '+case+'_non_scale.txt')
            os.system('samtools view -c -F 260 '+cont+N_TAIL+' > '+cont+'_non_scale.txt')
        elif options.process == 'compare':
            if not os.path.exists(case+'_non_scale.txt'):
                print('error')
                break
            with open(case+'_non_scale.txt') as f:
                a=f.readline()
            with open(cont+'_non_scale.txt') as f:
                b=f.readline()
            os.system('bamCompare -p 10 --scaleFactorsMethod readCount -bs 25 '+'--bamfile1 '+case+N_TAIL+' --bamfile2 '+cont+N_TAIL+' -o '+N_FILE)
            os.system('bamCompare -p 10 --scaleFactors '+str(int(b))+':'+str(int(a))+' -bs 25 '+'--bamfile1 '+case+O_TAIL+' --bamfile2 '+cont+O_TAIL+' -o '+O_FILE)
            os.system('bigwigCompare -p 10 --operation add --bigwig1 '+N_FILE+' --bigwig2 '+O_FILE+' -o chip_case_cont_'+str(i)+'_'+str(rep)+'_norm.bigWig')
            os.system('bamCompare -p 10 --pseudocount 1 --scaleFactorsMethod readCount -bs 25 '+'--bamfile1 '+case+TAIL+' --bamfile2 '+cont+TAIL+' -o '+FILE)

            os.system('multiBigwigSummary bins -p 10 -bs 25 -b '+' '.join(['chip_case_cont_'+str(i)+'_'+str(rep)+('_'+type if type != '' else '')+'.bigWig' for i in range(len(genotype_list)*len(stage_list)) for rep in range(2)])+
                                        ' --outRawCounts '+'chip_case_cont_merged_'+type+'.tab -o test.npz')
        else:
            pass
else:
    sp = ''
    stage = ['G1', 'HU45', 'HU90']
    if os.path.exists('chip_case_cont_merged_'+sp+'.tab'):
        df = pd.read_csv('chip_case_cont_merged'+sp+'.tab', header=0, sep='\t')
        df.columns = [x.replace('\'', '').replace('#', '') for x in df.columns]
        columns = df.columns[df.columns.str.startswith('chip_case_cont')]
        new_df = df.loc[:,~df.columns.str.startswith('chip_case_cont')]
        # aggregate replicates
        for i in range(9):
            print(df.columns.str.startswith(tuple(['chip_case_cont_'+str(i)+'_'+str(j) for j in range(2)])))
            genotype = genotype_list[int(i/3)]
            stage    = stage_list[int(i%3)]
            new_df.loc[:, genotype+'_'+stage] = df.loc[:, df.columns.str.startswith(tuple(['chip_case_cont_'+str(i)+'_'+str(j) for j in range(2)]))].mean(axis=1).values
        new_df.to_csv('chip_case_cont_mean'+sp+'.tsv', sep='\t', index=False)
        # keep replicates separately
        new_df = df.loc[:,~df.columns.str.startswith('chip_case_cont')]
        for i in range(9):
            for j in range(3):
                genotype = genotype_list[int(i/3)]
                stage    = stage_list[int(i%3)]
                new_df.loc[:, genotype+'_'+stage+'_'+str(j)] = df.loc[:, df.columns.str.startswith('chip_case_cont_'+str(i)+'_'+str(j))].mean(axis=1).values
        new_df.to_csv('chip_case_cont_mean'+sp+'_rep.tsv', sep='\t', index=False)
