import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', type=argparse.FileType('r'))
args = parser.parse_args()

df = pd.read_csv(args.filename.name, sep='\t', header=None)
df_list = []

min_allele_depth = 10
for index, line in df.iterrows():
    
	info_line = line[10]
	allele_depth = info_line.split(':')[1]
	ref = int(allele_depth.split(',')[0])
	alt = int(allele_depth.split(',')[1])
	sum_allele_depth = ref + alt
	if sum_allele_depth > 9:
		corr_ref = (ref / sum_allele_depth) *min_allele_depth 
		corr_alt = (alt / sum_allele_depth) *min_allele_depth
		if corr_ref == 0:
			continue
		elif corr_alt > 2.5:
			df_list.append(line)

pd.DataFrame(df_list).to_csv(args.filename.name + 'second_percent.vcf', sep='\t', index=False, header=False, mode='a')


