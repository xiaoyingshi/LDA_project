import sys
import pandas as pd
import argparse
import re
import subprocess
import os

def formatted(prefixes, contact_path, output):
	os.system("mkdir %s/formatted_file"%(output))
	for matrix in prefixes:
		os.system('''awk '{$2=$2" 0 0 "$1;$1="0 "$1;$3=$3" 1";print $0}' %s > %s/formatted_file/%s_formatted'''%(matrix, output, matrix))

def juicer(resolution, juicer_path, max_chr, output):
	os.system("mkdir %s/hic_file"%(output))
	os.system("mkdir %s/eigenvector"%(output))
	formatted_file = os.listdir("%s/formatted_file/"%(output))
	for matrix in formatted_file:
		# os.system("java -Xmx4g -jar %s pre %s/%s %s/hic_file/%s.hic %s"%(juicer_path, output, matrix, output, matrix, species))
		for i in list(range(1, max_chr+1))+['X','Y']:
			os.system("java -jar %s eigenvector NONE %s/hic_file/%s.hic %s -p BP %s %s/eigenvector/%s_eigen_chr%s.txt"%(juicer_path, output, matrix, i, resolution, output, matrix, i))

def coordinate(resolution, prefixes, output):
	eigenvector_file = os.listdir("%s/eigenvector/"%(output))
	os.system("mkdir %s/coordinate"%(output))

	for ei_file in eigenvector_file:
		chr_num = re.findall(".*_eigen_chr(.*).txt",ei_file)[0] #不一定是num
		os.system('''awk -v a=chr%s '{print a,((NR-1)*%s),(NR*%s),$0}' %s/eigenvector/%s > %s/coordinate/%s_tmp'''%(chr_num, resolution, resolution, output, ei_file, output, ei_file))
		os.system('''grep -v 'NaN'  %s/coordinate/%s_tmp >  %s/coordinate/%s_tmp_nona && mv %s/coordinate/%s_tmp_nona %s/coordinate/%s_tmp'''%(output, ei_file,output, ei_file, output, ei_file, output, ei_file))
		os.system("sed 's/ /\t/g' %s/coordinate/%s_tmp | sort -V -k1 > %s/coordinate/%s_tmp1 && mv %s/coordinate/%s_tmp1 %s/coordinate/%s_tmp"%(output, ei_file, output, ei_file, output, ei_file, output, ei_file))

	prefixes = [re.findall("(.*)_eigen.*",i)[0] for i in prefixes]

	# Merge files together
	for i in prefixes:
		os.system("cat %s/coordinate/%s* > %s/coordinate/%s_all_bed "%(output, i, output, i))
	os.system("rm %s/coordinate/*nona"%(output))
	
def check_direction(species, output):
	os.system("mkdir %s/cold_TAD"%(output))
	os.system("mkdir %s/cold_TAD_avg"%(output))

	coordinate_file = os.listdir("%s/coordinate/"%(output))
	coordinate_file_nobed = [i for i in coordinate_file if '_all_bed' not in i]

	for i in list(range(1, max_chr+1))+['X','Y']:
		os.system('''cat /mnt/Storage2/home/zhengrongbin/project/TAD_promoter/data/TAD_outlier_%s/%s_cold_TAD.xls | grep "chr%s$(printf '\t')" > %s/cold_TAD/TAD_cold_chr%s.txt'''%(species, species, i, output, i))
		
		for nobed in coordinate_file_nobed:
			chr_num = re.findall(".*_eigen_chr(.*).txt_tmp$",nobed)[0]
			if chr_num == str(i):
				os.system("bedtools intersect -wao  -b %s/coordinate/%s -a %s/cold_TAD/TAD_cold_chr%s.txt | grep -v 'interTAD'> %s/cold_TAD/%s_chr%s_cold.txt"%(output, nobed, output, i, output, nobed, i))
				os.system('''awk 'NR>1{arr9[$4]  += $9;arr2[$4]  += $2;arr3[$4]  += $3;arr5[$4]  += $5; arr7[$4]  += $7;arr8[$4]  += $8;count[$4] += 1}END{ for (a in arr9) { print a "\t" $1  "\t" arr2[a] / count[a] "\t" arr3[a] / count[a] "\t" arr5[a] / count[a] "\t" arr7[a] / count[a] "\t" arr8[a] / count[a]  "\t" arr9[a] / count[a]}}' %s/cold_TAD/%s_chr%s_cold.txt > %s/cold_TAD_avg/%s_cold_avg'''%(output, nobed, i, output, nobed))

def change_direction(output):
	os.system("mkdir %s/final_bin"%(output))
	change_list = {}  
	all_file = os.listdir("%s/cold_TAD_avg/"%(output))
	pre_list = list(set([re.findall("(.*)_eigen.*",i)[0] for i in all_file]))
	for pre in pre_list:
		change_list[pre] = []

	for file in all_file:
		prefix = re.findall("(.*)_eigen.*",file)[0]
		if os.stat("%s/cold_TAD_avg/"%(output)+file).st_size != 0:
			file1 = pd.read_csv("%s/cold_TAD_avg/"%(output)+file, sep='\t', header=None)
			l1 = file1.iloc[:,7].tolist()
			l2 = [i for i in l1 if i >0]
			l_ratio = len(l2)/len(l1)
			if l_ratio>0.5:
				chr_num = re.findall(".*chr(.*).txt.*",file)[0]
				# print(chr_num)
				change_list[prefix].append(chr_num)
	           

	beds1 = os.listdir("%s/coordinate/"%(output))
	beds = [i for i in beds1 if '_all_bed' in i]

	for bed in beds:
		if bed != '.DS_Store':
			gsm_bed = re.findall("(.*)_all_bed",bed)[0]
			print(bed)
			bed1 = pd.read_csv("%s/coordinate/"%(output)+bed, sep='\t',names=['chr','start','end','PC1'])
			for index,row in bed1.iterrows():
				for i in change_list[gsm_bed]:
					if row['chr'] == 'chr'+ i:
						bed1.loc[index,'PC1'] = - row['PC1']
			# print(bed1.head)
			bed1.to_csv("%s/final_bin/%s_new"%(output, bed), sep='\t', header=False, index = False)
                    
    
def final(TAD_domain, output):
	os.system("mkdir %s/result"%(output))

	new_beds = os.listdir("%s/final_bin/"%(output))
	for newbed in new_beds:
		os.system("bedtools intersect -a %s/final_bin/%s -b %s -wao -f 1.00 |sort|uniq> %s/result/%s_compart.bin.TAD"%(output, newbed, TAD_domain, output, newbed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perdict A/B compartment.')
    parser.add_argument("-c", "--contact_path",dest="contact_path", 
                    help="Input the directory that contains all contact matrix file if you want to start analysis from contact matrix level", metavar="FILE",
                    type=str)
    parser.add_argument("-e", "--eigenvector", dest="eigenvector",
    	help='Input the directory that contains all eigenvector value juicer get if you want to start analysis from eigenvector level',type=str)
    parser.add_argument("-s", "--species",dest="species", required=True,
                    help="Input species (eg:mm10)", metavar="FILE",
                    type=str)    
    parser.add_argument("-r", "--resolution", dest="resolution", required=True,
    	help='Input the resolution for perdict compartment(eg:50000)',
    	type=int)
    parser.add_argument("-t", "--TAD_domain", dest="TAD_domain", nargs='?',
    	default = "/mnt/Storage/home/shixiaoying/Projects/AB_predict/AB_predict/mm10_TAD_domain_withGap_new.xls",
    	help='Input the TAD_domain path',
    	type=str)
    parser.add_argument("-j", "--juicer_path", nargs='?',
    	default="/mnt/Storage/home/shixiaoying/software/juicer/juicer_tools_1.14.08.jar", 
    	help='Input the juicer tool path',type=str)
    parser.add_argument("-o", "--output", dest="output", required=True,
    	help='Input the result directory path',type=str)

    args = parser.parse_args()
    contact_path = args.contact_path
    eigenvector = args.eigenvector
    species = args.species
    resolution = args.resolution
    TAD_domain = args.TAD_domain
    juicer_path = args.juicer_path
    output = args.output

    if re.findall(r'[A-Za-z]+',species)[0] == "mm":
    	max_chr = 19
    else:
    	max_chr = 22

    if eigenvector:
    	os.system("mkdir %s/eigenvector"%(output))
    	os.system("cp %s/* %s/eigenvector"%(eigenvector, output))
    	prefixes = os.listdir(eigenvector)

    	coordinate(resolution, prefixes, output)
    	check_direction(species, output)
    	change_direction(output)
    	final(TAD_domain, output)

    if contact_path:
    	os.system("gunzip %s/*"%(contact_path))
    	prefixes = os.listdir(contact_path)

    	formatted(prefixes, contact_path, output)
    	juicer(resolution, juicer_path, max_chr, output)
    	coordinate(resolution, prefixes, output)
    	check_direction(species, output)
    	change_direction(output)
    	final(TAD_domain, output)


