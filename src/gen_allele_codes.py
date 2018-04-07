#!python3

'''
  Plotting clustering results
'''

import gzip as gz
import numpy as np
import pandas as pd
import click


'''
  Generate a set of data types which will 
'''
def data_type_dict(freq_file):
  type_dict = {}
  with gz.open(freq_file, 'rt') as f:
    header = f.readline()
    categories = header.split()
    for i in categories:
      if (i == 'A1') | (i == 'A2'):
        type_dict[i] = 'category'
      elif (i == 'CHR'):
        type_dict[i] = 'category'
      elif (i == 'SNP'):
        type_dict[i] = 'uint32'
      else:
        type_dict[i] = 'uint32'
  return(type_dict)


'''
  Reading in a allele count table and output it to a reasonable format
  Input:
    - freq_file : gzipped file of minor allele counts
    - total_file : gzipped file of total individuals genotyped
    - singletons : filtering or retention of singletons
    - dist : output distance to cluster centroid
  Output:
    - X2 : matrix of raw allele frequencies
    - snplist:  list of snps
    - singleton_snplist : dataframe of singleton snps (and their labels)
'''
def filter_file(freq_file, total_file):
  type_dict_mac = data_type_dict(freq_file) 
  type_dict_total = data_type_dict(total_file) 
  mac_table = pd.read_table(freq_file, dtype=type_dict_mac, low_memory=True)
  total_table = pd.read_table(total_file, dtype=type_dict_total, low_memory=True)
  default = ['CHR', 'SNP', 'A1','A2']
  tru_colnames = mac_table.columns[np.isin(mac_table.columns, default, invert=True)]
  mac = np.array(mac_table[tru_colnames], dtype=np.uint32)
  total_mac = np.sum(mac, axis=1)
  totals = np.array(total_table[tru_colnames], dtype=np.uint32)
  X2 = np.array(mac / totals)
  np.nan_to_num(X2, copy=False)
  snplist = None 
  snplist = mac_table.loc[(total_mac > 0), default]
  X2 = X2[(total_mac > 0),:]
  return(X2, snplist)


'''
  Generate the code for a particular allele
'''
def generate_codes(x,bins):
  code = np.zeros(len(x), dtype=np.uint8)
  zeros = np.where(x == 0)
  code[zeros] = str(0)
  for i in range(1,len(bins)):
    idx = np.where((x <= bins[i]) & (x > bins[i-1]))
    code[idx] = str(i)
  code = np.array(code, dtype=str)
  return(''.join(list(code)))


@click.command()
@click.option('--mac', help='Multidimensional SFS File (MAC)')
@click.option('--total', help='Multidimensional SFS (Total Individauls Genotyped)')
@click.option('--bins', help='Bins for labeling of allele frequencies', default=[0,0.05,1])
@click.option('--out', help='Output Text File')
def main(mac, total, bins, out):
  #Reading in files
  X2, snplist = filter_file(mac, total)
  codes = [generate_codes(X2[i,:], bins) for i in range(X2.shape[0])]
  snplist['Codes'] = codes
  snplist.to_csv(out, sep='\t', compression='gzip', index=False)

if __name__ =='__main__':
    main()

