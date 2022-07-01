import re
from multiprocessing import Pool
from typing import Tuple

import pyranges as pr
from tabulate import tabulate
import pandas as pd
import tqdm


coriell_file = 'coriell_chr1.vcf'
refFlat_file = 'refFlat.txt'

refFlat_cols = [
    'geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
]


def get_vcf_names(vcf_path: str):
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line[1:].split('\t')]
                break
    return vcf_names


vcf_ranges: pr.PyRanges


def calc_overlapped(gene_and_group: Tuple[str, pd.DataFrame]) -> Tuple[str, int]:
    gene, group = gene_and_group
    gene_ranges = pr.PyRanges(group[['Chromosome', 'Start', 'End']])
    overlaps = len(vcf_ranges.overlap(gene_ranges))
    return gene, overlaps


def main():
    vcf_col_names = get_vcf_names(coriell_file)
    vcf_df = pd.read_csv(
        coriell_file, comment='#', delim_whitespace=True,
        header=None, names=vcf_col_names,
    )
    compiled_regex = re.compile('DP=(.*?);')
    vcf_df['DP'] = vcf_df['INFO'].apply(lambda x: int(re.search(compiled_regex, x).group(1)))
    vcf_df['POS_END'] = vcf_df['POS'] + vcf_df['DP']
    vcf_df.rename(columns={'CHROM': 'Chromosome', 'POS': 'Start', 'POS_END': 'End'}, inplace=True)
    global vcf_ranges
    vcf_ranges = pr.PyRanges(vcf_df[['Chromosome', 'Start', 'End']])

    refFlat_df = pd.read_csv(
        refFlat_file, sep='\t', header=None,
        names=refFlat_cols, usecols=['geneName', 'chrom', 'txStart', 'txEnd'],
    )
    refFlat_df.rename(columns={'chrom': 'Chromosome', 'txStart': 'Start', 'txEnd': 'End'}, inplace=True)
    refFlat_df_grouped = list(refFlat_df.groupby(['geneName']))

    with Pool() as p:
        results = list(tqdm.tqdm(p.imap(calc_overlapped, refFlat_df_grouped), total=len(refFlat_df_grouped)))

    print(tabulate(results))


if __name__ == '__main__':
    main()
