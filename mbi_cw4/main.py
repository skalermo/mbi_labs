import pandas as pd
from tqdm import tqdm
from tabulate import tabulate


def is_overlapping(x1: int, x2: int, y1: int, y2: int) -> bool:
    return max(x1, y1) <= min(x2, y2)


def is_overlapping_by_80(x1: int, x2: int, y1: int, y2: int) -> bool:
    overlap = min(x2, y2) - max(x1, y1)
    return overlap >= 0.8 * (y2 - y1)


def main():
    finalcall_path = './finalcall.csv'
    cnv = pd.read_csv(finalcall_path, sep=',')
    cnv = cnv[['sample_name', 'chr', 'cnv', 'st_bp', 'ed_bp']]

    dgv_path = 'GRCh37_hg19_variants_2020-02-25.txt'
    dgv = pd.read_csv(dgv_path, sep='\t', header=0,
                      dtype={'chr': str, 'start': int, 'end': int, 'variantsubtype': str, 'samples': str})
    dgv = dgv[['chr', 'start', 'end', 'variantsubtype', 'samples']]
    dgv = dgv[dgv['chr'] == '20']

    del_overlaps_count = 0
    dup_overlaps_count = 0
    overlaps_by_80_count = 0

    for _, row in tqdm(list(cnv.iterrows())):
        x1, x2 = row['st_bp'], row['ed_bp']
        for _, row2 in dgv.iterrows():
            y1, y2 = row2['start'], row2['end']
            match row2['variantsubtype']:
                case 'deletion':
                    del_overlaps_count += is_overlapping(x1, x2, y1, y2)
                case 'duplication':
                    dup_overlaps_count += is_overlapping(x1, x2, y1, y2)
            overlaps_by_80_count += is_overlapping_by_80(x1, x2, y1, y2)

    print(tabulate({
        'Overlapping DGV\'s deletions': del_overlaps_count,
        'Overlapping DGV\'s duplications': dup_overlaps_count,
        'Overlapping by 80% or more': overlaps_by_80_count,
    }.items()))


if __name__ == '__main__':
    main()
