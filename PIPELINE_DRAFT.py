import hail as hl
import argparse
import hailtop.batch as hb
import pandas as pd
import os
from typing import List


def bytes_to_gb(in_file: str):
    """
    Convert the size from bytes to GiB
    :param in_file: path to file, str
    :return: file size in GiB
    """

    file_info = hl.utils.hadoop_stat(in_file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def format_input(filepath, SNP, A1, A2, BETA, P, outdir):

    print(f'formatting {filepath}')
    basename = os.path.basename(filepath)
    root, extension = os.path.splitext(basename)

    cols = [SNP, A1, A2, BETA, P]
    summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)
    summstats = summstats.rename(columns={SNP: 'SNP', A1: 'A1', A2: 'A2', BETA: 'BETA', P: 'P'}, inplace=False)
    new_order = ['SNP', 'A1', 'A2', 'BETA', 'P']
    df = summstats[new_order]

    outfile_name = f'{root}_formatted.txt'
    df.to_csv(f'{outdir}/formated_sst_files/{outfile_name}', sep='\t', index=False)


def run_prscsx(b: hb.batch.Batch,
               image: str,
               bim_file: str,
               summary_stats: List,
               N: List,
               pops: List,
               chr: int,
               meta):

    j = b.new_job(name=f'run-prscsx-{chr}')

    j.image(image)
    j.cpu(4)

    # get bfile
    #input_bfile = b.read_input(bim_file)
    bim_basename = os.path.basename(bim_file)
    bim_root, _ = os.path.splitext(bim_basename)

    sst_files = ''

    for file in summary_stats:
        input_sst = b.read_input(file)
        sst_files += f'{input_sst},'

    # format populations for input
    pop_list = ','.join(pops)

    # format sample sizes for input
    n_list = ','.join(map(str, N))

    j.command(f'''
    mkdir tmp_prscsx_output
    python3.8 PRS-CSx.py \
        --ref_dir=ref_panels \
        --bim_prefix={bim_root} \
        --sst_file={sst_files} \
        --n_gwas={n_list} \
        --pop={pop_list} \
        --out_dir=tmp_prscsx_output \
        --out_name=tmp \
        --chrom={chr} \
        --meta={meta}''')


def main(args):
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref',
                                bucket='ukb-diverse-pops')

    sst_list = []
    pop_list = []
    n_list = []
    with open(args.input_file, 'rt') as f:
        rows = f.readlines()

        for row in rows:
            sst_path, pop, n = row.strip().split('\t')
            sst_list.append(sst_path)
            pop_list.append(pop)
            n_list.append(n)

    print(sst_list)
    print(pop_list)
    print(n_list)

    #format_image = hb.docker.build_python_image('gcr.io/ukbb-diversepops-neale/prs-csx-python',
    #                                            requirements=['pandas', 'fsspec', 'gcsfs'])
    format_b = hb.Batch(backend=backend, name='prscsx-formatting',
                        default_python_image='gcr.io/ukbb-diversepops-neale/prs-csx-python')

    for sst in sst_list:
        j = format_b.new_python_job(name=f'formatting-{sst}')
        sst_size = bytes_to_gb(sst)
        disk_size = round(4.0 + 2.0 * sst_size)
        j.storage(disk_size)
        j.cpu(4)
        j.call(format_input, filepath=sst, SNP=args.SNP_col, A1=args.A1_col, A2=args.A2_col, BETA=args.A1_BETA_col,
                     P=args.P_col, outdir=args.out_dir)

    format_b.run()

    # format summstats string names for input
    sst_file_paths = hl.utils.hadoop_ls(f'{args.out_dir}/formated_sst_files')
    sst_list = []
    for i in sst_file_paths:
        sst_list.append(i['path'])

    # get ref panels and untar
    b = hb.Batch(backend=backend, name='prscsx')
    image = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscsx'

    refpanels = {}
    anc_list = ['afr', 'amr', 'eas', 'eur', 'sas']
    for anc in anc_list:
        file = os.path.join(args.ref_path, f'ldblk_ukbb_{anc}.tar.gz')
        refpanels[anc] = b.read_input(file)

    get_refpanels = b.new_job(name='get-ref-panels')
    get_refpanels.command('mkdir ref_panels')

    for ancestry, input_file in refpanels.items():
        get_refpanels.command(f'tar -zwvf {input_file} --directory ref_panels')

    for chrom in range(1, 23):
        run_prscsx(b, image, bim_file=args.bfile_path, summary_stats=sst_list, N=n_list, pops=pop_list, chr=chrom,
                   meta=args.meta)

    b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--SNP_col', required=True)
    parser.add_argument('--A1_col', required=True)
    parser.add_argument('--A2_col', required=True)
    parser.add_argument('--A1_BETA_col', required=True)
    parser.add_argument('--P_col', required=True)
    parser.add_argument('--bfile_path', required=True)
    parser.add_argument('--ref_path', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--meta', action="store_true")
    arguments = parser.parse_args()

    main(arguments)

    # python prs-csx/PIPELINE_DRAFT.py --input_file tmp_prscsx_inputfile_trialAFR.txt --SNP_col rsid --A1_col ALT --A2_col REF --A1_BETA_col inv_var_meta_beta --P_col inv_var_meta_p --bfile_path 'gs://ukb-diverse-pops/ktsuo_unrelateds_tmp/AFR' --ref_path 'gs://ukb-diverse-pops/prs-csx' --out_dir 'gs://ukb-diverse-pops/prs-csx'
