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


def format_input(filepath: str = None, SNP: str = None, A1: str = None, A2: str = None, BETA: str = None, P: str = None,
                 outdir: str = None):

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
               refpanels: dict,
               image: str,
               bim_file: str,
               summary_stats: List,
               N: List,
               pops: List,
               chrom: int,
               refs_size: int,
               snp_info_file: str,
               out_dir: str,
               meta):

    j = b.new_job(name=f'run-prscsx-{chrom}')

    # get bfile
    input_bfile = b.read_input_group(bim=bim_file)
    snp_info = b.read_input(snp_info_file)

    sst_files_list = []
    sst_files_sizes = 0

    for file in summary_stats:
        input_sst = b.read_input(file)
        sst_files_list.append(input_sst)
        sst_files_sizes += round(10.0 + 2.0 * bytes_to_gb(file))

    # format sst for input
    sst_files = ','.join(sst_files_list)

    # format populations for input
    pop_list = ','.join(pops)

    # format sample sizes for input
    n_list = ','.join(map(str, N))

    job_storage = refs_size + sst_files_sizes + 10

    j.image(image)
    j.cpu(8)
    j.storage(f'{job_storage}Gi')

    j.command('mkdir ref_panels')
    j.command('mkdir tmp_prscsx_output')
    j.command(f'cp {snp_info} ref_panels')

    for ancestry, input_file in refpanels.items():
        j.command(f'tar xf {input_file} -C ref_panels')
    j.command('ls ref_panels')

    j.command(f'''
    python3 PRScsx.py \
        --ref_dir=ref_panels \
        --bim_prefix={input_bfile} \
        --sst_file={sst_files} \
        --n_gwas={n_list} \
        --pop={pop_list} \
        --out_dir=tmp_prscsx_output \
        --out_name=tmp \
        --chrom={chrom} \
        --meta={meta}''')

    j.command(f'mv tmp_prscsx_output {j.scores}')
    b.write_output(j.scores, f'{out_dir}/prs_csx_scores')


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
    prs_img = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscsx'

    refpanels_dict = {}
    anc_list = ['afr', 'amr', 'eas', 'eur', 'sas']
    ref_file_sizes = 0
    for anc in anc_list:
        file = os.path.join(args.ref_path, f'ldblk_ukbb_{anc}.tar.gz')
        refpanels_dict[anc] = b.read_input(file)
        ref_size = bytes_to_gb(file)
        ref_file_sizes += round(10.0 + 2.0 * ref_size)

    for chrom in range(1, 23):
        run_prscsx(b=b, image=prs_img, refpanels=refpanels_dict, bim_file=args.bfile, summary_stats=sst_list,
                   N=n_list, pops=pop_list, chrom=chrom, meta=args.meta, refs_size=ref_file_sizes,
                   snp_info_file=args.snp_info, out_dir=args.out_dir)

    b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--SNP_col', required=True)
    parser.add_argument('--A1_col', required=True)
    parser.add_argument('--A2_col', required=True)
    parser.add_argument('--A1_BETA_col', required=True)
    parser.add_argument('--P_col', required=True)
    parser.add_argument('--bfile', required=True)
    parser.add_argument('--ref_path', required=True)
    parser.add_argument('--snp_info', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--meta', action="store_true")
    arguments = parser.parse_args()

    main(arguments)

    # python prs-csx/PIPELINE_DRAFT.py --input_file tmp_prscsx_inputfile_trialAFR.txt --SNP_col rsid --A1_col ALT --A2_col REF --A1_BETA_col inv_var_meta_beta --P_col inv_var_meta_p --bfile_path 'gs://ukb-diverse-pops/ktsuo_unrelateds_tmp/AFR' --ref_path 'gs://ukb-diverse-pops/prs-csx' --out_dir 'gs://ukb-diverse-pops/prs-csx'
