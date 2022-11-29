import argparse
import hailtop.batch as hb
import pandas as pd
import hail as hl
import os
from typing import Dict, Union
import numpy as np


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


def format_input(filepath: str = None,
                 SNP: str = None,
                 A1: str = None,
                 A2: str = None,
                 BETA: str = None,
                 P: str = None,
                 pheno: str = None,
                 outdir: str = None,
                 snp_info_filepath: str = None):
    """
    Format summary statistics for PRS-CS(x) input
    :param filepath: path to file to be formatted
    :param SNP: path to file to be formatted
    :param A1: path to file to be formatted
    :param A2: path to file to be formatted
    :param BETA: path to file to be formatted
    :param P: path to file to be formatted
    :param pheno: path to file to be formatted
    :param outdir: path to file to be formatted
    :param snp_info_filepath: path to file to be formatted
    :return:
    """

    print(f'1. FORMATTING {filepath}')
    basename = os.path.basename(filepath)
    root, extension = os.path.splitext(basename)

    cols = [SNP, A1, A2, BETA, P]
    summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)
    summstats = summstats.rename(columns={SNP: 'SNP', A1: 'A1', A2: 'A2', BETA: 'BETA', P: 'P'}, inplace=False)
    new_order = ['SNP', 'A1', 'A2', 'BETA', 'P']
    df = summstats[new_order]

    snp_info = pd.read_table(snp_info_filepath, header=0, sep='\t', usecols=['SNP'])
    df_filtered = df.merge(snp_info, how='inner', on='SNP')

    outfile_name = f'{root}'
    # TODO: for prs-csx formatted files are in '{outdir}/formatted_sst_files_for{target_pop}/
    df_filtered.to_csv(f'{outdir}/formatted_sst_files_for{pheno}/{outfile_name}', sep='\t', index=False)


def run_format(args,
               backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               infile_dict: Dict = None,
               method: str = None):

    if method == 'prscs':
        format_b = hb.Batch(backend=backend, name='prsCS-formatting',
                            default_python_image='gcr.io/ukbb-diversepops-neale/prs-csx-python')

        for i in range(0, len(infile_dict['phenos'])):
            print(f'''Analyzing {infile_dict['phenos'][i]}''')
            pheno = infile_dict['phenos'][i]
            sst = infile_dict['ssts'][i]
            SNP_name = infile_dict['SNP_names'][i]
            A1_name = infile_dict['A1_names'][i]
            A2_name = infile_dict['A2_names'][i]
            beta = infile_dict['beta_names'][i]
            pval = infile_dict['pval_names'][i]

            j = format_b.new_python_job(name=f'Formatting: {sst}')
            sst_size = bytes_to_gb(sst)
            disk_size = round(4.0 + 2.0 * sst_size)
            j.storage(disk_size)
            j.cpu(4)
            j.call(format_input, sst, SNP_name, A1_name, A2_name, beta, pval, pheno, args.out_dir, args.snp_info)

        format_b.run()

    if method == 'prscsx':
        format_b = hb.Batch(backend=backend, name='prsCSX-formatting',
                            default_python_image='gcr.io/ukbb-diversepops-neale/prs-csx-python')
        for i in range(0, len(infile_dict['phenos'])):
            print(f'''Analyzing {infile_dict['phenos'][i]}''')
            pheno = infile_dict['phenos'][i]
            sst_list = "".join(infile_dict['ssts'][i]).split(",")
            SNP_names_list = "".join(infile_dict['SNP_names'][i]).split(",")
            A1_names_list = "".join(infile_dict['A1_names'][i]).split(",")
            A2_names_list = "".join(infile_dict['A2_names'][i]).split(",")
            beta_names_list = "".join(infile_dict['beta_names'][i]).split(",")
            pval_names_list = "".join(infile_dict['pval_names'][i]).split(",")

            for sst_i in range(0, len(sst_list)):
                print(sst_list[sst_i], SNP_names_list[sst_i], A1_names_list[sst_i], A2_names_list[sst_i],
                      beta_names_list[sst_i],
                      pval_names_list[sst_i])
                j = format_b.new_python_job(name=f'Formatting: {sst_list[sst_i]}')
                sst_size = bytes_to_gb(sst_list[sst_i])
                disk_size = round(4.0 + 2.0 * sst_size)
                j.storage(disk_size)
                j.cpu(4)
                j.call(format_input, sst_list[sst_i], SNP_names_list[sst_i], A1_names_list[sst_i], A2_names_list[sst_i],
                       beta_names_list[sst_i], pval_names_list[sst_i], pheno, args.out_dir, args.snp_info)

        format_b.run()


def main(args):
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket='ukb-diverse-pops')

    file = pd.read_csv(args.input_file)
    file = file.fillna(np.nan).replace([np.nan], [None])

    steps_list = args.stages.split(',')
    steps = [x.lower() for x in steps_list]

    #############
    # 1. PRS-CS #
    ############
    if args.prs_cs:
        prscs_lists = {'phenos': [], 'ssts': [], 'refpanel_pops': [], 'ns': [], 'SNP_names': [], 'A1_names': [],
                       'A2_names':[], 'beta_names': [], 'pval_names': [], 'target_cohorts': [], 'dup_rsid_files': []}

        for index, row in file.iterrows():
            prscs_lists['phenos'].append(row['phenos'])
            prscs_lists['ssts'].append(row['sum_stats_paths'])
            prscs_lists['refpanel_pops'].append(row['populations'])
            prscs_lists['ns'].append(row['sample_sizes'])
            prscs_lists['SNP_names'].append(row['snp_names'])
            prscs_lists['A1_names'].append(row['A1_names'])
            prscs_lists['A2_names'].append(row['A2_names'])
            prscs_lists['beta_names'].append(row['Beta_names'])
            prscs_lists['pval_names'].append(row['p-value_label_names'])
            prscs_lists['target_cohorts'].append(row['target'])
            prscs_lists['dup_rsid_files'].append(row['dup_ids_file'])

        # 1.1. Formatting
        if 'format' in steps:
            # from format_input import run_format
            run_format(args=args, backend=backend, infile_dict=prscs_lists, method='prscs')
        else:
            print('Summstats already formatted!')

        # 1.2. PRS-CS
        if 'prs' in steps:
            from prscs import run_prscs
            run_prscs(args=args, backend=backend, infile_dict=prscs_lists)

        # 1.3. PLINK
        if 'plink' in steps:
            from plink import run_plink
            run_plink(args=args, backend=backend, infile_dict=prscs_lists, method='prscs')

    ##############
    # 2. PRS-CSX #
    #############
    elif args.prs_csx:
        prscsx_lists = {'phenos': [], 'ssts': [], 'gwas_pops': [], 'ns': [], 'SNP_names': [], 'A1_names': [],
                        'A2_names': [], 'beta_names': [], 'pval_names': [], 'target_pops': [], 'dup_rsid_files': []}

        for index, row in file.iterrows():
            prscsx_lists['phenos'].append(row['phenos'])
            prscsx_lists['ssts'].append([row['sum_stats_paths']])
            prscsx_lists['gwas_pops'].append([row['populations']])
            prscsx_lists['ns'].append([row['sample_sizes']])
            prscsx_lists['SNP_names'].append([row['snp_names']])
            prscsx_lists['A1_names'].append([row['A1_names']])
            prscsx_lists['A2_names'].append([row['A2_names']])
            prscsx_lists['beta_names'].append([row['Beta_names']])
            prscsx_lists['pval_names'].append([row['p-value_label_names']])
            prscsx_lists['target_pops'].append(row['target'])
            prscsx_lists['dup_rsid_files'].append(row['dup_ids_file'])

        # 2.1. Formatting
        if 'format' in steps:
            # from format_input import run_format
            run_format(args=args, backend=backend, infile_dict=prscsx_lists, method='prscsx')
        else:
            print('Summstats already formatted!')

        # 2.2. Run PRS-CSX
        if 'prs' in steps:
            from prscsx import run_prscsx
            run_prscsx(args=args, backend=backend, infile_dict=prscsx_lists)

        # 2.3 PLINK
        if 'plink' in steps:
            from plink import run_plink
            run_plink(args=args, backend=backend, infile_dict=prscsx_lists, method='prscsx')

    else:
        raise SystemExit('No PRS method found. Specify by either including --prs_cs or --prs_csx in your command')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True, help='path to csv-delimited file - SEE EXAMPLE')
    parser.add_argument('--prs_cs', action="store_true", help='include flag if running PRS-CS')
    parser.add_argument('--prs_csx', action="store_true", help='include flag if running PRS-CSx')
    parser.add_argument('--bfile_path', required=True, help='path to bfile of target population, bfiles need to be in format "gs://path/target_cohort_str.bed,bim,fam for PRS-CS or "gs://path/target_pop_str.bed,bim,fam for PRS-CSx; input "gs://path" here')
    parser.add_argument('--ref_path', required=True, help='path to reference panels directory for PRS-CSx')
    parser.add_argument('--snp_info', required=True, help='full path to SNP info file, including filename')
    parser.add_argument('--out_dir', required=True, help='path to output directory')
    parser.add_argument('--meta', action="store_true", help='include flag if you want PRS-CSx meta-analysis option')
    parser.add_argument('--stages', type=str, default='format,prs,plink', help='stage(s) of the pipeline to run')

    arguments = parser.parse_args()

    main(arguments)
