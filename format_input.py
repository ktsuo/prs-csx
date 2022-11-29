import hailtop.batch as hb
import hail as hl
import os
import pandas as pd
from typing import Dict, Union


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
                       beta_names_list[sst_i], pval_names_list[sst_i], pheno, args.out_dir)

        format_b.run()

