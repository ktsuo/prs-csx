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


def format_input(filepath: str = None, SNP: str = None, A1: str = None, A2: str = None, BETA: str = None, P: str = None, target_pop: str = None, outdir: str = None):
    """
    Format summary statistics for PRS-CSx input
    :return: formatted summary statistics for PRS-CSx
    """

    print(f'formatting {filepath}')
    basename = os.path.basename(filepath)
    root, extension = os.path.splitext(basename)

    cols = [SNP, A1, A2, BETA, P]
    summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)
    summstats = summstats.rename(columns={SNP: 'SNP', A1: 'A1', A2: 'A2', BETA: 'BETA', P: 'P'}, inplace=False)
    new_order = ['SNP', 'A1', 'A2', 'BETA', 'P']
    df = summstats[new_order]

    #outfile_name = f'{root}_formatted.txt'
    outfile_name = f'{root}'
    df.to_csv(f'{outdir}/formatted_sst_files_for{target_pop}/{outfile_name}', sep='\t', index=False)


def run_prscsx(b: hb.batch.Batch,
               depends_on_j,
               image: str,
               bfile: str,
               target_pop: str,
               summary_stats: List,
               N: List,
               pops: List,
               chrom: int,
               ref_panels_dir: str,
               refs_size: int,
               out_dir: str,
               meta):
    """
    Run PRS-CSx using formatted summary statistics
    :return: PRS-CSx output files
    """

    j = b.new_job(name=f'run-prscsx-{chrom}')
    j.depends_on(depends_on_j)

    # get bfile
    input_bim_file = b.read_input_group(bim=f'{bfile}/{target_pop}/{target_pop}.bim')
    #snp_info = b.read_input(snp_info_file)

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

    job_storage = refs_size + sst_files_sizes + 20

    j.image(image)
    j.cpu(8)
    j.storage(f'{job_storage}Gi')
    j.memory('highmem')
    j.command('mkdir tmp_prscsx_output')

    j.command(f'''
    python3 PRScsx.py \
        --ref_dir={ref_panels_dir} \
        --bim_prefix={input_bim_file} \
        --sst_file={sst_files} \
        --n_gwas={n_list} \
        --pop={pop_list} \
        --out_dir=tmp_prscsx_output \
        --out_name=tmp \
        --chrom={chrom} \
        --meta={meta}''')

    j.command(f'mv tmp_prscsx_output {j.scores}')
    b.write_output(j.scores, f'{out_dir}/prs_csx_output_for{target_pop}')

    return j


def run_plink(b: hb.batch.Batch,
            #   depends_on_j_list: list,
              bfile_size: int,
              image: str,
              source_pop: str,
              bfile: hb.ResourceGroup,
              target_pop: str,
              dup_ids_file: str,
              out_dir: str):
    """
    Using PRS-CSx output compute PRS for target population from each input summary statistic in PLINK
    :return: PLINK score files
    """

    j = b.new_job(name=f'run-plink-{source_pop}')
    # j.depends_on(*depends_on_j_list)
    j.image(image)
    j.cpu(8)

    j.memory('highmem')

    job_storage = bfile_size + 20
    j.storage(f'{job_storage}Gi')

    j.command('mkdir tmp_input')
    input_files = hl.hadoop_ls(f'{out_dir}/prs_csx_output_for{target_pop}/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_chr*')
    for file in input_files:
        file = b.read_input(file['path'])
        j.command(f'mv {file} tmp_input')
    j.command('ls tmp_input')

    j.command(f'cat tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_chr* > tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_ALLchr.txt')

    j.command('mkdir tmp_output')
    # exclude duplicate SNPs
    if dup_ids_file:
        dup_ids_input = b.read_input(dup_ids_file)
        j.command(f'plink \
            --bfile {bfile} \
                --score tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum center \
                    --exclude {dup_ids_input} \
                        --out tmp_output/from_{source_pop}')
        j.command(f'mv tmp_output {j.output}')
        b.write_output(j.output, f'{out_dir}/{target_pop}_plink_scores')

    if dup_ids_file is None:
        j.command(f'plink \
            --bfile {bfile} \
                --score tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum center \
                    --out tmp_output/from_{source_pop}')
        j.command(f'mv tmp_output {j.output}')
        b.write_output(j.output, f'{out_dir}/{target_pop}_plink_scores')


def main(args):
    global prscsx_j
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
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

    # print(sst_list)
    # print(pop_list)
    # print(n_list)

    # format_image = hb.docker.build_python_image('gcr.io/ukbb-diversepops-neale/prs-csx-python',
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
                     P=args.P_col, target_pop=args.target_pop, outdir=args.out_dir)

    format_b.run()

    # format summstats string names for input
    sst_file_paths = hl.utils.hadoop_ls(f'{args.out_dir}/formatted_sst_files_for{args.target_pop}')
    sst_list2 = []
    sst_list_basenames = []
    for i in sst_file_paths:
        sst_list2.append(i['path'])
        basename = os.path.basename(i['path'])
        sst_list_basenames.append(basename)
    # order summstats string names based on original input
    sst_list_input_basenames = []
    for i in sst_list:
        basename = os.path.basename(i)
        root, extension = os.path.splitext(basename)
        sst_list_input_basenames.append(root)

    sst_list_basenames_sorted = sorted(sst_list_basenames,key=sst_list_input_basenames.index)

    # WILL CHANGE THE LINES ABOVE: ALL WE NEED IS THESE LAST TWO LINES, BUT TAKE THE FORMATTED_SST_PATH AND APPEND TO THE ORIGINAL SST_LIST AND THAT KEEPS THE ORDER
    formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{args.target_pop}/'
    final_sst_list = ['{}{}'.format(formatted_sst_path,i) for i in sst_list_basenames_sorted]

    print(final_sst_list)
    print(pop_list)

    b = hb.Batch(backend=backend, name='prscsx')
    prs_img = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscsx'

    # read in snp info file
    snp_info = b.read_input(args.snp_info)

    # get ref panels
    refpanels_dict = {}
    anc_list = ['afr', 'amr', 'eas', 'eur', 'sas']
    ref_file_sizes = 0
    for anc in anc_list:
        file = os.path.join(args.ref_path, f'ldblk_ukbb_{anc}.tar.gz')
        refpanels_dict[anc] = b.read_input(file)
        ref_size = bytes_to_gb(file)
        ref_file_sizes += round(10.0 + 2.0 * ref_size)

    refpanels_j = b.new_job(name=f'tar_refpanels')
    refpanels_j.cpu(8)
    tar_refpanels_job_storage = ref_file_sizes + 20
    refpanels_j.storage(f'{tar_refpanels_job_storage}Gi')
    refpanels_j.command(f'mkdir {refpanels_j.ref_panels}')
    refpanels_j.command(f'mv {snp_info} {refpanels_j.ref_panels}')

    for ancestry, input_file in refpanels_dict.items():
        refpanels_j.command(f'tar xf {input_file} -C {refpanels_j.ref_panels}')
    refpanels_j.command(f'ls {refpanels_j.ref_panels}')

    # run PRS-CSx
    prscsx_jobs = []
    for chrom in range(1,23):
        prscsx_j = run_prscsx(b=b, image=prs_img, depends_on_j=refpanels_j, bfile=args.bfile_path,
                              target_pop=args.target_pop, summary_stats=final_sst_list, N=n_list, pops=pop_list, chrom=chrom, ref_panels_dir=refpanels_j.ref_panels, meta=args.meta, refs_size=ref_file_sizes,
                              out_dir=args.out_dir)
        prscsx_jobs.append(prscsx_j)

    # run PLINK
    plink_img = 'hailgenetics/genetics:0.2.37'
    input_bfile = b.read_input_group(bed=f'{args.bfile_path}/{args.target_pop}/{args.target_pop}.bed',
                                     bim=f'{args.bfile_path}/{args.target_pop}/{args.target_pop}.bim',
                                     fam=f'{args.bfile_path}/{args.target_pop}/{args.target_pop}.fam')
    for pop in pop_list:
        bed_size = bytes_to_gb(f'{args.bfile_path}/{args.target_pop}/{args.target_pop}.bed')
        run_plink(b=b, depends_on_j_list=prscsx_jobs, image=plink_img, bfile_size=bed_size, target_pop=args.target_pop,
                  bfile=input_bfile, source_pop=pop, dup_ids_file=args.dup_ids, out_dir=args.out_dir)

    b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True, help='path to tab-delimited file with first column path to summstats, second column POP str, third column sample size')
    parser.add_argument('--SNP_col', required=True, help='name of rsID column in summstats')
    parser.add_argument('--A1_col', required=True, help='name of effect allele column in summstats')
    parser.add_argument('--A2_col', required=True, help='name of non-effect allele column in summstats')
    parser.add_argument('--A1_BETA_col', required=True, help='name of beta column in summstats')
    parser.add_argument('--P_col', required=True, help='name of P-value column in summstats')
    parser.add_argument('--bfile_path', required=True, help='path to bfile of target population, bfiles need to be in format "gs://path/target_pop_str/target_pop_str.bed,bim,fam" but input "gs://path" here')
    parser.add_argument('--target_pop', required=True, help='POP str of target population')
    parser.add_argument('--ref_path', required=True, help='path to reference panels directory for PRS-CSx')
    parser.add_argument('--snp_info', required=True, help='full path to SNP info file, including filename')
    parser.add_argument('--dup_ids', help='full path to file with list of duplicated SNP IDs to exclude, including filename')
    parser.add_argument('--out_dir', required=True, help='path to output directory')
    parser.add_argument('--meta', action="store_true", help='include flag if you want PRS-CSx meta-analysis option')
    arguments = parser.parse_args()

    main(arguments)