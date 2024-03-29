import hail as hl
import argparse
import hailtop.batch as hb
import pandas as pd
import os
from hailtop.batch.job import Job
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

    print(f'formatting {filepath}')
    basename = os.path.basename(filepath)
    root, extension = os.path.splitext(basename)

    cols = [SNP, A1, A2, BETA, P]
    summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)
    summstats = summstats.rename(columns={SNP: 'SNP', A1: 'A1', A2: 'A2', BETA: 'BETA', P: 'P'}, inplace=False)
    new_order = ['SNP', 'A1', 'A2', 'BETA', 'P']
    df = summstats[new_order]

    snp_info = pd.read_table(snp_info_filepath, header=0, sep='\t', usecols=['SNP'])
    df_filtered = df.merge(snp_info, how='inner', on='SNP')

    #outfile_name = f'{root}_formatted.txt'
    outfile_name = f'{root}'
    #TODO: for prs-csx formatted files are in '{outdir}/formatted_sst_files_for{target_pop}/
    df_filtered.to_csv(f'{outdir}/formatted_sst_files_for{pheno}/{outfile_name}', sep='\t', index=False)


def run_prscs(b: hb.batch.Batch = None,
          depends_on_j: Job = None,
          bfile: str = None,
          target_cohort: str = None,
          pheno: str = None,
          summary_stats: str = None,
          N: int = None,
          chrom: int = None,
          ref_panels_dir: str = None,
          refs_size: int = None,
          out_dir: str = None,
          image: str = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscs'):
    """
    Run PRS-CS using formatted summary statistics
    :param b: Batch object to add jobs to
    :param depends_on_j: job that the created job should only run after
    :param image: Batch object to add jobs to
    :param bfile: Batch object to add jobs to
    :param target_cohort: Batch object to add jobs to
    :param pheno: job that the created jobs should only run after
    :param summary_stats: Batch object to add jobs to
    :param N: Batch object to add jobs to
    :param chrom: job that the created jobs should only run after
    :param ref_panels_dir: Batch object to add jobs to
    :param refs_size: Batch object to add jobs to
    :param out_dir: directory path to write output files to
    :return: PRS-CS output files
    """

    j = b.new_job(name=f'run-prscs-{chrom}')
    j.depends_on(depends_on_j)

    # get bfile
    input_bim_file = b.read_input_group(bim=f'{bfile}/{target_cohort}.bim')

    input_sst = b.read_input(summary_stats)
    sst_file_size = round(10.0 + 2.0 * bytes_to_gb(summary_stats))

    job_storage = refs_size + sst_file_size + 20

    j.image(image)
    j.cpu(8)
    j.storage(f'{job_storage}Gi')
    # j.memory('highmem')
    j.command('mkdir tmp_prscs_output')

    j.command(f'''
    python3 PRScs.py \
        --ref_dir={ref_panels_dir} \
        --bim_prefix={input_bim_file} \
        --sst_file={input_sst} \
        --n_gwas={N} \
        --out_dir=tmp_prscs_output/tmp \
        --chrom={chrom}''')

    j.command(f'mv tmp_prscs_output {j.scores}')
    b.write_output(j.scores, f'{out_dir}/{pheno}/prs_cs_output_for{target_cohort}')

    return j


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
               pheno: str,
               meta):
    """
    Run PRS-CSx using formatted summary statistics
    :return: PRS-CSx output files
    """

    j = b.new_job(name=f'run-prscsx-{chrom}')
    j.depends_on(depends_on_j)

    # get bfile
    input_bim_file = b.read_input_group(bim=f'{bfile}/{target_pop}/{target_pop}.bim')

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
    # j.memory('highmem')
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
    b.write_output(j.scores, f'{out_dir}/{pheno}/prs_csx_output_for{target_pop}')

    return j


# Changed from 'sum center' to 'sum' only
def run_plink(b: hb.batch.Batch,
              depends_on_j_list: list,
              bfile_size: int,
              image: str,
              bfile: hb.ResourceGroup,
              target: str, # either target_cohort for PRS-CS or target_pop for PRS-CSx
              pheno: str,
              dup_ids_file: str,
              out_dir: str,
              prs_method: str,
              source_pop = "None",
              meta = "None"): # optional argument depending on whether using PRS-CSx
    """
    Using PRS-CS(x) output compute PRS for target population from each input summary statistic in PLINK
    :return: PLINK score files
    """

    if prs_method == 'prscs':
        input_files = hl.hadoop_ls(f'{out_dir}/{pheno}/prs_cs_output_for{target}/tmp_pst_eff_a1_b0.5_phiauto_chr*')

        j = b.new_job(name=f'run-plink')
        j.depends_on(*depends_on_j_list)
        j.image(image)
        j.cpu(4)

        j.memory('highmem')

        job_storage = bfile_size + 20
        j.storage(f'{job_storage}Gi')

        j.command('mkdir tmp_input')
        for file in input_files:
            file = b.read_input(file['path'])
            j.command(f'mv {file} tmp_input')
        j.command('ls tmp_input')

        j.command('cat tmp_input/tmp_pst_eff_a1_b0.5_phiauto_chr* > tmp_input/tmp_pst_eff_a1_b0.5_phiauto_ALLchr.txt')

        j.command('mkdir tmp_output')
        # exclude duplicate SNPs
        if dup_ids_file:
            dup_ids_input = b.read_input(dup_ids_file)
            j.command(f'plink \
                --bfile {bfile} \
                    --score tmp_input/tmp_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                        --exclude {dup_ids_input} \
                            --out tmp_output/{target}_scores')
            j.command(f'mv tmp_output {j.output}')
            b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')

        if dup_ids_file is None:
            j.command(f'plink \
                --bfile {bfile} \
                    --score tmp_input/tmp_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                        --out tmp_output/{target}_scores')
            j.command(f'mv tmp_output {j.output}')
            b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')

    if prs_method == 'prscsx':
        if meta:
            input_files = hl.hadoop_ls(f'{out_dir}/{pheno}/prs_csx_output_for{target}/tmp_META_pst_eff_a1_b0.5_phiauto_chr*') 
            j = b.new_job(name=f'run-plink-prscsx-meta')
            j.depends_on(*depends_on_j_list)
            j.image(image)
            j.cpu(4)

            j.memory('highmem')

            job_storage = bfile_size + 20
            j.storage(f'{job_storage}Gi')

            j.command('mkdir tmp_input')
            for file in input_files:
                file = b.read_input(file['path'])
                j.command(f'mv {file} tmp_input')
            j.command('ls tmp_input')

            j.command(f'cat tmp_input/tmp_META_pst_eff_a1_b0.5_phiauto_chr* > tmp_input/tmp_META_pst_eff_a1_b0.5_phiauto_ALLchr.txt')

            j.command('mkdir tmp_output')
            # exclude duplicate SNPs
            if dup_ids_file:
                dup_ids_input = b.read_input(dup_ids_file)
                j.command(f'plink \
                    --bfile {bfile} \
                        --score tmp_input/tmp_META_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                            --exclude {dup_ids_input} \
                                --out tmp_output/from_META')
                j.command(f'mv tmp_output {j.output}')
                b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')

            if dup_ids_file is None:
                j.command(f'plink \
                    --bfile {bfile} \
                        --score tmp_input/tmp_META_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                            --out tmp_output/from_META')
                j.command(f'mv tmp_output {j.output}')
                b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')

        else:
            input_files = hl.hadoop_ls(f'{out_dir}/{pheno}/prs_csx_output_for{target}/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_chr*')

            j = b.new_job(name=f'run-plink-{source_pop}')
            j.depends_on(*depends_on_j_list)
            j.image(image)
            j.cpu(4)

            j.memory('highmem')

            job_storage = bfile_size + 20
            j.storage(f'{job_storage}Gi')

            j.command('mkdir tmp_input')
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
                        --score tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                            --exclude {dup_ids_input} \
                                --out tmp_output/from_{source_pop}')
                j.command(f'mv tmp_output {j.output}')
                b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')

            if dup_ids_file is None:
                j.command(f'plink \
                    --bfile {bfile} \
                        --score tmp_input/tmp_{source_pop}_pst_eff_a1_b0.5_phiauto_ALLchr.txt 2 4 6 sum \
                            --out tmp_output/from_{source_pop}')
                j.command(f'mv tmp_output {j.output}')
                b.write_output(j.output, f'{out_dir}/{pheno}/{target}_plink_scores')


def main(args):
    global prscsx_j
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket='ukb-diverse-pops')

    file = pd.read_csv(args.input_file)

    if args.prs_cs:
        phenos, ssts, refpanel_pops, ns, SNP_names, A1_names, A2_names, beta_names, pval_names, target_cohorts = ([] for i in range(10))

        for index, row in file.iterrows():
            phenos.append(row["phenos"])
            ssts.append(row["sum_stats_paths"])  
            refpanel_pops.append(row["populations"]) # either one refpanel_pop for PRS-CS or multiple populations corresponding to input GWAS for PRS-CSx
            ns.append(row["sample_sizes"])
            SNP_names.append(row["snp_names"])
            A1_names.append(row["A1_names"]) 
            A2_names.append(row["A2_names"])
            beta_names.append(row["Beta_names"])
            pval_names.append(row["p-value_label_names"])
            target_cohorts.append(row["target"])

        if args.format:
            format_b = hb.Batch(backend=backend, name='prscs_x-formatting',
                default_python_image='gcr.io/ukbb-diversepops-neale/prs-csx-python')
            for i in range(0, len(phenos)):
                print(f'Analyzing {phenos[i]}')
                pheno = phenos[i]
                sst = ssts[i]
                refpanel_pop = refpanel_pops[i]
                n = ns[i]
                SNP_name = SNP_names[i]
                A1_name = A1_names[i]
                A2_name = A2_names[i]
                beta = beta_names[i]
                pval = pval_names[i]
                target_cohort = target_cohorts[i]

                j = format_b.new_python_job(name=f'Formatting: {sst}')
                sst_size = bytes_to_gb(sst)
                disk_size = round(4.0 + 2.0 * sst_size)
                j.storage(disk_size)
                j.cpu(4)
                j.call(format_input, sst, SNP_name, A1_name, A2_name, beta, pval, pheno, args.out_dir, args.snp_info)            

            format_b.run()
        else:
            print('Summstats already formatted!')

        prscs_b = hb.Batch(backend=backend, name='run-prsCS')
        for i in range(0, len(phenos)):
            pheno = phenos[i]
            refpanel_pop = refpanel_pops[i]
            target_cohort = target_cohorts[i]
            n = ns[i]
            sst = ssts[i]

            formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{pheno}/'
            basename = os.path.basename(sst)
            root, extension = os.path.splitext(basename)
            final_sst = f'{formatted_sst_path}{root}'

            prs_img = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscs' 

            # read in ref panel
            ref_filename = os.path.join(args.ref_path, f'ldblk_ukbb_{refpanel_pop}.tar.gz')
            ref_panel = prscs_b.read_input(ref_filename)
            ref_size = bytes_to_gb(ref_filename)
            ref_file_size = round(10.0 + 2.0 * ref_size)

            refpanel_j = b.new_job(name=f'tar_refpanel')
            refpanel_j.cpu(4)
            tar_refpanel_job_storage = ref_file_size + 20
            refpanel_j.storage(f'{tar_refpanel_job_storage}Gi')
            refpanel_j.command(f'mkdir {refpanel_j.ref_panel}')
            refpanel_j.command(f'tar -zxvf {ref_panel} -C {refpanel_j.ref_panel}')
            refpanel_j.command(f'ls {refpanel_j.ref_panel}')

            refpanel_dir = f'{refpanel_j.ref_panel}/ldblk_ukbb_{refpanel_pop}'

            # run PRS-CS
            prscs_jobs = []
            for chrom in range(1,23):
                prscs_j = run_prscs(b=prscs_b, image=prs_img, depends_on_j=refpanel_j, bfile=args.bfile_path,
                                    pheno=pheno, target_cohort=target_cohort, summary_stats=final_sst, N=n, chrom=chrom, ref_panels_dir=refpanel_dir, refs_size=ref_file_size,
                                    out_dir=args.out_dir)
                prscs_jobs.append(prscs_j)

            # run PLINK
            plink_img = 'hailgenetics/genetics:0.2.37'
            input_bfile = b.read_input_group(bed=f'{args.bfile_path}/{target_pop}/{target_cohort}.bed',
                                            bim=f'{args.bfile_path}/{target_pop}/{target_cohort}.bim',
                                            fam=f'{args.bfile_path}/{target_pop}/{target_cohort}.fam')

            bed_size = bytes_to_gb(f'{args.bfile_path}/{target_cohort}.bed')
            run_plink(b=b, depends_on_j_list=prscs_jobs, image=plink_img, bfile_size=bed_size,
                    bfile=input_bfile, target=target_cohort, pheno=pheno, dup_ids_file=args.dup_ids, out_dir=args.out_dir, prs_method="prscs")
            # run_plink(b=b, image=plink_img, bfile_size=bed_size,
                    # bfile=input_bfile, target=target_cohort, pheno=pheno, dup_ids_file=args.dup_ids, out_dir=args.out_dir, prs_method="prscs")
        prscs_b.run()

        prscs_plink_b = hb.Batch(backend=backend, name='run-plink-prsCS')
        for i in range(0, len(phenos)):
            pheno = phenos[i]
            refpanel_pop = refpanel_pops[i]
            target_cohort = target_cohorts[i]
            n = ns[i]
            sst = ssts[i]

            # run PLINK
            plink_img = 'hailgenetics/genetics:0.2.37'
            input_bfile = prscs_plink_b.read_input_group(bed=f'{args.bfile_path}/{target_pop}/{target_cohort}.bed',
                                             bim=f'{args.bfile_path}/{target_pop}/{target_cohort}.bim',
                                             fam=f'{args.bfile_path}/{target_pop}/{target_cohort}.fam')

            bed_size = bytes_to_gb(f'{args.bfile_path}/{target_cohort}.bed')
            run_plink(b=prscs_plink_b, depends_on_j_list=prscs_jobs, image=plink_img, bfile_size=bed_size,
                      bfile=input_bfile, target=target_cohort, pheno=pheno, dup_ids_file=args.dup_ids,
                      out_dir=args.out_dir, prs_method="prscs")

        prscs_plink_b.run()

    if args.prs_csx:
        phenos, ssts, gwas_pops, ns, SNP_names, A1_names, A2_names, beta_names, pval_names, target_pops, dup_rsid_files = ([] for i in range(11))

        for index, row in file.iterrows():
            phenos.append(row["phenos"])
            ssts.append([row["sum_stats_paths"]])  
            gwas_pops.append([row["populations"]]) # either one refpanel_pop for PRS-CS or multiple populations corresponding to input GWAS for PRS-CSx
            ns.append([row["sample_sizes"]])
            SNP_names.append([row["snp_names"]])
            A1_names.append([row["A1_names"]]) 
            A2_names.append([row["A2_names"]])
            beta_names.append([row["Beta_names"]])
            pval_names.append([row["p-value_label_names"]])
            target_pops.append(row["target"])
            dup_rsid_files.append(row["dup_ids_file"])

        if args.format:
            format_b = hb.Batch(backend=backend, name='prscs_x-formatting', default_python_image='gcr.io/ukbb-diversepops-neale/prs-csx-python')
            for i in range(0, len(phenos)):
                print(f'Analyzing {phenos[i]}')
                pheno = phenos[i]
                sst_list = "".join(ssts[i]).split(",")
                SNP_names_list = "".join(SNP_names[i]).split(",")
                A1_names_list = "".join(A1_names[i]).split(",")
                A2_names_list = "".join(A2_names[i]).split(",")
                beta_names_list = "".join(beta_names[i]).split(",")
                pval_names_list = "".join(pval_names[i]).split(",")
                
                for sst_i in range(0, len(sst_list)):
                    print(sst_list[sst_i], SNP_names_list[sst_i], A1_names_list[sst_i], A2_names_list[sst_i], beta_names_list[sst_i],
                                pval_names_list[sst_i])
                    j = format_b.new_python_job(name=f'Formatting: {sst_list[sst_i]}')
                    sst_size = bytes_to_gb(sst_list[sst_i])
                    disk_size = round(4.0 + 2.0 * sst_size)
                    j.storage(disk_size)
                    j.cpu(4)
                    j.call(format_input, sst_list[sst_i], SNP_names_list[sst_i], A1_names_list[sst_i], A2_names_list[sst_i], beta_names_list[sst_i],
                                pval_names_list[sst_i], pheno, args.out_dir) # instead of creating a folder of formatted summstats for each target pop, create a folder for each pheno

            format_b.run()
        else:
            print('Summstats already formatted!')

        b = hb.Batch(backend=backend, name='prscsx')
 
        for i in range(0, len(phenos)):
            pheno = phenos[i]
            sst_list = "".join(ssts[i]).split(",")
            target_pop = target_pops[i]
    
            sst_list_basenames = []
            for sst_i in range(0, len(sst_list)):
                basename = os.path.basename(sst_list[sst_i])
                root, extension = os.path.splitext(basename)
                sst_list_basenames.append(root)

            formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{target_pop}/'
            final_sst_list = ['{}{}'.format(formatted_sst_path,i) for i in sst_list_basenames]

            gwas_pops_list = "".join(gwas_pops[i]).split(",")
            ns_list = "".join(ns[i]).split(",")
            target_pop = target_pops[i]

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
            prscsx_jobs = [] # this is j.scores
            for chrom in range(1,23):
                prscsx_j = run_prscsx(b=b, image=prs_img, depends_on_j=refpanels_j, bfile=args.bfile_path,
                                    target_pop=target_pop, summary_stats=final_sst_list, N=ns_list, pops=gwas_pops_list, chrom=chrom, ref_panels_dir=refpanels_j.ref_panels, 
                                    meta=args.meta, refs_size=ref_file_sizes, out_dir=args.out_dir, pheno=pheno)
                prscsx_jobs.append(prscsx_j)

            # run PLINK
            # TODO: create another batch and get rid of depends_on_j_list
            b = hb.Batch(backend=backend, name='prscsx')

            plink_img = 'hailgenetics/genetics:0.2.37'
            input_bfile = b.read_input_group(bed=f'{args.bfile_path}/{target_pop}/{target_pop}.bed',
                                            bim=f'{args.bfile_path}/{target_pop}/{target_pop}.bim',
                                            fam=f'{args.bfile_path}/{target_pop}/{target_pop}.fam')
            
            if args.meta:
                bed_size = bytes_to_gb(f'{args.bfile_path}/{target_pop}/{target_pop}.bed')
                # TODO: submit input_files to run_plink() instead of inside the function 
                input_files = hl.hadoop_ls(f'{args.out_dir}/{pheno}/prs_csx_output_for{target_pop}/tmp_META_pst_eff_a1_b0.5_phiauto_chr*') 

                run_plink(b=b, depends_on_j_list=prscsx_jobs, image=plink_img, bfile_size=bed_size, target=target_pop,
                        bfile=input_bfile, dup_ids_file=dup_rsid_files[i], out_dir=args.out_dir, prs_method='prscsx', pheno=pheno,
                        meta=args.meta)
                # run_plink(b=b, image=plink_img, bfile_size=bed_size, target=target_pop,
                #         bfile=input_bfile, dup_ids_file=dup_rsid_files[i], out_dir=args.out_dir, prs_method='prscsx', pheno=pheno,
                #         meta=args.meta)
                        
            else:
                for pop in gwas_pops_list:
                    bed_size = bytes_to_gb(f'{args.bfile_path}/{target_pop}/{target_pop}.bed')
                    run_plink(b=b, depends_on_j_list=prscsx_jobs, image=plink_img, bfile_size=bed_size, target=target_pop,
                            bfile=input_bfile, dup_ids_file=dup_rsid_files[i], out_dir=args.out_dir, prs_method='prscsx', pheno=pheno, source_pop=pop)
                    # run_plink(b=b, image=plink_img, bfile_size=bed_size, target=target_pop,
                    #         bfile=input_bfile, dup_ids_file=dup_rsid_files[i], out_dir=args.out_dir, prs_method='prscsx', pheno=pheno, source_pop=pop)

        b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True, help='path to csv-delimited file - SEE EXAMPLE')
    parser.add_argument('--prs_cs', action="store_true", help='include flag if running PRS-CS')
    parser.add_argument('--prs_csx', action="store_true", help='include flag if running PRS-CSx')
    parser.add_argument('--bfile_path', required=True, help='path to bfile of target population, bfiles need to be in format "gs://path/target_cohort_str.bed,bim,fam for PRS-CS or "gs://path/target_pop_str.bed,bim,fam for PRS-CSx; input "gs://path" here')
    parser.add_argument('--ref_path', required=True, help='path to reference panels directory for PRS-CSx')
    parser.add_argument('--snp_info', required=True, help='full path to SNP info file, including filename')
    # parser.add_argument('--dup_ids', help='full path to file with list of duplicated SNP IDs to exclude, including filename')
    parser.add_argument('--out_dir', required=True, help='path to output directory')
    parser.add_argument('--meta', action="store_true", help='include flag if you want PRS-CSx meta-analysis option')
    parser.add_argument('--format', action="store_true", help='include flag if you need to format summstats for PRS-CS(x) first')

    arguments = parser.parse_args()

    main(arguments)
