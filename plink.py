import hail as hl
import hailtop.batch as hb
from typing import Dict, Union


# Changed from 'sum center' to 'sum' only
def plink(b: hb.batch.Batch = None,
          bfile_size: int = None,
          image: str = 'docker.io/hailgenetics/genetics:0.2.67', # tried .37 last time
          bfile: hb.ResourceGroup = None,
          target: str = None,
          pheno: str = None,
          dup_ids_file: str = None,
          out_dir: str = None,
          prs_method: str = None,
          source_pop: str = "None",
          meta: bool = False):  # optional argument depending on whether using PRS-CSx
    """
    Using PRS-CS(x) output compute PRS for target population from each input summary statistic in PLINK
    :param b: Batch object to add jobs to
    :param bfile: target cohort or target pop PLINK files
    :param bfile_size: size of the PLINK files
    :param image: image to use for the job
    :param target: either target_cohort for PRS-CS or target_pop for PRS-CSx
    :param pheno: phenotype
    :param dup_ids_file: file with duplicate SNPs to be removed
    :param out_dir: output directory to write outputs to
    :param prs_method: PRS method used
    :param source_pop: either target_cohort for PRS-CS or target_pop for PRS-CSx. Optional
    :param meta:
    :return:
    """

    if prs_method == 'prscs':
        input_files = hl.hadoop_ls(f'{out_dir}/{pheno}/prs_cs_output_for{target}/tmp_pst_eff_a1_b0.5_phiauto_chr*')

        j = b.new_job(name=f'run-plink')
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


def run_plink(args,
              backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
              infile_dict: Dict = None, method: str = None):
    print('3. RUNNING PLINK')
    from format_input import bytes_to_gb

    if method == 'prscs':
        plink_prscs_b = hb.Batch(backend=backend, name='run-plink-prsCS')

        for i in range(0, len(infile_dict['phenos'])):
            pheno = infile_dict['phenos'][i]
            target_cohort = infile_dict['target_cohorts'][i]

            input_bfile = plink_prscs_b.read_input_group(bed=f'{args.bfile_path}/{target_cohort}.bed',
                                                         bim=f'{args.bfile_path}/{target_cohort}.bim',
                                                         fam=f'{args.bfile_path}/{target_cohort}.fam')

            bed_size = bytes_to_gb(f'{args.bfile_path}/{target_cohort}.bed')
            plink(b=plink_prscs_b, bfile_size=bed_size, bfile=input_bfile, target=target_cohort, pheno=pheno,
                  dup_ids_file=infile_dict['dup_rsid_files'][i], out_dir=args.out_dir, prs_method="prscs")

        plink_prscs_b.run()

    if method == 'prscsx':
        plink_prscsx_b = hb.Batch(backend=backend, name='run-plink-prsCSX')

        for i in range(0, len(infile_dict['phenos'])):
            pheno = infile_dict['phenos'][i]
            target_pop = infile_dict['target_pops'][i]
            gwas_pops_list = "".join(infile_dict['gwas_pops'][i]).split(",")

            input_bfile = plink_prscsx_b.read_input_group(bed=f'{args.bfile_path}/{target_pop}.bed',
                                                          bim=f'{args.bfile_path}/{target_pop}.bim',
                                                          fam=f'{args.bfile_path}/{target_pop}.fam')

            if args.meta:
                bed_size = bytes_to_gb(f'{args.bfile_path}/{target_pop}.bed')

                plink(b=plink_prscsx_b, bfile_size=bed_size, target=target_pop, bfile=input_bfile,
                      dup_ids_file=infile_dict['dup_rsid_files'][i], out_dir=args.out_dir, prs_method='prscsx',
                      pheno=pheno, meta=True)

            else:
                for pop in gwas_pops_list:
                    bed_size = bytes_to_gb(f'{args.bfile_path}/{target_pop}.bed')
                    plink(b=plink_prscsx_b, bfile_size=bed_size, target=target_pop, bfile=input_bfile,
                          dup_ids_file=infile_dict['dup_rsid_files'][i], out_dir=args.out_dir, prs_method='prscsx',
                          pheno=pheno, source_pop=pop)

        plink_prscsx_b.run()
