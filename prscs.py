import hailtop.batch as hb
import os
from hailtop.batch.job import Job
from format_input import bytes_to_gb
from typing import Dict, Union


def prscs(b: hb.batch.Batch = None,
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


def run_prscs(args,
              backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
              infile_dict: Dict = None):

    print('2. RUNNING PRS-CS')
    prscs_b = hb.Batch(backend=backend, name='run-prsCS')

    for i in range(0, len(infile_dict['phenos'])):
        pheno = infile_dict['phenos'][i]
        refpanel_pop = infile_dict['refpanel_pops'][i]
        target_cohort = infile_dict['target_cohorts'][i]
        n = infile_dict['ns'][i]
        sst = infile_dict['ssts'][i]

        formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{pheno}/'
        basename = os.path.basename(sst)
        root, extension = os.path.splitext(basename)
        final_sst = f'{formatted_sst_path}{root}'

        # read in ref panel
        ref_filename = os.path.join(args.ref_path, f'ldblk_ukbb_{refpanel_pop}.tar.gz')
        ref_panel = prscs_b.read_input(ref_filename)
        ref_size = bytes_to_gb(ref_filename)
        ref_file_size = round(10.0 + 2.0 * ref_size)

        refpanel_j = prscs_b.new_job(name=f'tar_refpanel')
        refpanel_j.cpu(4)
        tar_refpanel_job_storage = ref_file_size + 20
        refpanel_j.storage(f'{tar_refpanel_job_storage}Gi')
        refpanel_j.command(f'mkdir {refpanel_j.ref_panel}')
        refpanel_j.command(f'tar -zxvf {ref_panel} -C {refpanel_j.ref_panel}')
        refpanel_j.command(f'ls {refpanel_j.ref_panel}')

        refpanel_dir = f'{refpanel_j.ref_panel}/ldblk_ukbb_{refpanel_pop}'

        # run PRS-CS
        for chrom in range(21, 23):
            prscs(b=prscs_b, depends_on_j=refpanel_j, bfile=args.bfile_path, pheno=pheno, target_cohort=target_cohort,
                  summary_stats=final_sst, N=n, chrom=chrom, ref_panels_dir=refpanel_dir, refs_size=ref_file_size,
                  out_dir=args.out_dir)

    prscs_b.run()
