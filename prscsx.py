import hailtop.batch as hb
import os
from hailtop.batch.job import Job
from format_input import bytes_to_gb
from typing import List, Dict, Union


def prscsx(b: hb.batch.Batch,
           depends_on_j: Job = None,
           bfile: str = None,
           target_pop: str = None,
           summary_stats: List = None,
           N: List = None,
           pops: List = None,
           chrom: int = None,
           ref_panels_dir: str = None,
           refs_size: int = None,
           out_dir: str = None,
           pheno: str = None,
           meta = None,
           image: str = 'gcr.io/ukbb-diversepops-neale/ktsuo-prscsx'):
    """
    Run PRS-CSx using formatted summary statistics
    :return:
    """

    j = b.new_job(name=f'run-prscsx-{chrom}')
    j.depends_on(depends_on_j)

    # get bfile
    input_bim_file = b.read_input_group(bim=f'{bfile}/{target_pop}.bim')

    sst_files_list = []
    sst_files_sizes = 0

    for file in summary_stats:
        input_sst = b.read_input(file)
        # TRIAL 8/1
        j.command(f'head -n10 {input_sst}')
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
    b.write_output(j.scores, f'{out_dir}/{pheno}/prs_csx_output_for{target_pop}')


def run_prscsx(args,
              backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
              infile_dict: Dict = None):

    print('2. RUNNING PRS-CSX')
    prscsx_b = hb.Batch(backend=backend, name='run-prsCSX')

    for i in range(0, len(infile_dict['phenos'])):
        pheno = infile_dict['phenos'][i]
        sst_list = "".join(infile_dict['ssts'][i]).split(",")
        target_pop = infile_dict['target_pops'][i]

        sst_list_basenames = []
        for sst_i in range(0, len(sst_list)):
            basename = os.path.basename(sst_list[sst_i])
            root, extension = os.path.splitext(basename)
            sst_list_basenames.append(root)

        # formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{target_pop}/'
        formatted_sst_path = f'{args.out_dir}/formatted_sst_files_for{pheno}/'
        final_sst_list = ['{}{}'.format(formatted_sst_path, i) for i in sst_list_basenames]

        gwas_pops_list = "".join(infile_dict['gwas_pops'][i]).split(",")
        ns_list = "".join(infile_dict['ns'][i]).split(",")

        # read in snp info file
        snp_info = prscsx_b.read_input(args.snp_info)

        # get ref panels
        refpanels_dict = {}
        anc_list = ['afr', 'amr', 'eas', 'eur', 'sas']
        ref_file_sizes = 0
        for anc in anc_list:
            file = os.path.join(args.ref_path, f'ldblk_ukbb_{anc}.tar.gz')
            refpanels_dict[anc] = prscsx_b.read_input(file)
            ref_size = bytes_to_gb(file)
            ref_file_sizes += round(10.0 + 2.0 * ref_size)

        refpanels_j = prscsx_b.new_job(name=f'tar_refpanels')
        refpanels_j.cpu(8)
        tar_refpanels_job_storage = ref_file_sizes + 20
        refpanels_j.storage(f'{tar_refpanels_job_storage}Gi')
        refpanels_j.command(f'mkdir {refpanels_j.ref_panels}')
        refpanels_j.command(f'mv {snp_info} {refpanels_j.ref_panels}')

        for ancestry, input_file in refpanels_dict.items():
            refpanels_j.command(f'tar xf {input_file} -C {refpanels_j.ref_panels}')
        refpanels_j.command(f'ls {refpanels_j.ref_panels}')

        # run PRS-CSx
        for chrom in range(1, 23):
            prscsx(b=prscsx_b, depends_on_j=refpanels_j, bfile=args.bfile_path, target_pop=target_pop,
                   summary_stats=final_sst_list, N=ns_list, pops=gwas_pops_list, chrom=chrom,
                   ref_panels_dir=refpanels_j.ref_panels, meta=args.meta, refs_size=ref_file_sizes,
                   out_dir=args.out_dir, pheno=pheno)

    prscsx_b.run()
