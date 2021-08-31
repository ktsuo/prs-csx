import hail as hl
import argparse
import hailtop.batch as hb
import pandas as pd

def format_input(filepath, SNP, A1, A2, BETA, P):
    if hl.hadoop_exists(filepath):
        
        basename = os.path.basename(filepath)
        root, extension = os.path.splitext(basename)
        
        cols = [SNP, A1, A2, BETA, P]
        summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)[cols]
        summstats = summstats.rename(columns={SNP:'SNP', A1:'A1', A2:'A2', BETA:'BETA', P:'P'})
                
        outfile_name = f'{root}_formatted.txt'
        summstats.to_csv(outfile_name, sep='\t')

        return outfile_name
    else:
        print('summstats file not found!')
     
def run_prscsx(b: hb.batch.Batch,
        image: str,
        ref_dir: str, 
        bfile: hb.resource.ResourceFile, 
        summary_stats: list,
        N: list,
        pops: list,
        chr:int,
        meta,
        out_name: str,
        out_dir: str):

    j = b.new_job(name=f'run-prscsx-{chr}')
    
    j.image(image)
    j.cpu(4)
        
    j.command(f'''python3.8 PRS-CSx.py \
        --ref_dir={ref_dir} \ # this has to be a path to the reference panels
        --bim_prefix={bfile} \ # this has to be a path to the bim file plsu the bim prefix (without .bim)
        --sst_file={summary_stats} \
        --n_gwas={N} \
        --pop={pops} \
        --out_dir={out_dir} \
        --out_name={out_name} \
        --chrom={chr} \ 
        --meta={meta}''')

    return j

def main(args):
    backend = hb.ServiceBackend(billing_project='', bucket='')

    sst_list = []
    pop_list = []
    n_list = []
    with open(args.input_file, 'rt') as f:
        rows = f.readlines()
    
        for row in rows:   ### pheno col?
            sst_path, pop, n = row.split('\t')
            sst_list.append(sst_path)
            pop_list.append(pop)
            n_list.append(n)

    b = hb.Batch(backend=backend, name='prscsx')
    #jobs = []
    format_sst_output = [] 
    for sst in sst_list:
        j = b.new_python_job(name=f'formatting-{sst}')
        res = j.call(format_input, sst, args.SNP_col, args.A1_col, args.A2_col, args.A1_BETA_col, args.P_col) 
        #jobs.append(j)
        format_sst_output.append(res)
        # format summstats string names for input 
    sst_list = ','.join(format_sst_output)
    
    # format populations for input
    pop_list = ','.join(pop_list)

    # format sample sizes for input
    N_list = ",".join(map(str, n_list))

    image='gcr.io/ukbb-diversepops-neale/ktsuo-prscsx'

    ref_dir = b.read_input_group(**{'tar.gz': args.ref_path})
    open_refdir = b.new_job(name='open-ref-panels')
    open_refdir.command(f'tar -zwvf {ref_dir}')
    open_refdir.run()

    input_bfile = b.read_input_group(bed=f'{args.bfile_path}.bed', bim=f'{args.bfile_path}.bim', fam=f'{args.bfile_path}.fam')

    out_dir='tmp/hail/dir' # how to specify an out dir in hail batch environment?
    out_name='kristin_trial'

    chrom_list=list(map(str,range(1,22)))

    for chr in chrom_list:
        
        j = run_prscsx(b, image, refdir, input_bfile, sst_list, N_list, pop_list, chr, out_name, out_dir, args.meta)

    b.run()
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--SNP_col', required=True)
    parser.add_argument('--A1_col', required=True)
    parser.add_argument('--A2_col', required=True)
    parser.add_argument('--A1_BETA_col', required=True)
    parser.add_argument('--P_col', required=True)
    parser.add_argument('--bfile_path')
    parser.add_argument('--ref_path', required=True) 
    parser.add_argument('--meta', action="store_true")
                
    args = parser.parse_args()

    main(args)


