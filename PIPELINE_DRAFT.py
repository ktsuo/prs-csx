import hail as hl
import argparse
import hailtop.batch as hb
import pandas as pd

def format_input(filepath, SNP, A1, A2, BETA, P):
    if hl.hadoop_exists(filepath):
        
        basename = os.path.basename(filepath)
        root, extension = os.path.splitext(basename)
        
        cols = [SNP, A1, A2, BETA, P]
        summstats = pd.read_table(filepath, header=0, sep='\t', compression='gzip', usecols=cols)
        summstats = summstats.rename(columns={SNP:'SNP', A1:'A1', A2:'A2', BETA:'BETA', P:'P'})
                
        outfile_name = f'{root}_formatted.txt'
        summstats.to_csv(outfile_name, sep='\t')

        return outfile_name
    else:
        print('summstats file not found!')

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

 
## RUN PRS-Csx
def run_prs(batch: hb.batch.Batch, 
            image: str,
            ref_dir: str, 
            bim: hb.resource.ResourceFile, 
            summary_stats: List,
            N: int,
            chr:int,
            pops: List,
            out_name: str,
            out_dir: str)-> hb.job.Job:
    
    
   
    prs.image(image)
    prs.cpu(4)
            
    prs.command = f '" python3 {ref_dir}PRS-CSx.py \
        --ref_dir={ref_dir} \ #from input
        --bim_prefix={bim} \ #finput
        --sst_file={sst_list} \ #from above
        --n_gwas={n_list} \ #from above
        --pop={pop_list} \ #from above
        --out_dir={out_dir} \ #input
        --out_name={pheno} \ #from above
        --chrom={chrom} \ 
        --meta=TRUE \
        " '
    
    prs.depends_on(*jobs)
     
    for chr in range (1,22):chrom = f'chrom{chr}
    
    prs = batch.new_job(name = f'{pheno}-{chrom}-PRS-Csx')
    prs = run_prs(ref_dir, bim, out_dir)
            
    
    
    
    return prs


    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--SNP_col', required=True)
    parser.add_argument('--A1_col', required=True)
    parser.add_argument('--A2_col', required=True)
    parser.add_argument('--A1_BETA_col', required=True)
    parser.add_argument('--P_col', required=True)
    parser.add_argument('--chr', default='1-22')
    parser.add_argument('--ref_dir', required=True)
    parser.add_argument('--bim', required = True)
    parser.add_argument('--out_dir', required = True)                   
    args = parser.parse_args()

    main(args)


