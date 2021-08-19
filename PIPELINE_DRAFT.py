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
    
        for row in rows:
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


    sink = b.new_job(name=f'run-prscsx')
    sink.call(run_prscsx) # run prscsx function here?


    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--SNP_col', required=True)
    parser.add_argument('--A1_col', required=True)
    parser.add_argument('--A2_col', required=True)
    parser.add_argument('--A1_BETA_col', required=True)
    parser.add_argument('--P_col', required=True)
    args = parser.parse_args()

    main(args)


