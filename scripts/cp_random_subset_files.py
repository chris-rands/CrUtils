import os
import random
import shutil
import argparse

def main(in_dir,out_dir,n):
    if not os.path.exists(in_dir):
        print('Input dir does not exist')
        raise SystemExit
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for file_name in random.sample(os.listdir(in_dir),n):
            print('Copying: {}/{} TO {}/{}'.format(in_dir, file_name, out_dir, file_name))
            shutil.copy2('{}/{}'.format(in_dir, file_name), '{}/{}'.format(out_dir, file_name))

if __name__ == "__main__":

    # Options
    parser = argparse.ArgumentParser(description='Randomly sub-sample files from a directory @Chris Rands 2016')
    parser.add_argument('-i', '--in_dir',required=True, help='Input directory')
    parser.add_argument('-o', '--out_dir', required=True, help='Output directory')
    parser.add_argument('-n', '--number', type=int, required=True, help='Number of files to sub-sample')
    args = parser.parse_args()

    # Run
    main(args.in_dir,args.out_dir,args.number)
        
    
