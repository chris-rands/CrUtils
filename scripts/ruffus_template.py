'''
A minimalist Python ruffus (http://ruffus.readthedocs.io/en/latest/)
pipeline template.

This script could be executed as:
python ruffus_template.py --verbose 6 --target_tasks test_function2

Note this assumes you have ruffus version >=2.4, `ruffus.cmdline` does not exist
in earlier versions and a manual work around
(https://stackoverflow.com/a/41834283/6260170) is required.
'''

from ruffus import *

parser = cmdline.get_argparse(description='Example pipeline')
options = parser.parse_args()

@originate('test_out.txt')
def test_function(out_file):
        with open(out_file,'w') as f:
            f.write('it is working!\n')

@follows(test_function)
@transform(test_function, suffix('_out.txt'), '_out2.txt')
def test_function2(in_file, out_file):
    with open(in_file) as in_f, open(out_file,'w') as out_f:
        for line in in_f:
                out_f.write(line)    
        out_f.write('next file generated\n')

cmdline.run(options)
