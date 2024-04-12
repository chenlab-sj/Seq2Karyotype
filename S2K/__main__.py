import argparse
import configparser
import os
import shutil
import subprocess
import sys

import pandas as pd
from S2K.Run import Solution

from S2K import WGS
from S2K import Models

_description = "Scan chromosomes in search for non-HE segments. Assigns copy numbers and karytotypes if can."


__version__ = '0.1.0'



def main():
    parser = argparse.ArgumentParser (prog="S2K", description = _description)
    subparsers = parser.add_subparsers(title="actions", dest="action")
    subparsers.required = True

    ### Analyze subparser ###
    parser_analyze = subparsers.add_parser("analyze", description="runs the analysis")
    parser_analyze.add_argument ('-s', '--sample_name', required = False, 
                                 type = str, default = '',
                                 help = 'Input sample name. Default: from file name.')
    parser_analyze.add_argument ('-n', '--no_processes', required = False, default = 1, type = int,
                                 help = 'Number of processes. Default: 1')
    parser_analyze.add_argument ('-i', '--input_file', required = True,
                                 type = str, default = '',
                                 help = "Name of the input file with alleles' counts")
    parser_analyze.add_argument ('-c', '--config', required = False, default = 'config.ini',
                                 help = 'INI file with parameters')
    parser_analyze.add_argument ('-l', '--level', default = 'INFO', 
                                 choices = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'NOTSET'],
                                 help = 'Level of verbosity for std.err logger.')
    parser_analyze.add_argument ('-r', '--report_solutions', help = 'Generate report with all solutions.',
                                 action = 'store_true')
    parser_analyze.add_argument ('-m0', '--coverage_diploid', required = False, type = float,
                                 help = 'Coverage of diploid.', default = 0)
    parser_analyze.add_argument ('-mc', '--merge_coeff', required = False, type = float,
                                 help = 'Coefficient that controls the merging process. This coefficient sets the threshold for merging segments: \
                                 a higher value facilitates segment merging by increasing the similarity threshold required for merging. \
                                 Essentially, it makes the merging process more lenient, allowing segments with greater differences to be combined. \
                                 If set to 0, no merging will be performed, preserving all segments as initially detected. The default value is set to 6.0.'\
                                 , default = 6.0)
    parser_analyze.add_argument ('-v', '--version', help = 'Print version', action = 'version',
                                 version = 'S2K v. {version}'.format(version = __version__))
    parser_analyze.add_argument ('-m', '--models', choices=Models.model_presets.keys(),
                                 nargs='+', help = 'Specify which of models should be included.',
                                 default = ['(AB)(2+n)','(AB)(2-n)','A','AA','AAB','AAAB','AAA','AAAA','A+AA','AAB+AAAB','AA+AAA','AA+AAB','AAB+AABB','AAA+AAAA'])

    parser_analyze.add_argument ('-sf', '--skip_filtering', help = 'Do not filter input data through SG list.',
                                  action = 'store_true')

    parser_analyze.set_defaults (func=analyze)

    ### Viewer subparser ###
    parser_viewer = subparsers.add_parser ("viewer", description="launches the viewer")
    parser_viewer.add_argument ('--remote', action="store_true",
                                help="use a if running from a remote machine, for example a compute cluster")
    parser_viewer.add_argument ('-p', '--port', type=str, help="specific port to use")
    parser_viewer.set_defaults (func=viewer)

    ### Get Config subparser ###
    get_config = subparsers.add_parser ("getconfig", description="copies default S2K config to current dir")
    get_config.add_argument ("-d", "--directory", default=os.getcwd(),
                             help="copies config to this dir. Default is current working directory.")
    get_config.set_defaults (func=get_s2k_config)

    args = parser.parse_args()
    args.func(args)

def analyze(args):
    """ runs the analysis """
    ini = configparser.ConfigParser ()
    ini.read (args.config)
    output_filename = args.sample_name + '_m{}_mc{:.0f}'.format(args.coverage_diploid, args.merge_coeff)
    
    sample = WGS.WGS (args.input_file,  sample_name = args.sample_name, parameters = ini,
                      no_processes = args.no_processes, models = args.models,
                      output_filename = output_filename,
                      verbosity = args.level, skip_filtering = args.skip_filtering)

    sample.analyze (m0 = args.coverage_diploid, mc = args.merge_coeff)
    

    with open (output_filename + '.bed', 'w') as bed:
        bed.writelines (sample.report(report_type = 'bed'))

    with open (output_filename + '.par', 'w') as params:
        params.writelines (sample.report (report_type = 'params'))

    if args.report_solutions:
        with open (output_filename + '.solutions', 'w') as full:
            full.writelines (sample.report(report_type = 'solution'))

    keys = sample.genome.chromosomes.keys()
    data = pd.concat ([sample.genome.chromosomes[k].data for k in keys])

    data.to_csv (output_filename + '.dat.gz', index = None, sep = '\t', 
                 compression = 'gzip')

    print ('All done')

def viewer(args):
    """ launches the viewer - we get a port not in use"""
    import socket
    s = socket.socket()
    hostname = socket.gethostname()
    host_to_use = socket.gethostbyname(hostname)
    s.bind(("", 0))
    open_port = str(s.getsockname()[1])
    s.close()
    if args.port:
        open_port = args.port
    if not args.remote:
        host_to_use = "localhost"
    cmd = ["shiny", "run", "--port", open_port, "--host", host_to_use, "S2K.viewer.app"]
    print("**********")
    print(f"Access dashboard in browser via: http://{host_to_use}:{open_port}")
    print("**********")

    proc = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)    

def get_s2k_config(args):
    """ local helper to copy config to current working directory """
    configdir = os.path.dirname (os.path.realpath(__file__))
    shutil.copyfile (f'{configdir}/S2K.ini', f'{args.directory}/S2K.ini')
    

if __name__ == '__main__':
    sys.exit(main())
