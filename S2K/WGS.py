import sys
import logging
from S2K import Genome

class WGS:
    """Class to handle WGS read counts file and create the genome."""
    def __init__ (self, wgs_file_name,  sample_name, parameters, models,   
                  output_filename,
                  no_processes = 1, verbosity = 'INFO', skip_filtering = False):        

        self.sample_name = sample_name
        self.output_filename = output_filename
        self.no_processes = no_processes
        self.wgs_file = open (wgs_file_name, 'r')
        self.config = parameters
        self.models = models
        

        self.logger = self.create_logger (verbosity)
        
        if skip_filtering:
            self.SG_file = None
            self.logger.debug ("Skipping SG filtering.")    
        else:
            try:
                self.SG_file = open (self.config['Input']['SuperGood_filepath'], 'rb')
            except FileNotFoundError:
                sys.exit(f"SuperGood_filepath: {self.config['Input']['SuperGood_filepath']} should be the full real path to supergood file. Exiting")

        try:
            self.CB_file = open (self.config['Input']['CytoBand_filepath'], 'r')
        except FileNotFoundError:
            sys.exit(f"CytoBand_filepath: {self.config['Input']['CytoBand_filepath']} should be the full real path to cytoband file. Exiting")
        
        self.logger.debug ("WGS object created.")
        
    def create_logger (self, verbosity):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(self.output_filename + '.log', mode = 'w')
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter('%(asctime)s %(name)s: %(levelname)s: %(message)s', datefmt='%H:%M:%S')
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)
        
        sh = logging.StreamHandler ()
        sh.setLevel (verbosity)
        sh.setFormatter (fh_formatter)
        logger.addHandler (sh)
        
        return logger
            
    def analyze (self, m0 = 0, mc = 5.0):
        self.logger.debug ('Creating genome.')
        self.genome = Genome.Genome (self.sample_name, self.logger, self.config, self.CB_file, 
                                     self.models, self.no_processes)
        input_columns = [self.config['InputColumns']['chrom'],
                         self.config['InputColumns']['position'],
                         self.config['InputColumns']['ref_count'],
                         self.config['InputColumns']['alt_count'],
                         self.config['InputColumns']['Type']]
        
        self.genome.retrive_counts_create_chromosomes (data_file = self.wgs_file, SG_file = self.SG_file,
                                                       columns = input_columns)
        self.logger.debug ('Segmenting genome...')
        self.genome.segment_genome (m0, mc)
        self.logger.info ('Ready to report!')
    
    def report (self, report_type = 'bed'):
        return self.genome.report(report_type)
