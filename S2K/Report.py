from cmath import isnan
import numpy as np
import warnings as warn
import scipy.stats as sts

from S2K import Run
from S2K import Consts

class Report:
    """ Class that holds reports used by other objects in the program """
    def __init__(self, report_t):
        self._report_type = report_t

    def genome_report (self, genome):
        """ Generates a report for Genome objects """
        if self._report_type == 'bed':
            keys = list(genome.chromosomes.keys())

            keys.sort (key = Consts.CHROM_ORDER.index)
            report = '\n'.join([genome.chromosomes[key].report(report_type=self._report_type) for key in keys])
        elif self._report_type == 'params':
            report_list = ['m0\t'+ str(genome.genome_medians['m0']),
                           'l\t'+ str(genome.genome_medians['l']),
                           'v0\t'+str(genome.genome_medians['v0']),
                           'fb\t'+str(genome.genome_medians['fb']),
                           'm_ai\t'+str(genome.scorer.ai_param['m']),
                           's_ai\t'+str(genome.scorer.ai_param['s']),
                           'm_cn\t'+str(genome.scorer.cn_param['m']),
                           's_cn\t'+str(genome.scorer.cn_param['s']),
                           'm_d\t' +str(genome.scorer.dipl_dist['m']),
                           's_d\t' +str(genome.scorer.dipl_dist['s']),
                           'a_d\t' +str(genome.genome_medians['d_model']['a']),
                           'models\t'+str(genome.models),
                           'thr_model\t'+str(genome.genome_medians['thr_model']),
                           'thr_HE\t'+str(genome.genome_medians['thr_HE'])]
                        
            report = '\n'.join(report_list)
        else:
            report = ""
        return report

    def chromosome_report(self, segments):
        """ Generates a report for Chromosome objects """
        data = '\n'.join([s.report(report_type='bed') for s in segments])
        return data

    def segment_report (self, segment):

        """ Generates a report for Segment objects """
        if self._report_type == 'bed':    
            report = '\t'.join([str(p) for p in [segment.chrom, 
                                                 segment.start,
                                                 segment.end,
                                                 segment.end - segment.start,
                                                 segment.parameters['ai'], 
                                                 segment.parameters['n'],
                                                 segment.parameters['m'],

                                                 2*segment.parameters['m']/segment.genome_medians['m0'],
                                                 segment.parameters['d_HE'], 
                                                 segment.parameters['score_HE'], 
                                                 segment.parameters['model'],
                                                 segment.parameters['d_model'],
                                                 segment.parameters['AB'],
                                                 segment.parameters['(AB)(2+n)'],
                                                 segment.parameters['(AB)(2-n)'],
                                                 segment.parameters['A'],
                                                 segment.parameters['AA'],
                                                 segment.parameters['AAB'],
                                                 segment.parameters['AAAB'],
                                                 segment.parameters['AAA'],
                                                 segment.parameters['AAAA'],
                                                 segment.parameters['A+AA'],
                                                 segment.parameters['AAB+AAAB'],
                                                 segment.parameters['AA+AAA'],
                                                 segment.parameters['AA+AAB'],
                                                 segment.parameters['AAB+AABB'],
                                                 segment.parameters['AAA+AAAA'],
                                                 segment.parameters['score_model'],
                                                 segment.parameters['k'], 
                                                 segment.cytobands,
                                                 segment.centromere_fraction]])
                                                 #segment.parameters['call'], segment.parameters['call_FDR']]])
        else:
            report = ''
        return report

    def run_report(self, run):
        """ Generates a report for Run objects """
        fields = ['chi2', 'chi2_noO', 'positions', 'p_norm', 'merged_segments']
        lines = ['Run: ' + run.name]
        for solution in run.solutions:
            lines.append (str("\t")+'Solution')
            soldic = solution._asdict()
            for f in fields:
                lines.append (str("\t\t") + f + ': '+ str(soldic[f]))
        return '\n'.join(lines)