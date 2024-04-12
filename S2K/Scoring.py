import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
import sklearn.linear_model as slm

from S2K import Consts
from S2K import Models


class Scoring:

    def __init__(self, fb = None, m0 = None, window_size = None, initial_data = np.array([]), logger = False) -> None:
        
        if logger:
                self.logger = logger.getChild (f'{self.__class__.__name__}')
        else:
            class lg:
                def __init__ (self):
                    self.info = lambda x: None
                
            self.logger = lg()
                    
        try:
            self.median_size = np.median(initial_data[:, 2])
            scales = np.sqrt(initial_data[:,2])            
            self.ai_param = fit_QQgauss(initial_data[: ,0]*scales, fit_intercept = False)
            self.cn_param = fit_QQgauss(initial_data[: ,1]*scales, fit_intercept = False)
            
            self.theor_ai_std = fb*0.25/np.sqrt(m0*window_size)
            self.theor_cn_std = 2/np.sqrt (m0*window_size)
            
            greater = self.ai_param['s'] > self.theor_ai_std  
            self.logger.info (f"Estimated std of ai is {'greater' if greater else 'smaller'} than theoretical minimum {self.theor_ai_std}.")
            if not greater:
                self.logger.info ("Std of ai is replaced by theoretical minimum.")
                self.ai_param['s'] = self.theor_ai_std
                
            greater = self.cn_param['s'] > self.theor_cn_std  
            self.logger.info (f"Estimated std of cn is {'greater' if greater else 'smaller'} than theoretical minimum {self.theor_cn_std}.")
            if not greater:
                self.logger.info ("Std of cn is replaced by theoretical minimum.")
                self.cn_param['s'] = self.theor_cn_std
                
            wayoff = self.cn_param['s'] > 16 * self.ai_param['s']
            if wayoff:
                self.logger.info ("Estimated std of ai is way too small than estimated std of cn.")
                self.logger.info ("Std of ai is replaced by (std of cn) / 16.")
                self.ai_param['s'] = self.cn_param['s']/16
                            
            dds = scales[:, np.newaxis]*(initial_data[:,:-1]) - np.array([self.ai_param['m'], self.cn_param['m']])

            ds = dds/np.array([self.ai_param['s'],self.cn_param['s']])[np.newaxis, :]
            self.dipl_dist = fit_QQgauss (np.sqrt((ds**2).sum(axis = 1)))
            self.dipl_dist['alpha'] = Consts.SCORE_ALPHA
            self.dipl_dist['thr'] = sts.norm.ppf (1-self.dipl_dist['alpha'], 
                                                  self.dipl_dist['m'],
                                                  self.dipl_dist['s'])

            self.logger.info (f"Median segment size: {self.median_size}")
            self.logger.info (f"Distribution of distance to diploid (0,2): m = {self.dipl_dist['m']}, s = {self.dipl_dist['s']}")
            self.logger.info (f"Distribution of diploid allelic imbalance: m = {self.ai_param['m']}, s = {self.ai_param['s']}")
            self.logger.info (f"Distribution of diploid copy number: m = {self.cn_param['m']}, s = {self.cn_param['s']}")
        
        except (IndexError):
            self.ai_param = {'m' : 0, 's' : 0}
            self.cn_param = {'m' : 0, 's' : 0}
            self.dipl_dist = {'m' : np.nan, 's' : 0, 'thr' : 0, 'alpha': np.nan}
            self.median_size = 1
            self.logger.info ("Scorer created with no diploid reference.")
        
            
    def get_d_thr (self):
        return self.dipl_dist['thr']

    def score_dipl (self, segment): 

        scale = np.sqrt (segment.parameters['n'])
        ai = segment.parameters['ai']*scale
        m = segment.parameters['m']
        m0 = segment.genome_medians['m0']
        cn = scale*(2*m/m0-2)

        s_ai = self.ai_param['s']
        s_cn = self.cn_param['s']
                
        if (s_ai > 0)&(s_cn > 0):
            # d = np.sqrt ((ai/s_ai)**2 + ((cn/s_cn)**2))
            d = 0.5 * (np.sqrt ((ai/s_ai)**2 + ((cn/s_cn)**2))) ** 2 # New distance to diploid updated to match the distance to models
            p_d = sts.norm.sf (d, self.dipl_dist['m'], self.dipl_dist['s'])
        
            segment.parameters['d_HE'] = d

            segment.parameters['p_HE'] = p_d
            segment.parameters['score_HE'] = -np.log10(p_d)
        else:
            segment.parameters['d_HE'] = np.inf
            segment.parameters['p_HE'] = 0
            segment.parameters['score_HE'] = np.inf

    def analyze_segment (self, segment, models):
        m = segment.parameters['m']
        m0 = segment.genome_medians['m0']
        cn = 2*m/m0
        
        s_ai = self.ai_param['s']
        s_cn = self.cn_param['s']
        
        ai = segment.parameters['ai']
            
        try:
            if segment.chrom == 'chrY':
                segment.parameters.update ({'model' : 'A', 'd_model' : 0.0,
                                            'k': cn, 'p_model' : np.nan, 'AB':np.nan, '(AB)(2+n)': np.nan,
                                            '(AB)(2-n)': np.nan, 'A': np.nan, 'AA': np.nan, 'AAB': np.nan, 'AAAB': np.nan, 'AAA': np.nan, 'AAAA': np.nan,
                                            'A+AA': np.nan, 'AAB+AAAB': np.nan, 'AA+AAA': np.nan, 'AA+AAB': np.nan, 'AAB+AABB': np.nan, 'AAA+AAAA': np.nan})
            else:
                segment.parameters.update (Models.pick_model(ai,s_ai,cn,s_cn,models))
        except (IndexError, AssertionError):
            segment.parameters.update ({'model' : 'UN', 'd_model' : np.nan,
                                        'k': np.nan, 'p_model' : np.nan, 'AB':np.nan, '(AB)(2+n)': np.nan,
                                        '(AB)(2-n)': np.nan, 'A': np.nan, 'AA': np.nan, 'AAB': np.nan, 'AAAB': np.nan, 'AAA': np.nan, 'AAAA': np.nan,
                                        'A+AA': np.nan, 'AAB+AAAB': np.nan, 'AA+AAA': np.nan, 'AA+AAB': np.nan, 'AAB+AABB': np.nan, 'AAA+AAAA': np.nan})
                
    def set_thr (self, thr):
        self.dipl_dist['thr'] = thr

def fit_QQgauss (values, fit_intercept = True):
    x = sts.norm.ppf (np.linspace (0,1,len(values)+2)[1:-1])
    huber = slm.HuberRegressor (fit_intercept = fit_intercept)
    huber.fit (x[:, np.newaxis], np.sort(values))
    return {'m' : huber.intercept_, 's' : huber.coef_[0]}
    
def fit_smallest_gauss (values):
    current_values = np.sort(values)
    thr = current_values.max()+1
    previous_thr = thr + 1    
    
    while (previous_thr > thr):
        current_values = current_values[current_values < thr]
        
        popt, pcov = opt.curve_fit (sts.norm.cdf, current_values, np.linspace (0,1,len(current_values)),
                                    p0 = [np.mean (current_values), np.std (current_values)])
        previous_thr = thr
        thr = sts.norm.ppf (1-1/(5*len(current_values)), *popt)
        
    return {'m': popt[0], 's' : popt[1], 'thr' : thr, 'alpha' : 1/(5*len(current_values))}