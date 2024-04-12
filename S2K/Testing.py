import scipy.optimize as opt
import numpy as np
import multiprocessing as mpl
import pandas as pd
from math import ceil,floor
import scipy.stats as sts
from collections import namedtuple

from S2K import Consts


COV_results = namedtuple ('COV_results', ['m', 'dm', 'l', 'dl'])
VAF_results = namedtuple ('VAF_results', ['chi2', 'vaf', 'fb'])
HE_results = namedtuple ('HE_results', ['chi2', 'vaf', 'cov', 'b'])

class Testing:
    """
    Class that performs tests.
    """
    
    def __init__  (self, test_name, chromosomes, logger):
        """
        Class constructor        
        """
        assert test_name in ['COV', 'HE', 'VAF', 'HEn'], f"Unknown test: {test_name}!"
        self.test_name = test_name
        i = np.where([t == test_name for t in ['COV', 'HE', 'VAF', 'HEn']])[0][0]
        self.test = [COV_test, HE_test, VAF_test][i]
        self.chromosomes = chromosomes
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.test.__name__}')
        self.logger.debug (f'Object {self.test.__name__} created.')
        
    def run_test (self, *args, no_processes = 1, exclude_symbol = [], **kwargs):
        """
        Method that runs test
        """        
        if no_processes > 1:
            self.logger.debug (f'Runnig test in {no_processes} processes.')
            with mpl.Pool (processes = no_processes) as pool:
                results = pool.starmap (self.test, [(self.chromosomes[chrom].data, args, kwargs) 
                                                    for chrom in self.chromosomes.keys()])
        else:
            results = []
            index = []
            for chrom in self.chromosomes.keys():
                res = self.test (self.chromosomes[chrom].data, args, exclude_symbol = exclude_symbol)
                if res is not None:
                    results.append (res)
                    index.append (chrom)
                self.logger.debug (f'Runnig test for {chrom}.')
        self.results = pd.DataFrame.from_records (results, columns = results[0]._fields, 
                                                  index = index)
        self.logger.info (f'Test finished.')
        
    def analyze (self, parameters = {'alpha' : 0.01, 'r' : 0.5}, q = (5,99), outliers = [], skip_par = []):
        """ 
        Method to analyze results of the test.
        """        
        self.logger.debug (f'Starting to analyze test results.')
        self.normal_range = {}
        columns = self.results.columns
        status = {}
        for column in columns[[c not in skip_par for c in columns]]:
            try:
                alpha = float(parameters[column + '_alpha'])
                r = float(parameters[column + '_r'])
            except: 
                alpha, r = (0.01, 0.5)
            
            self.logger.debug (f'Parameter {column} being analyzed with alpha = {alpha} and r = {r}')
            res = self.results.loc[[c not in outliers for c in self.results.index.tolist()] & (self.results.notna().all(axis = 1)), column].values
            if len(res) == 0:
                self.logger.critical (f'No parameters to work on')
                exit(1)
            if all(np.isnan(np.unique(res))):
                self.logger.critical (f'Parameter {column} undetermined!')
                exit(1) 
            elif len (np.unique(res)) < 5:
                param_range = (res.min(), res.max())
                self.logger.warning(f'Parameter {column} range estimation based on normal approximation not possible. Min max used.')
            else:
                try:
                    param_range = get_outliers_thrdist (np.unique(res), alpha, r)
                except:
                    self.logger.exception (f'Test {self.test_name}: estimation of {column} distribution failed. Maybe BMT or bad measurement?')
                    exit (1)
            
            self.logger.info ('Estimated normal ragne of {} is from {} to {}'.format (column, *param_range))
            in_or_out = (self.results[column].values >= param_range[0]) & (self.results[column].values <= param_range[1])
            self.results[column + '_status'] = in_or_out
            
            status[column] = ['inlier' if iou else 'outlier' for iou in in_or_out]
            self.logger.info (f"{sum(in_or_out)} chromosomes' {column} within range.")
            self.normal_range[column] = param_range
        
        self.status = pd.DataFrame.from_dict (status)
        self.status['chrom'] = self.results.index
        self.status.set_index ('chrom', inplace = True)
        inlier_filter = (self.status == 'inlier').all(axis = 1) 
        self.medians = self.results.loc[inlier_filter, columns].median (axis = 0, numeric_only = True)

    def get_parameters (self, chromosome):
        
        try:
            test_results = self.results.T[chromosome].T
        except KeyError:
            test_results = self.medians
            self.logger.debug (f"No parameters for chromosome {chromosome}")
        return test_results
    
    def get_status (self, chromosome):
        #try:
        #    status = (self.status.T[chromosome] == 'inlier').T
        #except KeyError:
        #    status = False
        if chromosome not in self.status.index.values.tolist():
            #self.status.T[chromosome] = 'outlier'
            self.status.loc[chromosome] = 'outlier' 
        return (self.status.T[chromosome] == 'inlier').T
    
    def get_genome_medians (self):
        return self.medians
    
    def get_inliers (self):
        return self.status[(self.status == 'inlier').all(axis = 1)].index.values.tolist()
    
    def get_outliers (self):
        return self.status[(self.status == 'outlier').all(axis = 1)].index.values.tolist()
    
    def report_results (self) -> pd.DataFrame:
        """Method to report only test results without status."""
        all_columns = self.results.columns.tolist()
        indexes = np.where([c.find('status') ==  -1 for c in all_columns])[0]
        return self.results[[all_columns[i] for i in indexes]]
    
def COV_test (data, *args, **kwargs):
    
    initial_shape = Consts.COV_INITIAL_SHAPE if 'initial_shape' not in kwargs else kwargs['initial_shape']
    shape_range = Consts.COV_SHAPE_RANGE if 'shape_range' not in kwargs else kwargs['shape_range']
    exclude_symbol = [] if 'exclude_symbol' not in kwargs else kwargs['exclude_symbol']    
    percentiles = np.arange (0.01, 0.99, 0.01)
     
    if len(exclude_symbol):
        covs = data.loc[[s not in exclude_symbol for s in data.symbol.tolist()] , 'cov'].values
    else:
        covs = data['cov'].values
    
    try:
        y = np.percentile (covs, percentiles*100)
        cov_range = (floor(y[0]*0.5) , ceil (y[-1]*1.05))
        initial_m = np.median(covs)
             
        #here may be good place for try except
        popt, pcov = opt.curve_fit (lambda_ppf, percentiles, y, p0 = [initial_m, initial_shape],
                                    bounds = [(cov_range[0],shape_range[0]),
                                              (cov_range[1],shape_range[1])])
        m, l  = popt
        dm, dl  = np.sqrt(np.diag(pcov)) 
        return COV_results (m = m, dm = dm, l = l, dl = dl)
    except:
        return None
    
def Q (p,l):
    if l == 0:
        return np.log(p/(1-p))
    else:
        return (p**l-(1-p)**l)/l    

def lambda_ppf (p, m, lam):
    return Q (p, lam)*np.sqrt(m)+m

def HE_test (data, *args, **kwargs):
    #cov_perc_bounds = Consts.HE_COV_PERC_BOUNDS if 'cov_perc_bounds' not in kwargs else kwargs['cov_perc_bounds']
    vaf_bounds = Consts.HE_VAF_BOUNDS if 'vaf_bounds' not in kwargs else kwargs['vaf_bounds']
    fcov_bounds = Consts.HE_FCOV_BOUNDS if 'fcov_bounds' not in kwargs else kwargs['fcov_bounds']
    fN_bounds = Consts.HE_FN_BOUNDS if 'fN_bounds' not in kwargs else kwargs['fN_bounds']
    a_bounds = Consts.HE_A_BOUNDS if 'a_bounds' not in kwargs else kwargs['a_bounds']
    aN_bounds = Consts.HE_AN_BOUNDS if 'aN_bounds' not in kwargs else kwargs['aN_bounds']
    b_bounds = Consts.HE_B_BOUNDS if 'b_bounds' not in kwargs else kwargs['b_bounds']
    lerr_bounds = Consts.HE_LERR_BOUNDS if 'lerr_bounds' not in kwargs else kwargs['lerr_bounds']   
    
    #this failed 
    #cov_min, cov_max = np.percentile (data['cov'].values, q = cov_perc_bounds)
    
    mid, mid_sd = np.percentile (data['cov'].values, q = (50,66))
    out_min, out_max = get_outliers_thrdist (np.sort(data['cov'].values), alpha = 2.8e-7, r = 0.2)
    
    if (int(mid+5*(mid_sd-mid)))<(out_max):
        cov_min = int(np.max((mid-5*(mid_sd-mid), np.percentile(data['cov'].values, q = (0.1)))))
        cov_max = int(mid+5*(mid_sd-mid))
    else:
        cov_min = int(np.max((out_min, np.percentile(data['cov'].values,q = (0.1)))))
        cov_max = int(out_max)

    n = np.concatenate ([np.repeat(c, c+1) for c in np.arange (cov_min, cov_max +1)])
    a = np.concatenate ([np.arange(0, c+1) for c in np.arange (cov_min, cov_max +1)])
    c = np.zeros_like (a)
    
    data_of_interest = data.loc[(data['cov'] >= cov_min)&(data['cov'] <= cov_max)]
    counts = data_of_interest.groupby (by = ['cov', 'alt_count'])['chrom'].count()
    #this should run one time, so much quicker
    for i in np.arange(len(n)):
        try:
            c[i] = counts[(n[i],a[i])]
        except KeyError:
            pass
     
    fcov = (data['cov'].median() - cov_min) /  (cov_max - cov_min)
    aH = len(data.loc[(data['vaf'] > 0.1)&(data['vaf'] < 0.9)])/len(data)
    aN = len(data.loc[(data['vaf'] < 0.1)])/len(data)
    N = len(data_of_interest['cov'])#.sum()
    
    res = opt.minimize (fun_chi2, x0 = (0.5, fcov, 1.0, aH, aN, 1.3, 6, 6), args = (n,a,c,N),
                        bounds = (vaf_bounds, fcov_bounds, fN_bounds, a_bounds, aN_bounds, b_bounds, lerr_bounds, lerr_bounds),
                        options = {'maxiter' : 2000})
    
    if res.success:
        chi2 = res.fun
        vaf, fcov, fN, aH, aN, b, le, lf = res.x    
        cov = cov_min + fcov*(cov_max-cov_min)    
    else:
        chi2 = np.nan
        vaf = np.nan
        cov = np.nan
        b = np.nan
    
    return HE_results(chi2 = chi2, vaf = vaf, cov = cov, b = b)
    #return chi2, vaf, fcov, fN, aH, aN, b, le, lf 
    #return chi2, vaf, cov, b

def fun_chi2 (params, n, a, c, N):
        vaf, fcov, fN, aH, aN, b, le, lf = params
        fe = 10**(-le)
        ff = 10**(-lf)
        
        cov = n.min() + fcov*(n.max()-n.min())
        ns = cn2_cov_pdf (n, cov, b)
        
        nhe = cn2_vaf_pdf (vaf, n)
        nho = HO_vaf_pdf (a, n, fe)
        nno = NO_vaf_pdf (a, n, ff)
        
        pt = ns*(aH*nhe + (1-aH -aN)*nho + aN*nno)
        ct = 2*fN*N*pt/pt.sum()
#       ct = ns*(aH*nhe + (1-aH -aN)*nho + aN*nno)
        
        chi2i = (c - ct)**2/np.sqrt(ct*ct+1)
        
        chi2 = chi2i.sum()/len(chi2i)
        
        return chi2

def cn2_vaf_pdf (v,c):
    p = []
    for d in np.arange (c[0], c[-1]+1):
        sv = np.sqrt((v*(1-v))/d)    
        x = np.arange (-0.5, d+1)/d
        pc = sts.norm.cdf (x[1:], v, sv) - sts.norm.cdf (x[:-1], v, sv)
        p.append (pc)
    return np.concatenate (p)

def cn2_cov_pdf (n,c,b = 1):
    p = sts.norm.pdf (n, c, np.sqrt(b*c))
    return p

def HO_vaf_pdf (i, n, fe = 10**-6, b = 1):
    p = sts.binom.pmf(i, n, 1-b*fe)
    return p
    
def NO_vaf_pdf (i, n, fe = 10**-6, b = 1):
    p = sts.binom.pmf(i, n, b*fe)
    return p

def VAF_test (data, m, **kwargs):
    
    vaf_bounds = Consts.VAF_VAF_BOUNDS if 'vaf_bounds' not in kwargs else kwargs['vaf_bounds'] 
    n_thr = Consts.VAF_N_THR if 'n_thr' not in kwargs else kwargs['n_thr']
    run_fb = True if 'run_fb' not in kwargs else kwargs['run_fb']    
    
    def chi2 (v, counts):
        chi2 = 0
        cc = 0
        for c, h, n in counts:
            if n >= n_thr:
                cc += c + 1 
                na = n*sts.binom.pmf (np.arange(c+1),c,v)
                chi2 += sum((h-na)**2/np.sqrt(na*na+1))
        
        if cc == 0: 
            cc = 1        
        return chi2/cc

    counts = []
    if hasattr(m, '__len__'):
        m = m[0]
    cov_min = m - 2*np.sqrt(m)
    cov_max = m + 2*np.sqrt(m)
    filt = (data['symbol'] == Consts.E_SYMBOL) & (data['cov'] > cov_min) & (data['cov'] < cov_max)
    for c,d in data.loc[filt].groupby (by = 'cov'):
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        n = len(d)
        counts.append ((c, h[0], n))
        del (h)
    ###safe if counts is empty
    c_max = 0
    for c in counts:
        if c[2] > c_max:
            c_max = c[2]        
    if c_max > n_thr:
        res = opt.minimize (chi2, x0 = (0.5) , args = (counts), method = 'L-BFGS-B',
                        bounds = (vaf_bounds,))        
        if run_fb:
            fb = find_fb (np.sort(data.loc[data['symbol'] == Consts.E_SYMBOL, 'vaf'].values), m,
                                  f_max = Consts.FB_F_MAX, eps = Consts.FB_EPS)
        else:
            fb = np.nan    
        results = VAF_results (chi2 = res.fun, vaf = res.x[0], fb = fb)
    else:
        results = VAF_results (chi2 = 0, vaf = 0, fb = np.nan)
    return results

def find_fb (vafs, m, f_max = 1.4, eps = 1e-4):
    
    def two_gauss (v, dv, a):
        v0 = 0.5
        return a*sts.norm.cdf (v, v0 + dv, smin)+(1-a)*sts.norm.cdf (v,v0 - dv,smin)
    
    s = np.sqrt (0.25/m)
    f = 1.0
    df = 0.01
    smin = s*f
    cdf = np.linspace (0,1, len(vafs))
    popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0, 0.0),
                                      (1, 1)))
    a_prev = popt[1]
    f = f + df
    popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0, 0.0),
                                      (1, 1)))
    a_next = popt[1]

    while (a_next > eps)&(f < f_max):
        #new f
        f += df
        smin = s*f
        popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = popt,
                            bounds = ((0, 0.0),
                                      (1, 1)))
        a_prev = a_next
        a_next = popt[1]
    
    a_zero = a_next
    a_nzero = a_prev
    
    df = 0.001
    
    while (a_nzero > 2*eps)&(a_nzero - a_zero > eps)&(np.abs(df) > eps)&(f < f_max):
        f = f + df
        smin = s*f
        popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0.0, 0.0),
                                      (1.0, 1.0)))
        if popt[1] > eps:
            a_nzero = popt[1]
            df = np.abs(df/2)
        else:
            df = -np.abs(df/2)
        
    if f < f_max:
         ff = f - df
    else:
        ff = np.nan
    
    return ff

def get_outliers_thrdist (values, alpha = 0.01, r = 0.5):
    
    n = len (values)
    current = values
    core = values
    mc = np.mean (current)
    sc = np.std (current)
    for i in range (int(np.floor(n*r))):
        m = np.mean (current)
        s = np.std (current)
        R = np.abs ((current -m)/s)
        j = np.where (R == R.max())[0][0]
        
        ni = n - i
        p = 1 - alpha/(ni-1)
        t = sts.t.ppf (p, ni-1)
        l = (ni*t)/np.sqrt((ni - 1 +t*t)*(ni+1))
        
        current = np.concatenate ([current[:j], current[j+1:]])
        if R.max() > l:
            core = current.copy()
            mc = np.mean(current)
            sc = np.std (current)
            
    return (sts.norm.ppf (alpha, mc,sc), sts.norm.ppf (1-alpha, mc,sc))
