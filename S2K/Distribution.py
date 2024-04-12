import numpy as np
import scipy.optimize as opt
import scipy.stats as sts

from S2K import Consts
from S2K.Testing import get_outliers_thrdist


class Distribution:
    """
    Class design to describe the data with either single or double Gauss distributions with managing outliers.

    Atributes:
        parameters (dic): Parameters of estimated distribution(s)
        all_parameters (dic): Parameters estimated for both situations, if done
        key (string): Resolved to {single/double} distribution
    """
    
    def __init__ (self, values, p_thr = 0.3, thr_z = 1.0):
        """
        Class constructor. Estimates distribution(s)' parameters based on provided values.

        """
        assert len (values) > Consts.LENGTH_THRESHOLD, print ("Not enough data points to consider distribution.")
        
        try: 
            single_G_par = fit_single_G (np.sort(values), alpha = 0.01, r = 0.5)
            values_used = 'all'
        except:
            single_G_par = fit_single_G (np.sort(np.unique(values)), alpha = 0.01, r = 0.5)
            values_used = 'unique'
        
                   
        self.all_parameters = {}
        self.all_parameters['single'] = single_G_par
        self.all_parameters['single']['values_used'] = values_used
        
        z = np.abs(values-single_G_par['m'])/single_G_par['s']
        string = list ('O' * len(values))
        for i in np.where ((z < thr_z))[0]:
            string[i] = 'B'
        
        self.all_parameters['single']['string'] = string
        
        ##UGLY solution                
        if single_G_par['p'] < p_thr:
            try:
                double_G_par = fit_double_G (np.sort(values), alpha = 0.01, r = 0.5)
                values_used = 'all'
            except:
                double_G_par = fit_double_G (np.unique(np.sort(values)), alpha = 0.01, r = 0.5)
                values_used = 'unique'
                
            m1 = double_G_par['m'][0]
            m2 = double_G_par['m'][1]
            if (double_G_par['p'] > single_G_par['p'])&(np.abs(m1 - m2)/(m1+m2) > 0.01):

                self.key = 'double'
                self.parameters = double_G_par
                self.all_parameters ['double'] = double_G_par
                self.all_parameters ['double']['values_used'] = values_used
            
                z0 = np.abs(values-double_G_par['m'][0])/double_G_par['s'][0]
                z1 = np.abs(values-double_G_par['m'][1])/double_G_par['s'][1]
                #generate string
                string = list ('O' * len(values))
                #marking C population
                for i in np.where ((z0 < z1)&(z0 < thr_z))[0]:
                    string[i] = 'C'
                #marking D population
                for i in np.where ((z1 < z0)&(z1 < thr_z))[0]:
                    string[i] = 'D'
                self.parameters['string'] = ''.join(string)
                self.all_parameters['double']['string'] = string
                self.string = ''.join(string)
            else:
                self.key = 'single'
                self.parameters = single_G_par
                self.string = ''.join(string)
        else:
            self.key = 'single'
            self.parameters = single_G_par
            self.string = ''.join(string)
        
        
    def fail_normal (self):
        """
        Method to check if single Gauss failed.
        """
        return self.key == 'double'
        
    def combinations_of_params (self, size = 1, key = 'single', reverse = False):
        """Method to generate parameters in desired shape"""
        try:
            parameters_sets = self.all_parameters[key]
        except KeyError:
            raise (f'The key: {key} not in solutions. Go away!')
        
        try:
            m = parameters_sets[str(size)]['m']
            s = parameters_sets[str(size)]['s']
        except KeyError:
            raise (f'The size: {size} and key: {key} not available. Go away!')
        
        if reverse:
            return m[::-1], s[::-1]
        else:
            return m, s
        
            
    def to_string (self):
        return (self.parameters)
        
    def __str__ (self):
        return (self.parameters)
    
    def __repr__ (self):
        return (self.parameters)
    
    
def fit_single_G (values, alpha = 0.01, r = 0.5):
    """
    Function to fit Gauss to _values_
    
    Returns dictionary with parameters.
    """
    try:
        thr = get_outliers_thrdist (np.sort(values), alpha, r)
    except: 
        thr = (min(values), max(values))

    a = np.sort(values[(values >= thr[0])&(values <= thr[1])])
    popt, pcov = opt.curve_fit (sts.norm.cdf, a, np.linspace(0,1, len(a)), 
                                p0 = [np.mean(a), np.std (a)])
    ksp = sts.kstest (a, sts.norm.cdf, args = popt)
    return {'p' : ksp.pvalue, 
            'm': popt[0],
            's': popt[1],
            '1' : {'m': np.array(popt[0]),
                   's': np.array(popt[1])},
            '2' : {'m': np.array([popt[0], popt[0]]),
                   's': np.array([popt[1], popt[1]])},
            'thr' : thr, 'a' : np.ones(1)}

def fit_double_G (values_all, alpha, r = 0.5, initial_bounds = None, initial_p0 = None):
    """
    Function to fit two Gauss' to _values_
    
    Returns dictionary with parameters.
    """
    
    try:
        thr0 = get_outliers_thrdist (np.sort(values_all), alpha, r)   
    except:
        thr0 = (min(values_all),max(values_all))
 
    values = values_all[(values_all >= thr0[0]) & (values_all <= thr0[1])]
    
    if initial_p0 is None:
        p0 = (0.5, np.percentile (values, 25), np.percentile(values,40)-np.percentile(values,10),
              np.percentile (values, 75), np.percentile(values,90)-np.percentile(values,60))
    else:
        p0 = initial_p0
        
    
    #ValueError: Each lower bound must be strictly less than each upper bound

    if initial_bounds is None:
        bounds = [[0, np.percentile (values, 5), 0.1*np.percentile (values, 40)-0.1*np.percentile (values, 10),
                  np.percentile (values, 55), 0.1*np.percentile (values, 90)-0.1*np.percentile (values, 60)],
                 [1, np.percentile (values, 45), 2*np.percentile (values, 40)-2*np.percentile (values, 10),
                  np.percentile (values, 95), 2*np.percentile (values, 90)-2*np.percentile (values, 60)]]
    else:
        bounds = initial_bounds

    popti, pcovi = opt.curve_fit (gaus2, np.sort(values), np.linspace (0,1,len(values)), p0 = p0,
                                  bounds = check_bounds(bounds)) 

    if popti[3] > popti[1]:
        out_max = sts.norm.ppf (1-alpha, popti[3], popti[4])
        out_min = sts.norm.ppf (alpha, popti[1], popti[2])
    else:
        out_max = sts.norm.ppf (1-alpha, popti[1], popti[2])
        out_min = sts.norm.ppf (alpha, popti[3], popti[4])
        tmp = popti
        popti = (tmp[0], tmp[3], tmp[4], tmp[1], tmp[2])


    a = values[(values <= out_max)&(values >= out_min)]
    p0 = (0.5, np.percentile (a, 25), np.percentile(a,40)-np.percentile(a,10),
          np.percentile (a, 75), np.percentile(a,90)-np.percentile(a,60))
    
    bounds = [[0, np.percentile (a, 5), 0.1*np.percentile (a, 40)-0.1*np.percentile (a, 10),
                  np.percentile (a, 55), 0.1*np.percentile (a, 90)-0.1*np.percentile (a, 60)],
              [1, np.percentile (a, 45), 2*np.percentile (a, 40)-2*np.percentile (a, 10),
                  np.percentile (a, 95), 2*np.percentile (a, 90)-2*np.percentile (a, 60)]]
    popt, pcov = opt.curve_fit (gaus2, np.sort(a), np.linspace (0,1,len(a)), p0 = p0,
                                  bounds = check_bounds(bounds)) 
    
    a0,m0,s0,m1,s1 = popt
    da0,dm0,ds0,dm1,ds1 = np.sqrt(np.diag(pcov))
    ksp = sts.kstest (a, gaus2, args = popt)

    return {'p' : ksp.pvalue,
            'm' : np.array([m0,m1]),
            's' : np.array([s0,s1]),
            'a' : np.array([a0, 1-a0]),
            '2':{'m': np.array([m0,m1]), 
                 's': np.array([s0,s1]), 'a' : np.array([a0, 1-a0])},
            'thr' : np.array((out_min, out_max))}

def  check_bounds(bounds):
    lower = bounds[0].copy()
    upper = bounds[1].copy()
    check = [l>=u for l,u in zip(lower, upper)]
    if any (check):
        for i in np.where (check)[0]:
            lower[i] = lower[i]*0.5  
            upper[i] = upper[i]*1.5
        if upper[2] == 0.0:
            upper[2] = 0.05
        if upper[4] == 0.0:
            upper[4] = 0.05
    return [lower, upper]

def gaus2 (v, a, m0, s0, m1, s1):
    """
    Helper function returning cdf of two gauss'
    """
    return a*sts.norm.cdf (v, m0,s0) + (1-a)*sts.norm.cdf(v, m1, s1)

