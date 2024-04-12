##Dictonary with models is define now in one place

from collections import namedtuple
import numpy as np
import scipy.optimize as opt

def pen_scaling_sigmoid(std_cn, k=50, center=0.10428316240586664):
    return 1 - (0.5 / (1 + np.exp(-k * (std_cn - center))))

def find_min(d):
    min_key = None
    min_value = float('inf')

    for key, value in d.items():
        if isinstance(value, float) and not np.isnan(value):
            if value < min_value:
                min_value = value
                min_key = key
    return min_key

Preset = namedtuple ('Preset', ['k','m', 'ai'])

#Presets of models with 2 +/- 1 copies

model_presets = {
                 '(AB)(2+n)' : Preset(k = lambda m,dv,m0: np.abs(m/m0 - 1),
                                      m = lambda k,m0: (1+k)*m0,
                                      ai = lambda k,m0:  np.zeros_like(k)),

                 '(AB)(2-n)' : Preset(k = lambda m,dv,m0: np.abs(m/m0 - 1),
                                      m = lambda k,m0: (1-k)*m0,
                                      ai = lambda k,m0:  np.zeros_like(k)),
                   
                 'A'   : Preset(k = lambda m,dv,m0: 2*dv/(0.5+dv),
                                m = lambda k,m0: (2-k)*m0/2,
                                ai = lambda k,m0: k/(2*(2-k))),
                 
                 'AA'  : Preset(k = lambda m,dv,m0: 2*dv,
                                m = lambda k,m0: (np.zeros_like(k)+1)*m0,
                                ai = lambda k,m0: k/2),
                 
                 'AAB' : Preset(k = lambda m,dv,m0: 2*dv/(0.5-dv),
                                m = lambda k,m0: (2+k)*m0/2,
                                ai = lambda k,m0: k/(2*(2+k))),
    
                 'AAAB' : Preset (k = lambda m,dv,m0 : 2*dv/(1-2*dv),
                                  m = lambda k,m0 : (1+k)*m0,
                                  ai = lambda k,m0 : k/(2+2*k)),
                   
                 'AAA' : Preset (k = lambda m,dv,m0 : 4*dv/(3-2*dv),
                                 m = lambda k,m0 : (2+k)*m0/2,
                                 ai = lambda k,m0 : 3*k/(4+2*k)),
                    
                 'AAAA' : Preset (k = lambda m,dv,m0 : dv/(1-dv),
                                  m = lambda k,m0 : (1+k)*m0,
                                  ai = lambda k,m0 : k/(1+k)),
                 
                 'A+AA' : Preset (k = lambda m,dv,m0 : np.abs(2*m/m0 - 1), # TODO: new
                                  m = lambda k,m0: (1 + k)*m0/2,
                                  ai = lambda k,m0 : 0.5 * np.ones_like(k)),
                   
                 'AAB+AAAB' : Preset (k = lambda m,dv,m0 : (6*dv-1)/(1-2*dv),
                                      m = lambda k,m0 : (3+k)*m0/2,
                                      ai = lambda k,m0 : (1+k)/(6+2*k)),
                   
                 'AA+AAA' : Preset (k = lambda m,dv,m0 : np.abs(2*m/m0 - 2), #TODO: new
                                    m = lambda k,m0 : (1 + k/2)*m0,
                                    ai = lambda k,m0 : 0.5 * np.ones_like(k)),
                 
                 'AA+AAB' : Preset (k = lambda m,dv,m0 : 2*(1-2*dv)/(2*dv+1),
                                    m = lambda k,m0 : (2+k)*m0/2,
                                    ai = lambda k,m0 : (2-k)/(4+2*k)),
                   
                 'AAB+AABB' : Preset (k = lambda m,dv,m0 : (1-6*dv)/(2*dv+1),
                                      m = lambda k,m0 : (3+k)*m0/2,
                                      ai = lambda k,m0 : (1-k)/(6+2*k)),
                 
                 'AAA+AAAA' : Preset (k = lambda m,dv,m0 : np.abs(2*m/m0 - 3), #TODO: new
                                      m = lambda k,m0 : (3+k)*m0/2,
                                      ai = lambda k,m0 : 0.5 * np.ones_like(k))
                 }

def pick_model(ai, s_ai, cn, s_cn, models):
    assert not np.isnan (ai), "ai is nan"
    assert not np.isnan (cn), "cn in nan"
    
    l = pen_scaling_sigmoid(s_cn)
    
    model_weights = {
        '(AB)(2+n)': 2,
        '(AB)(2-n)': 2,
        'A': 2,
        'AA': 3,
        'AAB': 2,
        'AAAB': 3,
        'AAA': 4,
        'AAAA': 5,
        'A+AA': 4, 
        'AAB+AAAB': 4,
        'AA+AAA': 5,
        'AA+AAB': 4,
        'AAB+AABB': 4,
        'AAA+AAAA': 6
    }
    
    results = {'model' : None, 'd_model' : np.nan, 'k': np.nan, 
                'AB': np.nan, '(AB)(2+n)': np.nan, '(AB)(2-n)': np.nan, 'A': np.nan, 'AA': np.nan,
                'AAB': np.nan, 'AAAB': np.nan, 'AAA': np.nan, 'AAAA': np.nan, 'A+AA':np.nan,
                'AAB+AAAB': np.nan, 'AA+AAA':np.nan, 'AA+AAB': np.nan, 'AAB+AABB': np.nan, 'AAA+AAAA': np.nan}
    
    scaled_model_weights = model_weights.copy()
    for key in scaled_model_weights:
        if key in ['A', 'AA', 'A+AA']:
            scaled_model_weights[key] *= l
    
    for model in models:
        d, k = calculate_distance_minim(ai, s_ai, cn, s_cn, model_presets[model]) 
        nll = 0.5 * d ** 2
        d = nll + scaled_model_weights[model]
        results[model] = d
            
    da = ai/s_ai
    dc = (cn-2)/s_cn
    dd = 0.5 * (np.sqrt(da**2+dc**2)) ** 2
    
    results['AB'] = dd
    selected = find_min(results)
    results['model'] = selected
    results['d_model'] = results[selected]
    if selected == 'AB':
        results['k'] = 0
    else:
        _, k = calculate_distance_minim(ai, s_ai, cn, s_cn, model_presets[selected])
        results['k'] = k
        
    return results
    
def calculate_distance_minim (ai, s_ai, cn, s_cn, model):
    
    res = opt.minimize_scalar (dist, bounds = (0,1), args = ((ai, s_ai, cn, s_cn, model)), 
                               method = 'bounded')
    
    return res.fun, res.x

def dist (k, ai, s_ai, cn, s_cn, model):
    da = (ai - model.ai(k,2))/s_ai
    dc = (cn - model.m(k,2))/s_cn
    return np.sqrt(da**2+dc**2)