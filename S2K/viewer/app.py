from shiny import *
from .Plots import * #need a dot before Plots, like this: .Plots

from S2K import Models
from S2K import Consts
from S2K import Scoring


from S2K import Consts


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig


from collections import defaultdict
from matplotlib.patches import Ellipse
from sklearn.neighbors import KernelDensity

chromdic = {}
for c in Consts.CHROM_ORDER:
    chromdic[c] = c

def merge_records (all_records, chrom):
    
    filt_indexes = np.where([r.filt for r in all_records])[0]
    #records = [all_records[i] for i in filt_indexes]
    if len(filt_indexes) == 0:
        r = all_records[0]
        records = all_records
        status = 'NA'
        size = np.array ([r.end-r.start for r in records])
        m = np.array([r.m for r in records])
        cn = np.array([r.cn for r in records])
        k = np.abs(np.array([r.k for r in records]))
        #chi2 = (np.array([r.k_score for r in records])*np.log2(10))
        cyto = '-'.join([records[0].cyto.split('-')[0], records[-1].cyto.split('-')[-1]])
        model = 'NA'
        mr = (chrom, r.start, records[-1].end, 
              sum(m*size)/sum(size), sum(cn*size)/sum(size), model,
              sum(k*size)/sum(size), cyto, 'NA')
    else:
        
        records = [all_records[i] for i in filt_indexes]
        i = 0
        try:
            while not records[i].filt:
                i +=1 
        except:
            i = 0

        r = records[i]
        mr = ()
        if len (records) == 1:
            mr = (chrom, r.start, r.end, r.m, r.cn,
                r.model if r.status != 'norm' else 'AB',
                np.abs(r.k) if r.status != 'norm' else 0, r.cyto, r.k_score)
        else:
            status = r.status
            size = np.array ([r.end-r.start for r in records])
            m = np.array([r.m for r in records])
            cn = np.array([r.cn for r in records])
            k = np.abs(np.array([r.k for r in records]))
            chi2 = (np.array([r.k_score for r in records])*np.log2(10))
            cyto = '-'.join([records[0].cyto.split('-')[0], records[-1].cyto.split('-')[-1]])
            if status == 'norm':
                model = 'AB'
            else:
                model = r.model
            
            mr =  (chrom, records[0].start, records[-1].end,
                    sum(m*size)/sum(size), sum(cn*size)/sum(size), model,
                    sum(k*size)/sum(size), cyto, -np.log10(sts.chi2.sf(2*sum(chi2), 2*len(records))))
    
    return mr

  
app_ui = ui.page_fluid(
    ui.h2 ({"style" : "text-align: center;"}, "S2K results viewer."),
   
    ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Segments filtering:"),
                                       #ui.input_slider ('cent_thr', "Centromere fraction threshold",
                                       #                 value = 0.5, min = 0, max = 1),
                                       ui.input_slider ('k_thr', "Min clonality",
                                                        value = 0.05, min = 0, max = 1, step = 0.01),
                                       ui.input_slider ('size_thr', "Min non-centromere segment size (in Mb)",
                                                        value = 5, min = 0, max = 10),
                                       ui.h4 ("Display settings:"),
                                       ui.output_text ("auto_model_thr"),
                                       ui.input_slider ('model_thr', "Model score threshold",
                                                        value = 3, min = 0, max = 10, step = 0.1),

                                       ui.output_text ("auto_HE_thr"),
                                       ui.input_slider ('HE_thr', "Max HE score:",

                                                        value = 2, min = 0, max = 10),
                                       
                                       width = 2),
                      ui.panel_main(
                                    ui.navset_tab (
                                                   ui.nav("Genome-wide view",
                                                          ui.row(ui.input_file ('bed_file', "Choose BED file to upload:",
                                                                                 multiple = False, accept = '.bed'),),
                                                          ui.row(ui.column(8,
                                                                           ui.output_plot ('genome_plot'),),
                                                                 ui.column(4,
                                                                           ui.h6 ('Attempt of grouping CNVs into clones'),
                                                                           ui.output_plot ('clones_plot')),),
                                                          ui.row(ui.column(12,
                                                                           ui.h5 (''),
                                                                           ui.input_checkbox_group ('chroms_selected',
                                                                                                    "Select chromosomes to highlight",
                                                                                                    chromdic, inline = True)),),   
                                                          ui.row (ui.input_file ('par_file', "Choose PAR file to upload:",
                                                                                          multiple = False, accept = '.par')),
                                                          
                                                          ui.row(ui.h4('Diploid analysis'),
                                                                 ui.column(6,
                                                                           ui.h5 ('Diploid segments:'),
                                                                           ui.output_plot ('diploid_segments_plot'),
                                                                           ),
                                                                 ui.column(6,
                                                                           ui.h5 ('Diploid distance distribution:'),
                                                                           ui.output_plot ('diploid_distance_plot'),
                                                                           )),
                                                          
                                                          ui.row(ui.h4 ('Karyotyping review:'),
                                                                 ui.column(6, 
                                                                           ui.h5 ('Models review:'),
                                                                           ui.output_plot ('models_plot'),
                                                                           ),
                                                                 ui.column(6,
                                                                           ui.h5 ('Models distance distribution:'),
                                                                           ui.output_plot ('model_distance_plot'),
                                                                           ))
                                                         ),
                                                   
                                                   ui.nav("LOG", 
                                                          ui.row(ui.column(12,
                                                                 ui.input_file ('log_file', "Choose LOG file to screen:",
                                                                                multiple = False, accept = '.log'),
                                                                 ui.output_ui("dyn_log_ui")))),
 
                                                   ui.nav("CNVs",
                                                         ui.row(ui.column (2,
                                                                           ui.input_radio_buttons ('sort_CNV_by',
                                                                                                   "Sort list by:",
                                                                                                   {'position':'position', 'score':'score'}, inline = True))),
                                                         ui.row(ui.output_text ("number_CNVs"), ""),
                                                         ui.row(ui.output_table (id = 'CNVs'),)),
                                                   
                                                   ui.nav("Solution test",
                                                          ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Optimize settings:"),
                                                                                             ui.input_slider ('min_cn', "Min cov (relative):",
                                                                                                              value = 0.7, min = 0.1, max = 1, step = 0.05),
                                                                                             ui.input_slider ('max_cn', "Max cov (relative):",
                                                                                                              value = 1.1, min = 0.5, max = 1.5, step = 0.05),
                                                                                             ui.row(ui.column(6,ui.input_numeric ("step", 'Step:', 0.5, min = 0.5, max = 5, step = 0.5)),
                                                                                                    ui.column(6,#ui.h6 ("# steps"),
                                                                                                                ui.input_text ("number_points", '# points')),
                                                                                                   ),    
                                                                                             ui.h6("Coverage range:"),
                                                                                             ui.output_text_verbatim ("coverage_range"),
                                                                                             ui.input_action_button ('opt',
                                                                                                                     "Optimize solution"),
                                                                                             ui.input_checkbox_group ('models_selected',
                                                                                                    "Select models to use for refitting",
                                                                                                    choices = list(Models.model_presets.keys()),
                                                                                                    inline = False),
                                                                                             width = 2),
                                                                            ui.panel_main(ui.h6("Solution check plots:"),
                                                                                          ui.row(ui.column(4,
                                                                                                           ui.output_plot('solution_plot_dipl_opt', "Solution plot - diploid")                                                   ),
                                                                                                 ui.column(8, 
                                                                                                           ui.output_plot('solution_plot_opt', "Solution plot"))),
                                                                                          ui.h6("Total weighted distance to solution at coverage:"),
                                                                                          ui.output_plot('opt_plot', "Relative distance plot"),
                                                                                          ui.h6("Pick coverage to plot models:"),
                                                                                          ui.row(ui.input_slider ('m0_cov', "Diploid coverage",
                                                                                                               value = 0.5, min = 0, max = 1, width = '200%')),
                                                                                                 ui.output_text_verbatim ('solutions_info')))),
                                                                                          
                                                   ui.nav("Chromosome view",
                                                          ui.row(ui.column(12,
                                                                           ui.input_file ('data_file', "Choose data file to upload:",
                                                                                           multiple = False, accept = '.gz'),
                                                                           ui.input_radio_buttons ('chrom_view', "Choose chromosome to inspect",                                                                                    
                                                                                                   chromdic, inline = True),)),
                                                          ui.row(ui.input_radio_buttons('f_to_plot', "Which function plot to compare:", 
                                                                                         {'PDF' : 'Probability density function',
                                                                                          'CDF' : 'Cumulative density function'}, 
                                                                                         selected = 'CDF', inline = True)),                                                          
                                                          ui.row(ui.output_plot ('data_plot')),
                                                          ui.row(ui.column(6, 
                                                                           ui.output_plot ('compare_plot')),
                                                                 ui.column(6, 
                                                                           ui.output_plot ('cov_plot'))),
                                                          ui.row(ui.output_table (id = 'chrom_segments'))),

                                                    ui.nav("Report",
                                                           ui.row(ui.column(12,
                                                                            ui.output_plot('report_plot'),
                                                                            ui.input_checkbox_group ('report_models',
                                                                                                     "Pick balanced model types to include in the report:",
                                                                                                      {'(AB)(2-n)' : '(AB)(2-n)', '(AB)(2+n)' : '(AB)(2+n)'},
                                                                                                      inline = True),
                                                                            ui.output_table('report')))),

                                                    ui.nav("Publication",
                                                           ui.h6 ("AI still in school"))
                                                            
                                                   )
                                   )
          )
        ) 

def server(input, output, session):
    
    bed_full = reactive.Value(pd.DataFrame())
    bed = reactive.Value(pd.DataFrame())
    opt_bed = reactive.Value(pd.DataFrame())
    data = reactive.Value(pd.DataFrame())
    par = reactive.Value({})
    opt_solution = reactive.Value((np.array([]), np.array([]), np.array([])))
    chrom_sizes = reactive.Value(pd.Series())
    m0 = reactive.Value(np.nan)
    m0_opt = reactive.Value(np.nan)
    log_file = reactive.Value ([])
    bed_report = reactive.Value(pd.DataFrame())
    model_presets = reactive.Value (Models.model_presets)
    solutions_list = reactive.Value ({})    
          
    @output
    @render.text
    def solutions_info ():
        solutions = solutions_list()
        if len (solutions):
            header = "Paired model fit and diploid dist minima.\n"
            ms = np.array(list(solutions_list())) 

            d_total = np.array([solutions[m][1] for m in solutions.keys()])
            d_HE = np.array([solutions[m][0].dipl_dist['m'] for m in solutions.keys()])
                
            ind_d = sig.argrelmin (d_total, order = 2)[0]
            ind_HE = sig.argrelmin (d_HE, order = 1)[0]
            print (ind_d)
            print (ind_HE) 
                        
            dist = np.abs (ind_d[:,np.newaxis] - ind_HE[np.newaxis,:])
            index = np.argsort(dist, axis = 1)
            
            minims = []
            for i, id in enumerate(ind_d):
                idc = index[i,0]
                minims.append ((ms[id], ms[ind_HE[idc]], d_total[id], d_HE[ind_HE[idc]]))
             
            minims.sort (key = lambda x: x[2], reverse = False)
        
            return header + '\n'.join(['m_d = ' + str(m_d) + ' m_HE = ' + str(m_H) +'  d = ' + str(d) + ' d_HE = ' + str(HE)  for m_d, m_H,d,HE in minims])
  
   
    @output
    @render.ui
    def dyn_log_ui():
        return ui.TagList (ui.tags.textarea ([l.strip() for l in log_file()], 
                                             cols = "150", rows = "30"))
    
    @output
    @render.text
    def coverage_range ():
        if np.isnan(m0()):
            text = 'n.a.'
        else:
            text = '{:.2f}'.format(int(m0()*input.min_cn())) + ' - ' + '{:.2f}'.format(int(m0()*input.max_cn()))
        return text
    
    @reactive.Effect
    def _():
        try:
            text = str(len(np.arange(m0()*input.min_cn(), m0()*input.max_cn(), input.step())))
            ui.update_text ('number_points', value = text)
        except:
            pass
    
    
    @reactive.Effect
    @reactive.event (input.log_file)
    def _():
        file_input = input.log_file()
        if not file_input:
            return
        with open(file_input[0]['datapath']) as f:
            log = f.readlines()
        log_file.set(log)    
        
    @reactive.Effect
    @reactive.event(input.bed_file)
    def _():
        file_input = input.bed_file()
        
        if not file_input:
            return
        
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t', header = None, 
                    names = ['chrom', 'start', 'end', 'size', 'ai','n', 'm', 'cn',
                            'd_HE', 'score_HE', 'model', 'd_model', 
                            'AB', '(AB)(2+n)','(AB)(2-n)',
                            'A','AA','AAB','AAAB','AAA','AAAA',
                            'A+AA', 'AAB+AAAB', 'AA+AAA', 'AA+AAB','AAB+AABB', 'AAA+AAAA',                                                                                  
                            'score_model', 'k',
                            'cyto', 'cent'])
        
        df['size'] = (df['end'] - df['start'])/1e6
        
        data.set(pd.DataFrame())
        par.set({})
        bed_full.set(df)
        opt_solution.set ((np.array([]), np.array([]), np.array([])))
        log_file.set([])
        bed_report.set(pd.DataFrame())
        solutions_list.set ({})
            
    @reactive.Effect
    @reactive.Calc
    def _():
        tmp = bed_full()
        if len(tmp) > 0:
            chrom_sizes.set(tmp.groupby(by = 'chrom').agg({'end' : 'max'})['end'])
            tmp['filt'] = (False | ((1-tmp['cent'])*tmp['size'] >= input.size_thr())) & ((tmp['k'] >= input.k_thr()) | (tmp['model'] == 'AB'))
            #(tmp['cent'] <= input.cent_thr())
            bed.set (tmp.loc[tmp.filt])
            opt_bed.set(tmp.loc[tmp.filt])
    


    @reactive.Effect
    @reactive.Calc
    def _():
        #bf = bed_full()    
        b = opt_bed()
        if len(b):
            chrs = chrom_sizes().index.values.tolist()
            chrs.sort (key = Consts.CHROM_ORDER.index)
            
            report = b.loc[b['score_HE'] > input.HE_thr()].copy()
            report['model'] = [ m if score <= input.model_thr() else m+'*' for m,score in zip (report['model'].tolist(),
                                                                                              report['score_model'].tolist())]
            #report.loc[report['score_model'] > input.model_thr(), 'model'] = '--'
            #report.loc[report['score_model'] > input.model_thr(), 'k'] = np.nan
            
            #merging TBD
            
            bed_report.set(report)


    
    @reactive.Effect
    @reactive.event(input.par_file)
    def _():
        file_input = input.par_file()
       
        if not file_input:
            return
        pard = {}
        with open(file_input[0]['datapath'],'r') as f:
            for line in f.readlines():
                try:
                    key, value = line.split('\t')
                    pard[key] = float(value)
                except:
                    try:
                        key, values = line.split('\t')
                        pard[key] = [v.strip().strip("'[]") for v in values.split(',')]
                    except:
                        print ('Line: ' + line + 'not parsed.')
        par.set(pard)
        opt_solution.set ((np.array([]), np.array([]), np.array([])))
        ui.update_checkbox_group ('models_selected', selected = pard['models'])

        m0.set(pard['m0'])
        m0_opt.set(pard['m0'])
    
            
    @reactive.Effect
    @reactive.event(input.data_file)
    def _():
        file_input = input.data_file()
        
        if not file_input:
            return
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t')    
        data.set(df)
    
    @output
    @render.plot (alt = "Genomewide view")
    def genome_plot ():
        bed_data = bed()
        if  len(bed_data):
            fig, axs = plt.subplots (2, 1, figsize = (16,6), sharex = True)
            meerkat_plot (bed_data, axs, chrom_sizes(), cn_max=bed_data['cn'].max(),
                          model_thr = input.model_thr(), HE_thr = input.HE_thr())

            models = bed_data['model'].unique()
            for model in models:
                axs[0].plot ((),(), lw = 10, color = colorsCN[model], label = model)  
            
            low_cns = bed_data[(bed_data['cn']<0.6)]
            if len(low_cns)>0:
                axs[0].plot ((),(), lw = 10, color = 'purple', label = 'CN(0,1)')
            
            complexes = bed_data[(bed_data['cn']>4.5) & (bed_data['cn']<=10)]
            if len(complexes)>0:     
                axs[0].plot ((),(), lw = 10, color = 'yellow', label = 'CN[5,10]')
                
            high_cns = bed_data[(bed_data['cn']>=10) ]
            if len(high_cns)>0:
                axs[0].plot ((),(), marker='+', markersize=8, color='red', ls='', label='CN>10' )
                for i, row in high_cns.iterrows():
                    prev_chr = chr_order[chr_order.index(row['chrom'])-1]
                    start = chrom_sizes.loc[:prev_chr].sum()
                    annotate_loc = (start + row['start'], 5)
                    axs[1].plot(annotate_loc[0], annotate_loc[1], marker='+', markersize=8, ls='', color='red', label='')
        
            axs[0].legend (bbox_to_anchor = (0., 1.02, 1., 1.02),
                            ncol = int(len(model_presets())/2) + 1,
                            loc = 'lower left', mode = 'expand', fontsize = 8,
                            title = 'Models of mixed clones: normal (AB) and abnormal karyotypes:')
           
            return fig
    
    @output
    @render.plot (alt = 'Clones view')
    def clones_plot():
        bed_data = bed()
        if  len(bed_data):
            fig, ax = plt.subplots (1, 1, figsize = (6,6), sharex = True)
            tmp = bed_data.loc[bed_data['score_HE'] > input.HE_thr()]
            
            if len(tmp):
                k = np.sort(tmp['k'].values)[:, np.newaxis]
                order = np.argsort (tmp['k'].values)
                s = tmp['size'].values[order]
                c = [colorsCN[m] for m in tmp['model'].values[order]]
                        
                #how bandwidth?
                #median distance to the model? or similar
                kde = KernelDensity (kernel = 'gaussian', bandwidth = 0.018).fit (k)
                xt = np.linspace (min(-0.01, k[0,0]), k[-1,0]*1.1, 1000)[:, np.newaxis]
                ax.scatter (k, np.linspace (0,1, len(k)), s = np.sqrt(s)*10, c = c)
                ax.plot (k, np.linspace (0,1, len(k)), 'k:')
                axt = ax.twinx()
                axt.plot (xt, np.exp(kde.score_samples(xt)))
            
                ax.set_xlabel ('clonality')
                ax.set_ylabel ('clonaticty cdf')
                axt.set_ylabel ('clonaticty kde')
            
                return fig
    
    
        
    @output
    @render.plot (alt = "Genomewide view")
    def report_plot ():
        bed_data = bed_report ()
        if len (bed_data):
            fig, axs = plt.subplots (2, 1, figsize = (16,4), sharex = True)
            reporting_plot (bed_data, axs, chrom_sizes())
            for model in model_presets().keys():
                axs[0].plot ((),(), lw = 10, color = colorsCN[model], label = model)
            axs[0].plot ((),(), lw = 10, color = 'yellow', label = 'complex')
            #axs[0].plot ((),(), lw = 10, color = 'red', label = 'fail')
            axs[0].legend (bbox_to_anchor = (0.5, 2), ncol = len(model_presets())+2,
                           loc = 'upper center', title = 'Models of mixed clones: normal (AB) and abnormal karyotypes:')
            return fig
    
    @output
    @render.plot (alt = 'Diploid segments plot')
    def diploid_segments_plot():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))

            
            filt = (bed_data['ai'] < Consts.DIPLOID_AI_THR) &\
                   (np.abs(bed_data['cn']-2) < Consts.DIPLOID_dCN_THR)
            
            dip_bed_data = bed_data.loc[filt]            
            
            n_median = np.median (dip_bed_data['n'])
            
            if len(dip_bed_data):
                check_solution_plot_opt (dip_bed_data, ax, model_thr = input.model_thr(),
                                         highlight = input.chroms_selected())
                
            
            def ellipse (par, **kwargs):
                return Ellipse ((par['m_cn']+2, par['m_ai']), 2*d*par['s_cn'], 2*d*par['s_ai'], **kwargs)
                        
            p = 10**(-par_d['thr_HE'])
            d = sts.norm.ppf(1-p, par_d['m_d'], par_d['s_d'])/np.sqrt(n_median)
            if np.isfinite (d):
                ax.add_patch (ellipse(par_d, label = 'auto thr',
                              lw = 1, fill = False, color = 'r', ls = ':'))
            
            p = 10**(-input.HE_thr())
            d = sts.norm.ppf(1-p, par_d['m_d'], par_d['s_d'])/np.sqrt(n_median)
            if np.isfinite(d):
                ax.add_patch (ellipse(par_d, label = 'user thr', 
                              lw = 1, fill = False, color = 'b', ls = '--') )
            ax.legend()
            ymin, ymax = ax.get_ylim()
            ax.set_xlim (0.95*(2-Consts.DIPLOID_dCN_THR), 1.05*(2+Consts.DIPLOID_dCN_THR))
            ax.set_ylim (max((-0.01*Consts.DIPLOID_AI_THR, ymin)),
                         min((1.01*Consts.DIPLOID_AI_THR, ymax)))

            return fig
        
    @output
    @render.plot (alt = 'Diploid distance distribution plot')
    def diploid_distance_plot():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (figsize = (6,6))
            plot_cdf (bed_data['d_HE'].values, ax, par = (par_d['m_d'],par_d['s_d']),
                      all_colors = np.array([colorsCN[m] for m in bed_data['model']]), half = True)
            ax.set_ylabel ('cdf - HE distance')
            return fig

    @output
    @render.plot (alt = 'Model segments plot')
    def models_plot():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            
            k = np.linspace (0,1,100)
            m0 = par_d['m0']
            for model in model_presets().keys():
                if model in par()['models']:
                    ls = '-'
                else:
                    ls = ':'
                ax.plot (2*model_presets()[model].m(k, m0)/m0, model_presets()[model].ai(k, m0),  
                         lw = 1, linestyle = ls, color = colorsCN[model], alpha  = 1)

            check_solution_plot_opt (bed_data, ax, model_thr = input.model_thr(),
                                     highlight = input.chroms_selected())
            
            ax.set_xlim (2*0.9*bed_data.m.min()/m0, 2*1.1*bed_data.m.max()/m0)
            ax.set_ylim ((max(-0.02, -0.02*bed_data.ai.max()), bed_data.ai.max()*1.1))
            
                       
            return fig
        
    @output
    @render.plot (alt = 'Model distance distribution plot')
    def model_distance_plot():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            tmp_bed = bed_data.loc[[m not in ['AB', '(AB)(2+n)', '(AB)(2-n)'] for m in bed_data['model']]].sort_values (by = 'd_model')
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            ax.scatter (tmp_bed['d_model'].values, np.linspace (0,1, len(tmp_bed)),
                            c = np.array([colorsCN[m] for m in tmp_bed['model']]),
                            s = np.sqrt(tmp_bed['size']))
            x = np.linspace (tmp_bed['d_model'].min(), tmp_bed['d_model'].max(), 100)
            ax.plot (x , 1 - np.exp (par_d['a_d'] * x), 'r-')
            ax.set_ylabel ('cdf - Model distance')  
            return fig

    
    
    @output
    @render.plot (alt = "Scoring view")
    def scoring_plots ():
        bed_data = bed()
        par_d = par()
        
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, axs = plt.subplots (2, 1, figsize = (6,6))
                        
            plot_cdf (bed_data['d_HE'].values, axs[0], par = (par_d['m_d'],par_d['s_d']),
                      all_colors = np.array([colorsCN[m] for m in bed_data['model']]), half = True)
            axs[0].set_ylabel ('cdf - HE distance')

            tmp_bed = bed_data.loc[bed_data['model'] != 'AB'].sort_values (by = 'd_model')
            axs[1].scatter (tmp_bed['d_model'].values, np.linspace (0,1, len(tmp_bed)),
                            c = np.array([colorsCN[m] for m in tmp_bed['model']]),
                            s = np.sqrt(tmp_bed['size']))
            x = np.linspace (tmp_bed['d_model'].min(), tmp_bed['d_model'].max(), 100)
            axs[1].plot (x , 1 - np.exp (par_d['a_d'] * x), 'r-')
            axs[1].set_ylabel ('cdf - Model distance')  
            
            fig.tight_layout()
            
            return fig
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            k = np.linspace (0,1,100)
            m0 = par_d['m0']
            for model in model_presets().keys():
                if model in par()['models']:
                    ls = '-'
                else:
                    ls = ':'
                ax.plot (2*model_presets()[model].m(k, m0)/m0, model_presets()[model].ai(k, m0),  
                         lw = 1, linestyle = ls, color = colorsCN[model], alpha  = 1)
            
            check_solution_plot_opt (bed_data, ax, 
                                     highlight = input.chroms_selected(),
                                     model_thr = input.model_thr())
            
            ax.set_xlim (2*0.9*bed_data.m.min()/m0, 2*1.1*bed_data.m.max()/m0)
            ax.set_ylim ((max(-0.02, -0.02*bed_data.ai.max()), bed_data.ai.max()*1.1))
            
            return fig
    

    
    @output
    @render.plot(alt = "")
    def solution_plot_dipl_opt ():
        bed_data = opt_bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))

            
            filt = (bed_data['ai'] < Consts.DIPLOID_AI_THR) &\
                   (np.abs(bed_data['cn']-2) < Consts.DIPLOID_dCN_THR)
            
            dip_bed_data = bed_data.loc[filt]            
            
            n_median = np.median (dip_bed_data['n'])
            
            if len(dip_bed_data):
                check_solution_plot_opt (dip_bed_data, ax, model_thr = np.inf,
                                         highlight = input.chroms_selected())
                
            
            def ellipse (par, **kwargs):
                return Ellipse ((par['m_cn']+2, par['m_ai']), 2*d*par['s_cn'], 2*d*par['s_ai'], **kwargs)
                        
            p = 10**(-par_d['thr_HE'])
            d = sts.norm.ppf(1-p, par_d['m_d'], par_d['s_d'])/np.sqrt(n_median)
            if np.isfinite (d):
                ax.add_patch (ellipse(par_d, label = 'auto thr',
                              lw = 1, fill = False, color = 'r', ls = ':'))
            
            p = 10**(-input.HE_thr())
            d = sts.norm.ppf(1-p, par_d['m_d'], par_d['s_d'])/np.sqrt(n_median)
            if np.isfinite(d):
                ax.add_patch (ellipse(par_d, label = 'user thr', 
                              lw = 1, fill = False, color = 'b', ls = '--') )
            ax.legend()
            ymin, ymax = ax.get_ylim()
            ax.set_xlim (0.95*(2-Consts.DIPLOID_dCN_THR), 1.05*(2+Consts.DIPLOID_dCN_THR))
            ax.set_ylim (max((-0.01*Consts.DIPLOID_AI_THR, ymin)),
                         min((1.01*Consts.DIPLOID_AI_THR, ymax)))

            return fig

        
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot_opt ():
        opt_bed_data = opt_bed()
        par_d = par()
        if (len(opt_bed_data) != 0) & (len(par_d.keys()) != 0):
           
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            
            k = np.linspace (0,1,100)
            for model in model_presets().keys():

                if model in input.models_selected():

                    ls = '-'
                else:
                    ls = ':'
                ax.plot (model_presets()[model].m(k, m0_opt()), model_presets()[model].ai(k, m0_opt()), 
                         lw = 1, linestyle = ls, color = colorsCN[model], alpha = 0.6, label = model)
            
            check_solution_plot_opt (opt_bed_data, ax, model_thr = np.inf,
                                          highlight = [], xcol = 'm')
            
            ax.legend (bbox_to_anchor = (1.2,1), loc = 'upper center')
               
            return fig
    
    @reactive.Effect
    @reactive.event(solutions_list)
    def _():
        if len(solutions_list().keys()):
            ms = np.array(list(solutions_list()))
            diff = np.abs(ms - m0_opt())
            i = np.where(diff == diff.min())[0][0]
            value = ms[i]
            ui.update_slider ('m0_cov', value = value,
                              min = min(solutions_list().keys()),
                              max = max(solutions_list().keys()),
                              step = input.step() )
    
    @reactive.Effect
    @reactive.event(par)
    def _():
        if 'thr_model' in par().keys():

            if np.isfinite(par()['thr_model']):
                ui.update_slider ('model_thr', value = par()['thr_model'],
                                  min = 0,
                                  #max = 10, 
                                  max = max([10, int(par()['thr_model'])*5]),
                                  step = 0.1)


    @reactive.Effect
    @reactive.event(par)
    def _():
        if 'thr_HE' in par().keys():
            if np.isfinite(par()['thr_HE']):
                ui.update_slider ('HE_thr', value = par()['thr_HE'],
                                  min = 0,
                                  #max = 10,
                                  max = max([10, int(par()['thr_HE'])*5]), 
                                  step = 0.1)


    @output
    @render.table
    def chrom_segments ():
        bed_data = bed()
        if (len(bed_data) != 0):
            return bed_data.loc[bed_data.chrom == input.chrom_view()].sort_values(by = 'start')

    @output
    @render.table
    def report():
        report = bed_report().copy()
        if len(report) > 0:
            omit_models = ['AB', '(AB)(2+n)', '(AB)(2-n)']
            for m in input.report_models():
                omit_models.remove (m)
            return report.loc[[m not in omit_models for m in report.model.tolist()]][['chrom', 'start', 'end', 'cn',
                                                                                      'model', 'score_model', 'k', 'cyto']]
        
    @output
    @render.table
    def CNVs ():
        bed_data = bed()
        if (len(bed_data) != 0):
            
            if input.sort_CNV_by() == 'score':

                tmp_bed = bed_data.loc[bed_data['score_HE'] > input.HE_thr()].sort_values(by = 'score_HE', ascending = False)
            else:
                tmp_bed = bed_data.loc[bed_data['score_HE'] > input.HE_thr()]

            return tmp_bed

    @output
    @render.text
    def number_CNVs():
        bed_data = bed()
        if len(bed_data) > 0: 

            message = "Number of CNVs found: " + str(len(bed_data.loc[bed_data['score_HE'] > input.HE_thr()]))

        else:
            message = ''
        return  message

    @output
    @render.text
    def auto_model_thr():
        par_d = par()
        if 'thr_model' in par_d.keys():
            message = "Auto model thr: " + '{:.2f}'.format(par_d['thr_model']) 
        else:
            message = ''
        return message
        
    @output
    @render.text
    def auto_HE_thr():
        par_d = par()
        if 'thr_HE' in par_d.keys():
            message = "Auto HE thr: " + '{:.2f}'.format(par_d['thr_HE']) 
        else:
            message = ''
        return message
    
    
    @output
    @render.plot
    def data_plot ():
        bed_data = opt_bed()
        par_d = par()
        data_df = data()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0) & (len(data_df) != 0):
            fig, axs = plt.subplots (4, 1, figsize = (12,4), sharex = True)
            earth_worm_plot (data_df, bed_data, par_d, input.chrom_view(), axs, 
                             max_score_HE = input.HE_thr(), model_thr = input.model_thr())

            return fig
        
    @output
    @render.plot
    def compare_plot():
        bed_data = opt_bed()
        data_df = data()
        if (len(bed_data) != 0) & (len(data_df) != 0):
            CNV_bed = bed_data.loc[bed_data.chrom == input.chrom_view()]
            data_chrom = data_df.loc[data_df.chrom == input.chrom_view()]
            fig, ax = plt.subplots (figsize = (6,8))
            verification_plot_CNV (data_chrom, CNV_bed, ax, par(), input.f_to_plot(), column = 'vaf')
            return fig
    
    @output
    @render.plot
    def cov_plot():
        bed_data = opt_bed()
        data_df = data()
        if (len(bed_data) != 0) & (len(data_df) != 0):
            CNV_bed = bed_data.loc[bed_data.chrom == input.chrom_view()]
            data_chrom = data_df.loc[data_df.chrom == input.chrom_view()]
            fig, ax = plt.subplots (figsize = (6,8))
            verification_plot_CNV (data_chrom, CNV_bed, ax, par(), column = 'cov')
            return fig
    
        
    ###optimalization below    
    @reactive.Effect
    @reactive.event (input.opt)
    def _():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
           
            ms = np.arange (int(m0()*input.min_cn()),
                            int(m0()*input.max_cn())+input.step(),
                            input.step())
            
            ai = np.array(bed_data['ai'].values)

            l = len (ai)
            m_cov = np.array(bed_data['m'].values)
            sizes = np.array(bed_data['size'].values)
            n = np.array (bed_data['n'].values)

            
            ai_index = ai < Consts.DIPLOID_AI_THR
            solutions = {}
            
            with ui.Progress (min = ms[0], max = ms[-1]) as p:
                p.set(message = 'Optimizing solution...', )
                                            
                for m in ms:
                    p.set(m, message = 'Calculating')
                    cn = 2*m_cov/m
                    cn_index = np.abs(cn -2) < Consts.DIPLOID_dCN_THR
                    index = np.where(ai_index & cn_index)[0]
                    if len(index) > 2:

                        data_for_scoring = np.concatenate([ai[index], cn[index]-2, n[index]]).reshape (3,len(index)).T
                        print(len(index))
                        #par()['m0']
                        scorer = Scoring.Scoring(fb = par()['fb'], m0 = m, window_size = Consts.SNPS_IN_WINDOW, 
                                                 initial_data = data_for_scoring)
                        
                    
                        m_ai = scorer.ai_param['m']
                        s_ai = scorer.ai_param['s']
                        m_cn = scorer.cn_param['m']                   
                        s_cn = scorer.cn_param['s']
                    
                        d_HE = np.sqrt (((ai-m_ai)/s_ai)**2 + ((cn-2-m_cn)/s_cn)**2)
                        p_d = sts.norm.sf (d_HE, scorer.dipl_dist['m'], scorer.dipl_dist['s'])
                        thr = FDR (np.sort(p_d[np.isfinite(p_d)]), Consts.DIPLOID_ALPHA)
                        
                    else:
                        scorer = Scoring.Scoring()
                        p_d = np.zeros (l)
                        thr = 0
                                       
                    models = []
                    d_model = np.repeat(np.nan, l).astype(float)
                        
                    for i, pd in enumerate (p_d):
                        try:
                            sm = Models.pick_model(ai[i], s_ai, cn[i], s_cn, input.models_selected())
                        except (IndexError, AssertionError, UnboundLocalError):
                            sm = {'model' : 'UN', 'd_model' : np.nan,
                                  'k': 1, 'p_model' : np.nan,}
                        
                        models.append(sm['model'])
                        d_model[i] = sm['d_model']
                        
                    d_total = np.nansum((d_model*sizes))    
                    solutions[m] = (scorer, d_total/sizes.sum(), models)
            

            solutions_list.set (solutions)    
                    
       
    @output
    @render.plot (alt = 'Total distance to the model')
    def opt_plot ():
        
        if len(solutions_list().keys()) > 0:
            fig,ax = plt.subplots (figsize = (6,6))
            
            solutions = solutions_list()
            

            ms = solutions.keys()

            d_total = np.array([solutions[m][1] for m in solutions.keys()])
            d_HE = np.array([solutions[m][0].dipl_dist['m'] for m in solutions.keys()])
             
            fig, ax = plt.subplots(1,1, figsize = (6,6))
            ax.plot (ms, d_total, 'ro-')
            axt = ax.twinx()
            axt.plot (ms, d_HE, 'bo-')
            
            ax.set_xlabel ('Covearage')
            ax.set_ylabel ('Distance to the model')
            ax.yaxis.label.set_color ('r')
            ax.tick_params(axis='y', colors = 'r')
            
            
            axt.set_ylabel ('Diploid shift')
            axt.yaxis.label.set_color('b')
            axt.tick_params(axis='y', colors = 'b')
            
            # ax.legend()
            return fig
               
    @reactive.Effect
    @reactive.event (input.m0_cov)
    def _():
        if len(solutions_list()):
            
            m0_opt.set(input.m0_cov())
            bed = opt_bed().copy()
            bed['cn'] = 2*bed['m']/m0_opt()
            solution = solutions_list()[m0_opt()]
            print (solution)
            bed['model'] = solution[-1]
            opt_bed.set(bed)
            
app = App(app_ui, server, debug=True)

def FDR (p, alpha, score = False):
    k = np.arange (1, len(p)+1)
    index = np.where(p <= alpha*k/len(p))[0]
    try:
        if score:
            return -np.log10(p[np.max(index)])
        else:
            return p[np.max(index)]
    except:
        return np.inf
