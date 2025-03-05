import pandas as pd
import numpy as np
import math

import multiprocessing as mpl
import scipy.stats as sts
import sys
import sklearn.linear_model as slm
import hmmlearn.hmm as hmm
import ruptures as rpt
from ruptures.exceptions import BadSegmentationParameters

from S2K import Testing
from S2K import Chromosome
from S2K import Consts
from S2K import Report
from S2K import Scoring
from S2K import Consts

class Genome:
    """Class to run genome wide tests of HE and create chromosomes."""
    def __init__(self, sample_name, logger, config, CB_file, models, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        self.models = models
        self.CB = pd.read_csv (CB_file, sep = '\t')
        self.logger = logger.getChild (f'{self.__class__.__name__}')
        self.genome_medians = {}
        self.logger.debug ("Object genome created.")
    
    def retrive_counts_create_chromosomes (self, data_file, columns, SG_file = None):
        """Reads the data in and filters through SuperGood list, if not None"""
        alldata = pd.read_csv (data_file, sep = '\t', usecols = columns)[columns]
        alldata.columns = ['chrom','position','ref_count', 'alt_count', 'Type']
        alldata['chrom'] = np.where(alldata['chrom'].str.contains("chr") == False, "chr"+alldata['chrom'], alldata['chrom'])
        
        self.logger.debug ('Data retrived.') 
        if SG_file is not None:
            SGlist = pd.read_csv (SG_file, compression = 'gzip',
                                header = None, sep = '\t', names = ['chrom','position','Ref', 'Alt', 'Status'],
                                dtype = {'chrom' : str, 'position' : np.int32, 'Ref' : str, 'Alt' : str, 'Status' : str})
            self.data = alldata.loc [alldata['Type'] == 'SNP'].merge (SGlist[['chrom','position']], how = 'inner')
            del (SGlist)
            del (alldata)
            self.logger.info (f"Input file filtered through SG list: {len(self.data)} SNPs to analyze.")
            
        else:
            self.data = alldata
            self.logger.info (f"{len(self.data)} SNPs to analyze. No filtering applied.")

        self.data['cov'] = self.data['ref_count'] + self.data['alt_count']
        self.data['vaf'] = self.data['alt_count']/self.data['cov']
        
        self.chromosomes = {}
        self.sex_chromosomes = {}
        
        if len(self.data.loc[~self.data['vaf'].isna()]) == 0:
            sys.exit (f'No data read in. Please cross check column names: {columns} with input file.')
            
        self.logger.debug ("Creating chromosomes...")
        for chrom, data in self.data.loc[~self.data['vaf'].isna()].groupby (by = 'chrom'):
            if chrom not in Consts.SEX_CHROMS:
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                             self.config, self.logger,
                                                             self.genome_medians, 
                                                             self.CB.loc[self.CB['chrom'] == chrom],
                                                             self.models)
            else:
                self.sex_chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                             self.config, self.logger,
                                                             self.genome_medians, 
                                                             self.CB.loc[self.CB['chrom'] == chrom],
                                                             self.models)
            self.logger.debug (f"Chromosome {chrom} has {len(data)} markers.")
         
    def segment_genome (self, m0 = 0, mc = 5.0, fb_alpha = Consts.FB_ALPHA):
        
        self.logger.debug ('Starting testing ...')
        
        self.HE = Testing.Testing ('HE', 
                                   self.chromosomes,
                                   self.logger)
        
        self.HE.run_test(no_processes = self.no_processes)

        self.logger.debug ("Genomewide heterozygosity: " + "\n" + str(self.HE.report_results()))
                
        outliers = self.HE.results.loc[self.HE.results['chi2'] > float(self.config['HE']['max_chi2'])].index.tolist()
        outliers_fraction = len(outliers)/len(self.HE.results)
        if len(outliers):
            self.logger.info (f'HE outliers due to chi2: {outliers}.')
        if outliers_fraction == 1:
            self.logger.critical ('Sample failed HE model. All chromosomes chi2 above threshold.')
            exit (1)
        elif outliers_fraction > 0.5:
            self.logger.warning ('Half or more of the chromosomes above threshold. May be inaccurate.')
               
        self.HE.analyze (parameters = self.config['HE'], outliers = outliers, skip_par = ['cov'])
        self.logger.info ("Genome heterozygosity reference: "+ "\n" + str(self.HE.get_genome_medians()))
        self.genome_medians['HE'] = self.HE.get_genome_medians()
        
        self.logger.debug ('Running N/E marking.')
        
        for chrom in self.chromosomes.keys():
            status = self.HE.get_status (chrom)
            params = self.HE.get_parameters(chrom).copy()
            for  par in status.index.values:
                if ~status[par]:
                    params[par] = self.HE.get_genome_medians()[par]  
            
            self.chromosomes[chrom].markE_onHE (params, float(self.config['HE']['z_thr']))
            
        #for formatting
        for  par in status.index.values:
            params[par] = self.HE.get_genome_medians()[par]  

        
        for chrom in self.sex_chromosomes.keys():
            self.sex_chromosomes[chrom].markE_onHE (params, float(self.config['HE']['z_thr']))
        
        self.logger.debug ('Testing N/E marking.')    
        
        self.VAF = Testing.Testing ('VAF', 
                                    self.chromosomes,
                                    self.logger)
        
        self.VAF.run_test (self.HE.medians['cov'], no_processes = self.no_processes)
        
        self.logger.debug ("Genomewide VAF: " + "\n" + str(self.VAF.report_results()))
        
        
        VAFS_chi2 = self.VAF.results['chi2']
        thr_chi2 = float(self.config['VAF']['max_chi2'])

        outliers = self.VAF.results.loc[(VAFS_chi2 > thr_chi2)|(VAFS_chi2 == 0)].index.tolist()
        outliers_fraction = len(outliers)/len(self.HE.results)
        if len(outliers):
            self.logger.info (f'VAF outliers due to chi2: {outliers}.')
        if outliers_fraction == 1:
            self.logger.critical ('All chromosomes failed VAF model. All chromosomes chi2 above threshold.')
            exit (1)
        elif outliers_fraction > 0.5:
            self.logger.warning ('Half or more of the chromosomes above threshold. May be inaccurate.')
        self.VAF.analyze (parameters = self.config['VAF'], outliers = outliers)
        
        self.logger.info ("Genome VAF reference: "+ "\n" + str(self.VAF.get_genome_medians()))
        self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        for chrom in self.chromosomes.keys():
            status = self.VAF.get_status (chrom)
            #self.logger.debug (f'Chromosome {chrom} inlier: {status}')
            if ~status.T.all(axis = 0):
                params = self.VAF.get_parameters (chrom)
                self.chromosomes[chrom].mark_on_full_model (self.HE.medians['cov']) 

                self.logger.debug (f'Chromosome {chrom} marked on full model.')
                
        for chrom in self.sex_chromosomes.keys():
            self.logger.info (f'Analyzing {chrom}: ...')
            chrom_vaf_results = Testing.VAF_test (self.sex_chromosomes[chrom].data,
                                                     self.HE.medians['cov'])
            self.logger.info (f'VAF test results: {chrom_vaf_results}')
            self.VAF.results.append(pd.DataFrame.from_records ([chrom_vaf_results],
                                                                columns = chrom_vaf_results._fields, 
                                                                index = [chrom]))
        
            self.chromosomes[chrom] = self.sex_chromosomes[chrom]
            self.logger.info(f'Chromosome {chrom} added.')        

        
        self.COV = Testing.Testing ('COV', 
                                    self.chromosomes,
                                    self.logger)
        self.COV.run_test(no_processes = self.no_processes,
                          exclude_symbol = [Consts.N_SYMBOL, Consts.U_SYMBOL])
        self.logger.debug ('Genomewide COV: ' + "\n" + str(self.COV.report_results()))
        
        self.COV.analyze (parameters = self.config['COV'], outliers = self.VAF.get_outliers()+Consts.SEX_CHROMS)
        
        self.logger.info ("Genome COV reference: " + "\n" + str(self.COV.get_genome_medians()))
        self.genome_medians['COV'] = self.COV.get_genome_medians()

        
        inliers = self.VAF.get_inliers()
        inliers_fb = self.VAF.results.loc[inliers,'fb'].values
        
        if len (np.unique(inliers_fb)) < 4:
            self.genome_medians['fb'] = np.percentile (inliers_fb,1-fb_alpha)
            self.logger.warning(f'Widening parameter estimation based on normal approximation not possible. {1-fb_alpha} percentile used.')
        else:
            try:
                self.genome_medians['fb'] = Testing.get_outliers_thrdist (np.unique(inliers_fb), fb_alpha, r = 0.5)[1]
            except:
                self.logger.exception ('Estimation of fb failed.')
                exit (1)
  
        self.genome_medians['v0'] = self.genome_medians['VAF']['vaf']
        if m0 > 0:
            self.logger.info (f"Using user supplied m0 = {m0}, instead of estimated m0 = {self.genome_medians['COV']['m']}")
            self.genome_medians['m0'] = m0
        else:
            self.genome_medians['m0'] = self.genome_medians['COV']['m'] 

        self.genome_medians['l'] = self.genome_medians['COV']['l']
        
        if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
            self.logger.critical (f"Coverage is below threshold {self.COV.medians['m']} < {self.config['COV']['min_cov']}")
            exit (1)

        self.logger.debug ("Starting segmentation.")
        if self.no_processes > 1:            
            with mpl.Pool (processes = self.no_processes) as pool:
                segmented_chroms = pool.map (f, self.chromosomes.values()) # function at bottom of this script
                for sc in segmented_chroms:
                    self.chromosomes [sc.name] = sc
        else:
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].generate_segments (1e6)
        self.logger.debug ("Segmentation finished.") 
        self.logger.info (f"Median coverage for the sample: m = {str(self.genome_medians['COV']['m'])}")
        self.logger.info ("Scoring segments.")
        
        self.all_chrs = []
        self.all_segments = []
        self.all_wrongs = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                self.all_chrs.append(chrom)
                self.all_segments.append(seg)
            for wrong in self.chromosomes[chrom].wrong_segments:
                self.all_wrongs.append(wrong)
        
        self.logger.info('Total number of segments is {}'.format(len(self.all_wrongs)+len(self.all_segments)))
        self.logger.info('Number of valid segments is {}'.format(len(self.all_segments)))
        
        self.diploid_segment_filter = [(s.centromere_fraction < Consts.CENTROMERE_THR) &\
                                        (s.parameters['ai'] < Consts.DIPLOID_AI_THR) &\
                                        (np.abs(s.parameters['m']/self.genome_medians['m0']-1) < Consts.DIPLOID_dCN_THR/2)
                                        for s in self.all_segments] 
                        
        self.diploid_chrs = [self.all_chrs[i] for i in np.where(self.diploid_segment_filter)[0]]
        self.diploid_segments = [self.all_segments[i] for i in np.where(self.diploid_segment_filter)[0]]

        self.logger.info('Total number of diploid segments is {}'.format(len(self.diploid_segments)))
        
        filtered_segments = []
        for s in self.diploid_segments:
            cns = s.parameters['m']/s.genome_medians['m0']
            if cns > 0.5 and cns < 1.5:
                filtered_segments.append(s)
        self.logger.info(('Total number of diploid segments after copy number and # of SNVs filtering is {}'.format(len(filtered_segments))))
        
        """
        new debug stuff 
        """
        # # TODO newly added shit
        # self.diploid_segment_filter = [
        #         (s.centromere_fraction < Consts.CENTROMERE_THR) &  # Centromere condition must be met
        #         (
        #             (s.parameters['ai'] < Consts.DIPLOID_AI_THR) |  # Either AI condition
        #             (np.abs(s.parameters['m'] / self.genome_medians['m0'] - 1) < Consts.DIPLOID_dCN_THR / 2)  # Or dCN condition
        #         )
        #         for s in self.all_segments
        #     ]
        
        # # Convert to NumPy array for further processing
        # self.diploid_segment_filter = np.array(self.diploid_segment_filter)
        
        # # Step 2: Identify which segments pass the filter
        # passed_filter_indices = np.where(self.diploid_segment_filter)[0]
        
        # self.logger.info(f"Number of segments passing the filter: {len(passed_filter_indices)}")
        
        # self.diploid_chrs = [self.all_chrs[i] for i in passed_filter_indices]
        # self.diploid_segments = [self.all_segments[i] for i in passed_filter_indices]
        
        # self.logger.info('Total number of diploid segments is {}'.format(len(self.diploid_segments)))
        
        # # Step 3: Further filtering for copy number and SNVs
        # filtered_segments = []
        # for s in self.diploid_segments:
        #     cns = s.parameters['m'] / s.genome_medians['m0']
        #     cn_condition = 0.5 < cns < 1.5
        
        #     # Log details of copy number filtering
        #     self.logger.info(f"Segment {s}:")
        #     self.logger.info(f"  - Copy number condition: {cn_condition} (cns: {cns})")
        
        #     if cn_condition:
        #         filtered_segments.append(s)
        
        # self.logger.info('Total number of diploid segments after copy number and SNVs filtering is {}'.format(len(filtered_segments)))
        
        # # Step 4: Prepare data for scoring
        # if len(filtered_segments) > 0:
        #     data_for_scoring = np.array([(s.parameters['ai'], 
        #                                   2 * s.parameters['m'] / self.genome_medians['m0'] - 2, 
        #                                   s.parameters['n']) for s in filtered_segments])
        #     self.scorer = Scoring.Scoring(fb=self.genome_medians['fb'], 
        #                                   m0=self.genome_medians['m0'], 
        #                                   window_size=Consts.SNPS_IN_WINDOW, 
        #                                   initial_data=data_for_scoring, 
        #                                   logger=self.logger)
        # else:
        #     self.logger.warning("No segments passed all filtering criteria. Skipping scoring step.")
            
            
        #     self.scorer = None
        # # TODO end of newly added shit
        
        
        data_for_scoring = np.array([(s.parameters['ai'], 2*s.parameters['m']/self.genome_medians['m0']-2, s.parameters['n']) for s in filtered_segments])
        self.scorer = Scoring.Scoring (fb = self.genome_medians['fb'], m0 = self.genome_medians['m0'], window_size = Consts.SNPS_IN_WINDOW, 
                                       initial_data = data_for_scoring, logger = self.logger)
        
        ps = np.zeros (len(self.all_segments))
        for i, seg in enumerate (self.all_segments):
            self.scorer.score_dipl(seg)
            ps[i] = seg.parameters['p_HE']
        
        thr = FDR(np.sort(ps[np.isfinite(ps)]), alpha = Consts.DIPLOID_ALPHA, score = True)
        self.genome_medians['thr_HE'] = thr
        self.scorer.set_thr (thr)
        for seg in self.all_segments:
            self.scorer.analyze_segment(seg, self.models)
        self.score_model_distance ()
        
        """
        Merge process
        """
        if mc > 0:
            # Merging segments
            self.logger.info('Merging segments...')
            merged = self.merge_segments_cpd() 
            
            self.logger.debug ("Starting merging segmentation.")      
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].generate_merged_segments(merged[chrom])
            self.logger.debug ("Merging segmentation finished.") 
            
            self.all_merged_chrs = []
            self.all_merged_segments = []
            self.all_merged_wrongs = []
            for chrom in self.chromosomes.keys():
                for seg in self.chromosomes[chrom].merged_segments:
                    self.all_merged_chrs.append(chrom)
                    self.all_merged_segments.append(seg)
                for wrong in self.chromosomes[chrom].wrong_segments:
                    self.all_merged_wrongs.append(wrong)
            
            self.logger.info('Total number of merged segments is {}'.format(len(self.all_merged_wrongs)+len(self.all_merged_segments)))
            self.logger.info('Number of valid merged segments is {}'.format(len(self.all_merged_segments)))
                    
            ps = np.zeros (len(self.all_merged_segments))
            for i, seg in enumerate (self.all_merged_segments):
                self.scorer.score_dipl(seg)
                ps[i] = seg.parameters['p_HE']
            
            thr = FDR(np.sort(ps[np.isfinite(ps)]), alpha = Consts.DIPLOID_ALPHA, score = True)
            self.genome_medians['thr_HE'] = thr
            self.scorer.set_thr (thr)
            for seg in self.all_merged_segments:
                self.scorer.analyze_segment(seg, self.models)
            self.score_model_distance (merged=True)      
        
            # Postprocesing merged segments for new potential merging
            post_merge_indicator = 1
            merge_cnt = 0
            while post_merge_indicator > 0:
                new_merged = self.post_merge(mc = mc)
                self.logger.info('Merging round {}'.format(merge_cnt+1))
                self.logger.debug ("Starting postprocessing segmentation.")      
                for chrom in self.chromosomes.keys():
                    self.chromosomes[chrom].generate_merged_segments(new_merged[chrom])
                self.logger.debug ("Postprocessing finished.") 
        
                self.logger.info ("Scoring final segments.")
            
                self.all_merged_chrs = []
                self.all_merged_segments = []
                self.all_merged_wrongs = []
                for chrom in self.chromosomes.keys():
                    for seg in self.chromosomes[chrom].merged_segments:
                        self.all_merged_chrs.append(chrom)
                        self.all_merged_segments.append(seg)
                    for wrong in self.chromosomes[chrom].wrong_segments:
                        self.all_merged_wrongs.append(wrong)
                
                self.logger.info('Total number of merged segments is {}'.format(len(self.all_merged_wrongs)+len(self.all_merged_segments)))
                self.logger.info('Number of valid merged segments is {}'.format(len(self.all_merged_segments)))
                        
                ps = np.zeros (len(self.all_merged_segments))
                for i, seg in enumerate (self.all_merged_segments):
                    self.scorer.score_dipl(seg)
                    ps[i] = seg.parameters['p_HE']
                
                thr = FDR(np.sort(ps[np.isfinite(ps)]), alpha = Consts.DIPLOID_ALPHA, score = True)
                self.genome_medians['thr_HE'] = thr
                self.scorer.set_thr (thr)

                for seg in self.all_merged_segments:
                    self.scorer.analyze_segment(seg, self.models)
                self.score_model_distance (merged=True)
                
                post_merge_indicator = self.merge_check(mc) # TODO
                merge_cnt += 1
                if merge_cnt == 10:
                    break
            
            """
            Ad hoc fix for chromosome 19 low coverage issue
            """
            flag = 0
            for i, seg in enumerate(self.all_merged_segments):
                if seg.chrom != 'chr19':
                    continue
                else:
                    newm = seg.parameters['m'] + 2 * 3.178160862432513
                    if seg.parameters['model'] == '(AB)(2-n)':
                        self.all_merged_segments[i].parameters['m'] = newm
                        self.logger.info('replaced new coverage {} to {}, {}'.format(seg.parameters['m'], self.all_merged_segments[i].parameters['m'], newm))
                        flag = 1
                        
            if flag == 1:
                for seg in self.all_merged_segments:
                    self.scorer.analyze_segment(seg, self.models)
                self.score_model_distance (merged=True)
            
        else:
            self.logger.info('No merging is performed')
            self.all_merged_segments = self.all_segments
            self.all_merged_chrs = self.all_chrs
            self.all_merged_wrongs = self.all_wrongs
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].merged_segments = []
                for seg in self.chromosomes[chrom].segments:
                    self.chromosomes[chrom].merged_segments.append(seg)
                    
        """
        Final typing with clustering
        """

        bed = pd.DataFrame()
        for segment in self.all_merged_segments:
            if (1-segment.centromere_fraction)*(segment.end - segment.start) < 5e6 or segment.parameters['model'] == 'UN':
                continue
            else:
                attr_list = [segment.chrom, 
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
                             segment.parameters['AB']] + \
                            [segment.parameters[x] for x in segment.models] + \
                            [segment.parameters['model_confidence'],
                            segment.parameters['score_model'], 
                            segment.parameters['k'], 
                            segment.cytobands,
                            segment.centromere_fraction]
                bed = bed.append(pd.Series(attr_list), ignore_index=True)
                
        bed.columns = [
            'chrom', 'start', 'end', 'length', 'ai', 'n', 'm', 'cn', 'd_HE', 
            'score_HE', 'model', 'd_model', 'AB'
        ] + \
        segment.models + \
        [
            'model_confidence', 'score_model', 'k', 'cytobands', 'centromere_fraction'
        ]
        
        if len(bed) > 20:
            
            from sklearn.preprocessing import StandardScaler
            from sklearn.mixture import GaussianMixture
            
            self.logger.info('Calibration with clustering method')

            bed_filtered = bed[((1-bed['centromere_fraction']) * bed['length'] >= 5e6)]
            scaler = StandardScaler()
            X = scaler.fit_transform(bed_filtered[['cn', 'ai']].values)
            gmm = GaussianMixture(
                                    n_components=10,
                                    covariance_type='diag',
                                    init_params='kmeans',
                                    n_init=10,
                                    tol=1e-4,
                                    max_iter=500,
                                    random_state=0
                                ).fit(X)
            labels = gmm.predict(X)
            
            self.logger.info(f"Generated labels: {np.unique(labels)}")

            # Step 4: Add labels to the filtered dataframe
            self.logger.info('Adding labels to filtered dataframe...')
            bed_filtered = bed_filtered.copy()  # Ensure bed_filtered is independent of bed
            bed_filtered['cluster'] = labels
            self.logger.info(f"First few rows of filtered data with clusters:\n{bed_filtered[['cn', 'ai', 'cluster']].head()}")
            self.logger.info(f"Size of bed_filtered: {bed_filtered.shape}")
            
            # Step 5: Initialize cluster column in the original dataframe
            self.logger.info('Initializing cluster column in the original dataframe...')
            bed['cluster'] = np.nan
            self.logger.info(f"Initial 'cluster' column in bed:\n{bed['cluster'].head()}")
            self.logger.info(f"Size of bed: {bed.shape}")
            
            # Step 6: Map labels back to the original dataframe
            self.logger.info('Mapping labels back to the original dataframe...')
            bed.loc[bed_filtered.index, 'cluster'] = bed_filtered['cluster']
            self.logger.info(f"Updated 'cluster' column in bed:\n{bed[['centromere_fraction', 'length', 'cluster']].head()}")
            self.logger.info(f"Size of bed after mapping clusters: {bed.shape}")

            for cluster_label, group in bed.groupby('cluster'):
                # Add verbose logging
                self.logger.info(f"Processing cluster label: {cluster_label}")
                
                # Log the specific columns for the group
                specific_columns = ['ai', 'cn', 'model']
                if all(col in group.columns for col in specific_columns):
                    self.logger.info(f"Details for cluster {cluster_label} (specific columns):\n{group[specific_columns]}")
                else:
                    self.logger.warning(f"Some of the specific columns ({specific_columns}) are missing in the group for cluster {cluster_label}.")
                
                # Calculate model sums and log them
                models_sum = group[['AB'] + segment.models].sum()
                self.logger.info(f"Model sums for cluster {cluster_label}:\n{models_sum}")
                
                # Determine the calibrated model and log it
                calibrated_model = models_sum.idxmin()
                self.logger.info(f"Calibrated model for cluster {cluster_label}: {calibrated_model}")
                
                # Update the 'calibrated_model' column in the original dataframe
                bed.loc[bed['cluster'] == cluster_label, 'calibrated_model'] = calibrated_model

                
            for i, row in bed.iterrows():
                if pd.isna(row['cluster']):
                    continue
                if row['calibrated_model'] == 'A' and row['cn'] > 2:
                    continue
                chrom = row['chrom']
                start = row['start']
                end = row['end']
            
                # Find the index of the matching segment
                for idx, segment in enumerate(self.chromosomes[chrom].merged_segments):
                    if segment.start == start and segment.end == end:
                        # Update the model and k parameters in the original segment object
                        self.logger.info(f"{row['cluster']}")
                        self.logger.info(f"{segment.chrom} {segment.start} {segment.end}")
                        self.logger.info(f"{row['model']} to {row['calibrated_model']}")
                        self.chromosomes[chrom].merged_segments[idx].parameters['model'] = row['calibrated_model']
                        self.chromosomes[chrom].merged_segments[idx].parameters['k'] = Models.calc_k(
                            segment.parameters['ai'], 
                            self.scorer.ai_param['s'], 
                            2 * segment.parameters['m'] / segment.genome_medians['m0'], 
                            self.scorer.cn_param['s'], 
                            row['calibrated_model']
                        )
                        break
            
        # """
        # special cases with homozygous deletion in TSG
        # """
        # # locations of homozygous deletion
        # HDs = Consts.HD
        # # find segments within those locations
        # indices, hd_tsg = find_hds(self.all_merged_chrs, self.all_merged_segments, preprocess_coordinates(HDs))
        # if len(hd_tsg) > 0:
        #     self.logger.info('Found homozygous deletion at {}'.format(hd_tsg))
        #     for i, seg in zip(indices, hd_tsg):
        #         prev = self.all_merged_segments[i-1]
        #         succ = self.all_merged_segments[i+1]
        #         self.logger.error(prev)
        #         self.logger.error(succ)

        #         if prev.parameters['model'] == 'A' or prev.parameters['model'] == 'AA' or prev.parameters['model'] == 'AAA' or prev.parameters['model'] == 'AAAA':
        #             a = 1 - prev.parameters['k'] # AB
        #             x = (4 * a * seg.parameters['ai'])/(1 - 2 * seg.parameters['ai']) # A
        #             y = 1 - a - x # (AB)(2-n)
        #         elif succ.parameters['model'] == 'A' or succ.parameters['model'] == 'AA' or succ.parameters['model'] == 'AAA' or succ.parameters['model'] == 'AAAA':
        #             a = 1 - succ.parameters['k'] # AB
        #             x = (4 * a * seg.parameters['ai'])/(1 - 2 * seg.parameters['ai']) # A
        #             y = 1 - a - x # (AB)(2-n)
        #         else:
        #             self.logger.error('Cannot re-estimate homozygous deletion. something wrong with adjacent model selection')
        #         self.all_merged_segments[i].parameters['k'] = y
        
    def score_model_distance (self, merged=False):
        
        if merged:
            all_segments = self.all_merged_segments
        else:
            all_segments = self.all_segments
        
        size_filter = np.array([(seg.end - seg.start)/1e6 > Consts.SIZE_THR  for seg in all_segments])
        cent_filter = np.array([seg.centromere_fraction < Consts.CENTROMERE_THR  for seg in all_segments])
        model_filter = np.array([seg.parameters['model'] not in ['AB', '(AB)(2-n)', '(AB)(2+n)'] for seg in all_segments])
        finite_filter = np.array([np.isfinite(seg.parameters['d_model']) for seg in all_segments])
        
        filter = size_filter & cent_filter & model_filter & finite_filter
        models = self.models

        if sum(filter) >= 3:
            indexes = np.where(filter)[0]
            segments = [all_segments[i] for i in indexes]
            zs_ns = [seg.parameters['d_model'] for seg in segments]
            z_n = np.array(zs_ns)
            x = sts.expon.ppf (np.linspace (0,1,len(z_n)+2)[1:-1])
            huber = slm.HuberRegressor (fit_intercept = False)
            huber.fit (x[:, np.newaxis], np.sort(z_n))
            a = -1./huber.coef_[0]
            self.logger.info ('Distance from model /d/ distribution: FI(d) = exp({:.5f} d)'.format (a))
            ps = []
            for seg in all_segments:
                if seg.parameters['model'] != 'AB':
                    p = np.exp (a*seg.parameters['d_model'])
                    seg.parameters['score_model'] = -np.log10 (p)
                    ps.append(p)
                else:
                    seg.parameters['score_model'] = 0
                    
                model_scores = np.array([seg.parameters[m] for m in models])
                model_index = np.argmin(model_scores)
                confidences = softmax(model_scores)
                confidence_score = confidences[model_index]
                seg.parameters['model_confidence'] = confidence_score
                
            self.genome_medians['d_model'] = {'a' : a}
            ps = np.array(ps)
            self.genome_medians['thr_model'] = FDR (np.sort(ps[np.isfinite(ps)]), Consts.MODEL_APLHA, score = True)
        else:
            self.logger.info ('Not enough non diploid regions to perform meaningful scoring')
            self.genome_medians['d_model'] = {'a' : np.nan}
            for seg in all_segments:
                seg.parameters['score_model'] = 0
                seg.parameters['model_confidence'] = 0
            self.genome_medians['thr_model'] = np.nan
    
    def merge_segments_cpd(self):
        merged_segments = {}
        self.logger.info('Merging segments using CPD...')

        for c in Consts.CHROM_ORDER:
            merged_segments[c] = []
            self.logger.debug('currently merging chromosome {}'.format(c))
            try:
                tmp = self.chromosomes[c].segments
            except KeyError:
                continue
            tmp_x = [[s.parameters['m'], s.parameters['AB']] + [s.parameters[x] for x in self.models] for s in tmp]
                
            tmp_x = np.array(tmp_x)
            tmp_x = np.nan_to_num(tmp_x, nan=-1)
                        
            algo = rpt.KernelCPD(kernel='linear', min_size=1, jump=1).fit(tmp_x)
            # algo = rpt.BottomUp(model='l1', min_size=1, jump=1).fit(tmp_x)
            breaks = np.linspace(1, 20, 20)
            
            best_cost = float('inf')
            best_bkps = None
            for bkp in breaks:
                try:
                    result = algo.predict(n_bkps=bkp)
                    cost = algo.cost.sum_of_costs(result)
                    
                    if cost < best_cost:
                        best_cost = cost
                        best_bkps = bkp
                except BadSegmentationParameters:
                    continue
            self.logger.info('Best break points {} best score {}'.format(best_bkps, best_cost))
            try:
                result = algo.predict(n_bkps=best_bkps)
            except (BadSegmentationParameters, AssertionError):
                self.logger.error('BadSegmentationParameters error at merging chromosome {}'.format(c))
                if len(tmp) == 0:
                    continue
                else:
                    merged_segments[c].append((tmp[0].start, tmp[-1].end))
                continue
            result = np.array(result) - 1
            
            for i in range(len(result)):
                if i == 0:
                    tobemerged = (tmp[0].start, tmp[result[i]].end)
                else:
                    tobemerged = (tmp[result[i-1]+1].start, tmp[result[i]].end)
                merged_segments[c].append(tobemerged)
        
        return merged_segments
    
    def post_merge(self, mc):
        new_merges = {}
    
        self.logger.info('Postprocessing the CPD merged segments...')
    
        for c in Consts.CHROM_ORDER:
            new_merges[c] = []
    
            self.logger.debug('currently processing chromosome {}'.format(c))
            try:
                tmp = self.chromosomes[c].merged_segments
            except KeyError:
                continue
            tmp_x = [[s.start, 
                      s.end,
                      s.end - s.start,
                      s.parameters['ai'],
                      2 * s.parameters['m'] / s.genome_medians['m0'],  # CN
                      s.parameters['model']] for s in tmp]
    
            ai_th = min(mc * self.scorer.ai_param['s'], 0.1) # merging coefficient
            cn_th = min(mc * self.scorer.cn_param['s'], 0.5) # merging coefficient
    
            current_start = None
            current_end = None
            current_ai = None
            current_cn = None
    
            # Iterate over each row in the DataFrame
            for index, item in enumerate(tmp_x):
                start, end, _, ai, cn, model = item
                if current_start is None:
                    # Initialize current values on the first iteration
                    current_start = start
                    current_end = end
                    current_ai = ai
                    current_cn = cn
                else:
                    # Check if conditions are met to merge the current segment with the previous one
                    if abs(ai - current_ai) < ai_th and abs(cn - current_cn) < cn_th:
                        # Extend the current segment
                        current_end = end
                        current_ai = ai
                        current_cn = cn
                    else:
                        # Save the current segment and start a new one
                        new_merges[c].append((current_start, current_end))
                        current_start = start
                        current_end = end
                        current_ai = ai
                        current_cn = cn
    
            # Save the last segment if there is one
            if current_start is not None:
                new_merges[c].append((current_start, current_end))
    
        return new_merges
    
    def merge_check(self, mc):
        total_merge = []
        ai_th = min(mc * self.scorer.ai_param['s'], 0.1) # merging coefficient
        cn_th = min(mc * self.scorer.cn_param['s'], 0.5) # merging coefficient
        for c in Consts.CHROM_ORDER:
            merge_status = []
            try:
                tmp = self.chromosomes[c].merged_segments
            except KeyError:
                merge_status.append(False)
                continue
            
            tmp_x = [[s.start, 
                      s.end,
                      s.end - s.start,
                      s.parameters['ai'],
                      2 * s.parameters['m'] / s.genome_medians['m0'],  # CN
                      s.parameters['model']] for s in tmp]
            for i, s in enumerate(tmp_x):
                s1 = tmp_x[i + 1] if i < len(tmp_x) - 1 else None
                if s1 is None:
                    diff1 = False
                else:
                    ai = s[3]
                    cn = s[4] 
                    ai1 = s1[3]
                    cn1 = s1[4]  
                    diff1 = (abs(cn - cn1) < cn_th) & (abs(ai - ai1) < ai_th)
                
                merge_status.append(diff1)
            merge_cnt = sum(1 for item in merge_status if item == True)
            if merge_cnt > 0:
                total_merge.append(1)
            else:
                total_merge.append(0)
        return sum(total_merge)
    
    def report (self, report_type = 'bed'):
        return Report.Report(report_type).genome_report(self)

def merge_contiguous_numbers(states):
    merges = []
    start_index = 0
    for i in range(1, len(states)):
        if states[i] != states[i-1]:
            merges.append((start_index, i-1))
            start_index = i
    merges.append((start_index, len(states)-1))
    return merges

# Define the function to regularize the covariance matrix
def regularize_covariance(covs, min_regularization=1e-4):
    if covs.ndim == 3:  # full covariance matrices
        for matrix in covs:
            np.fill_diagonal(matrix, matrix.diagonal() + min_regularization)
    elif covs.ndim == 2:  # diag covariance matrices
        covs += min_regularization
    return covs

def find_elbow(scores):
    # Calculate the differences between consecutive scores
    diffs = np.diff(scores)

    # Calculate the rate of change between the differences
    rate_change = np.diff(diffs)

    # The elbow point is where the rate of change is maximized
    elbow_idx = np.argmax(rate_change) + 1  # +1 as diff reduces the length by 1

    return elbow_idx

#k, m, s must have correct dimentions
def score_double_gauss (k, m, s):
    z = (k - m)/s
    ps =np.concatenate((sts.norm.cdf (z[:,0])[:, np.newaxis], sts.norm.sf (z[:,1])[:,np.newaxis]), axis = 1)
    p = np.min(ps , axis = -1)
    return  np.sort(p)

def fit_huber (data, alpha):
    k = np.log10 (data[:,0])
    s = np.log10 (data[:,1])
    huber = slm.HuberRegressor(alpha = 0.0, epsilon = 1.35)
    huber.fit(s[:, np.newaxis], k)

    A = -huber.coef_[0]
    B = 1
    C = -huber.intercept_
    d = (A*s+B*k+C)/np.sqrt (A**2+B**2)
    
    
    down, up = Testing.get_outliers_thrdist (d, alpha = alpha)
    inlier_ds = d[(d > down)&(d < up)]
    m, std = sts.norm.fit (inlier_ds)
    std = std/Bolch_correction (len(inlier_ds))

    score_FDR = FDR (np.sort(sts.norm.sf (d, m, std)), alpha, score = True)

    return {'A' : A, 'B' : B, 'C' : C, 'down' : down, 'up' : up, 'm' : m,
            's' : std, 'score_FDR' : score_FDR}

def Bolch_correction (n):
    return 1 - 1/(4*n) - 7/(32*n**2) - 19/(128*n**3)


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

def f (c):
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
def preprocess_coordinates(genomic_coordinates):
    sorted_coordinates = {}
    for chr, ranges in genomic_coordinates.items():
        sorted_ranges = sorted([ranges] if not isinstance(ranges, list) else ranges)
        sorted_coordinates[chr] = sorted_ranges
    return sorted_coordinates

def binary_search(ranges, segment):
    low, high = 0, len(ranges) - 1
    while low <= high:
        mid = (low + high) // 2
        start, end = ranges[mid]
        if end < segment.start:
            low = mid + 1
        elif start > segment.end:
            high = mid - 1
        else:
            return True  # Overlap found
    return False  # No overlap

def find_hds(chrs, segments, sorted_coordinates):
    overlapping_segments = []
    overlapping_indices = []
    i = 0
    for chrom, segment in zip(chrs, segments):
        if chrom in sorted_coordinates:
            if binary_search(sorted_coordinates[chrom], segment) and 2*segment.parameters['m']/segment.genome_medians['m0']<0.9 and segment.parameters['n']>1:
                overlapping_segments.append(segment)
                overlapping_indices.append(i)
        i += 1
    return overlapping_indices, overlapping_segments
