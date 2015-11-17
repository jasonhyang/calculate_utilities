from .calculate_dependencies import *
from math import hypot
from copy import copy
class calculate_correlation():
    '''
    Find correlations in a dataset
    Based on the methods described in 

    Analysis of strain and regional variation in gene expression in mouse brain
    Genome Biology 2001
    Paul Pavlidis and William S Noble
    '''
    def __init__(self):
        self.data=[];

    def convert_profileStr2List(self,profile_I):
        '''Convert a string profile to a list
        INPUT:
        profile_I = 0-1-2-3
        OUPTPUT
        profile_O = [0,1,2,3]
        '''
        profile_list = profile_I.split("-");
        profile_O = [];
        for p in profile_list:
            profile_O.append(int(p));
        return profile_O;
    def convert_profileList2Str(self,profile_I):
        '''Convert a list profile to string representation
        INPUT:
        profile_I = [0,1,2,3]
        OUPTPUT
        profile_O = 0-1-2-3
        '''
        pstr = [str(x) for x in profile_I]
        profile_O = '-'.join(pstr);
        return profile_O;
    def convert_data2trend(self,data_I,
                data_stdev_I=[],
                data_lb_I=[],data_ub_I=[],
                tol_I=1e-4,
                criteria_I='difference'):
        '''Convert a data list to a trend
        INPUT:
        data_I = list of data of mean values
        data_stdev_I = list of standard deviations
        data_lb_I = list of lower bounds
        data_ub_I = list of upper bounds
        tol_I = tolerance to consider a step-difference
        criteria_I = "difference" use the mean difference to determine if two data points are different
                     "stdev" use the +/- stdev to determine if two data points are different
                     "lb/ub" use the lb/ub to determine if two data points are different
                     Default = difference
        OUTPUT:
        trend_0 = list representation of a trend (e.g. [0,1,2,3])
        '''

        profile_tmp = [];
        prev_profile = 0;
        prev_value = 0.0
        #generate the template

        #generate the profile
        for i,d in enumerate(data_I):
            if i==0:
                profile_tmp.append(prev_profile);
                prev_value = d;
                continue;
            if criteria_I == 'difference':
                difference = d-prev_value;
                if numpy.abs(difference)>tol_I:
                    if difference > 0: prev_profile += 1;
                    else: prev_profile -= 1;
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            elif criteria_I == 'stdev':
                difference = d-prev_value;
                if prev_value - data_stdev_I[i-1] > d + data_stdev_I[i]:
                    prev_profile += 1;
                    profile_tmp.append(prev_profile);
                elif prev_value + data_stdev_I[i-1] < d - data_stdev_I[i]:
                    prev_profile -= 1;
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            elif criteria_I == 'lb/ub':
                difference = d-prev_value;
                if data_lb_I[i-1] > data_ub_I[i]:
                    prev_profile += 1;
                    profile_tmp.append(prev_profile);
                elif data_ub_I[i-1] < data_lb_I[i]:
                    prev_profile -= 1;
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            else:
                print("criteria not recognized.");
            #re-initialize the previous value
            prev_value = d;

        #normalize the profile
        trend_O = self.normalize_trend(profile_tmp);  
    def normalize_trend(self,trend_I):
        '''Normalize a trend to the range [0,inf)
        INPUT:
        trend_I = list representation of a trend w/ or w/o negative numbers
        OUTPUT:
        trend_O = list representation of a trend w/o negative numbers'''

        trend_O = profile_I;

        #get the distance between the min and max values
        min_val = min(trend_I);
        max_val = max(trend_I);
        dist = numpy.sqrt(max_val*max_val+min_val*min_val);
        dist=int(dist);

        #add the distance to each value
        if dist>0:
            trend_O = [dist+x for x in trend_I];
        else:
            trend_O = [0 for x in trend_I];
        return trend_O;
    def convert_data2profile(self,data_I,
                data_stdev_I=[],
                data_lb_I=[],data_ub_I=[],
                tol_I=1e-4,
                criteria_I='difference'):
        '''Convert a data list to a profile
        INPUT:
        data_I = list of data of mean values
        data_stdev_I = list of standard deviations
        data_lb_I = list of lower bounds
        data_ub_I = list of upper bounds
        tol_I = tolerance to consider a step-difference
        criteria_I = "difference" use the mean difference to determine if two data points are different
                     "stdev" use the +/- stdev to determine if two data points are different
                     "lb/ub" use the lb/ub to determine if two data points are different
                     Default = difference
        OUTPUT:
        profile_0 = list representation of a profile (e.g. [0,1,2,3])
        '''

        profile_tmp = [];
        prev_value = 0.0
        prev_profile = 0;
        #generate the template scale
        scale = self.convert_data2scaleDict(data_I);
        #generate the profile
        for i,d in enumerate(data_I):
            if i==0:
                prev_profile=scale[d];
                profile_tmp.append(scale[d]);
                prev_value = d;
                continue;
            if criteria_I == 'difference':
                difference = d-prev_value;
                if numpy.abs(difference)>tol_I:
                    prev_profile=scale[d];
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            elif criteria_I == 'stdev':
                difference = d-prev_value;
                if prev_value - data_stdev_I[i-1] > d + data_stdev_I[i]:
                    prev_profile=scale[d];
                    profile_tmp.append(prev_profile);
                elif prev_value + data_stdev_I[i-1] < d - data_stdev_I[i]:
                    prev_profile=scale[d];
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            elif criteria_I == 'lb/ub':
                difference = d-prev_value;
                if data_lb_I[i-1] > data_ub_I[i]:
                    prev_profile=scale[d];
                    profile_tmp.append(prev_profile);
                elif data_ub_I[i-1] < data_lb_I[i]:
                    prev_profile=scale[d];
                    profile_tmp.append(prev_profile);
                else:
                    profile_tmp.append(prev_profile);
            else:
                print("criteria not recognized.");
            #re-initialize the previous value
            prev_value = d;

        #normalize the profile (should be normalized by default)
        profile_O = profile_tmp;        

        return profile_O;
    def convert_data2scaleDict(self,data_I):
        '''convert list of numeric data to an integer scale
        INPUT:
        data_I = list of numerics
        OUTPUT
        scale_O = {} of data[i]:int
        '''
        scale_O = [];
        data = copy(data_I);
        data.sort(); #lowest to highest
        #make the data_dict
        scale_O = {};
        for i,d in enumerate(data):
            scale_O[d]=i;
        return scale_O;
       
    def calculate_correlation_pearsonr(self,profile_I = [],data_I = []):
        '''Check that the profile length matches the data length
        INPUT:
        profile_I = list representation of a profile (e.g. [0,1,2,3])
        data_I = list of data to correlate with the profile
        OUTPUT:
        rho_O = Pearson correlation coefficient
        pval_O = two-sided p-value for the hypthesis that the two sets of data are uncorrelated
        '''
        rho_O, pval_O = None,None;
        if self.check_profileAndDataLength(profile_I,data_I):
            rho_O, pval_O = scipy.stats.pearsonr(profile_I,data_I);
        return rho_O, pval_O;

    def calculate_correlation_spearmanr(self,profile_I = [],data_I = []):
        '''Check that the profile length matches the data length
        INPUT:
        profile_I = list representation of a profile (e.g. [0,1,2,3])
        data_I = list of data to correlate with the profile
        OUTPUT:
        rho_O = Spearman correlation matrix or correlation coefficient
        pval_O = two-sided p-value for the hypthesis that the two sets of data are uncorrelated
        '''
        rho_O, pval_O = None,None;
        if self.check_profileAndDataLength(profile_I,data_I):
            rho_O, pval_O = scipy.stats.spearmanr(profile_I,data_I);
        return rho_O, pval_O;

    def check_profileAndDataLength(self,profile_I,data_I):
        '''Check that the profile length matches the data length
        INPUT:
        profile_I = list representation of a profile (e.g. [0,1,2,3])
        data_I = list of data to correlate with the profile
        OUTPUT:
        check_O = boolean, true if the lengths are the same
        '''
        check_O = False;
        if len(profile_I)==len(data_I):
            check_O = True;
        else:
            print('profile and data lengths do not match');
        return check_O;
