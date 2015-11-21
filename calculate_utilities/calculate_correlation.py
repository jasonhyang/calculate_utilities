from .calculate_dependencies import *
from math import hypot
from copy import copy
class calculate_correlation():
    '''
    Find correlations in a dataset
    Loosely based on the methods described in the following:
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
        
        #convert to lb/ub
        data_lb,data_ub = self.calculate_lbAndUb(data_I,
                data_stdev_I=data_stdev_I,
                data_lb_I=data_lb_I,
                data_ub_I=data_ub_I,
                tol_I=tol_I,
                criteria_I=criteria_I);
        #generate the profile
        profile_tmp=[0 for i,x in enumerate(data_I)];
        for i1,d1 in enumerate(data_I):
            if i1 == 0:
                profile_tmp[i1]=prev_profile;
                continue;
            # check if the profile is significantly different to any previous profile
            for i2,d2 in reversed(list(enumerate(data_I[0:i1]))):
                significant,distance,direction=self.calculate_LBUBDifference(data_lb[i2],data_ub[i2],data_lb[i1],data_ub[i1]);
                if significant and direction == '+':
                    prev_profile += 1;
                    profile_tmp[i1]=prev_profile;
                    break;
                elif significant and direction == '-':
                    prev_profile -= 1;
                    profile_tmp[i1]=prev_profile;
                    break;
                else:
                    profile_tmp[i1]=prev_profile;
        ##TODO: check in reverse order if it is significantly greater or less the previous sample
        #for i,d in enumerate(data_I):
        #    if i==0:
        #        profile_tmp.append(prev_profile);
        #        prev_value = d;
        #        continue;
        #    if criteria_I == 'difference':
        #        difference = d-prev_value;
        #        if numpy.abs(difference)>tol_I:
        #            if difference > 0: prev_profile += 1;
        #            else: prev_profile -= 1;
        #            profile_tmp.append(prev_profile);
        #        else:
        #            profile_tmp.append(prev_profile);
        #    elif criteria_I == 'stdev':
        #        difference = d-prev_value;
        #        if prev_value - data_stdev_I[i-1] > d + data_stdev_I[i]:
        #            prev_profile -= 1;
        #            profile_tmp.append(prev_profile);
        #        elif prev_value + data_stdev_I[i-1] < d - data_stdev_I[i]:
        #            prev_profile += 1;
        #            profile_tmp.append(prev_profile);
        #        else:
        #            profile_tmp.append(prev_profile);
        #    elif criteria_I == 'lb/ub':
        #        difference = d-prev_value;
        #        if data_lb_I[i-1] > data_ub_I[i]:
        #            prev_profile -= 1;
        #            profile_tmp.append(prev_profile);
        #        elif data_ub_I[i-1] < data_lb_I[i]:
        #            prev_profile += 1;
        #            profile_tmp.append(prev_profile);
        #        else:
        #            profile_tmp.append(prev_profile);
        #    else:
        #        print("criteria not recognized.");
            ##re-initialize the previous value
            #prev_value = d;

        #normalize the profile
        trend_O = self.normalize_trend(profile_tmp); 
        return trend_O; 
    def normalize_trend(self,trend_I):
        '''Normalize a trend to the range [0,inf)
        INPUT:
        trend_I = list representation of a trend w/ or w/o negative numbers
        OUTPUT:
        trend_O = list representation of a trend w/o negative numbers'''

        #get the distance between the min and max values
        min_val = min(trend_I);
        max_val = max(trend_I);
        dist = numpy.sqrt(max_val*max_val+min_val*min_val);
        dist=int(dist);

        #add the distance to each value
        trend_tmp = [];
        if dist>0:
            trend_tmp = [dist+x for x in trend_I];
        else:
            trend_tmp = [0 for x in trend_I];

        #check that the trend is normalized to zero
        trend_O = [];
        min_val = min(trend_tmp);
        if min_val == 0:
            trend_O = trend_tmp;
        else:
            trend_O = [x-min_val for x in trend_tmp];

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
        #convert to lb/ub
        data_lb,data_ub = self.calculate_lbAndUb(data_I,
                data_stdev_I=data_stdev_I,
                data_lb_I=data_lb_I,
                data_ub_I=data_ub_I,
                tol_I=tol_I,
                criteria_I="difference");
        #get the range rank dictionary
        range = self.convert_data2RangeDict(data_I,data_lb,data_ub);
        #convert to lb/ub
        data_lb,data_ub = self.calculate_lbAndUb(data_I,
                data_stdev_I=data_stdev_I,
                data_lb_I=data_lb_I,
                data_ub_I=data_ub_I,
                tol_I=tol_I,
                criteria_I=criteria_I);
        #generate the profile
        profile_tmp=[0 for i,x in enumerate(data_I)];
        for i1,d1 in enumerate(data_I):
            if i1 == 0:
                profile_tmp[i1]=range[d1];
                continue;
            # check if the profile is significantly different to any previous profile
            for i2,d2 in reversed(list(enumerate(data_I[0:i1]))):
                significant,distance,direction=self.calculate_LBUBDifference(data_lb[i2],data_ub[i2],data_lb[i1],data_ub[i1]);
                if significant:
                    #compare to the first data point
                    significant,distance,direction=self.calculate_LBUBDifference(data_lb[0],data_ub[0],data_lb[i1],data_ub[i1]);
                    if significant:
                        profile_tmp[i1]=range[d1];
                    else:
                        profile_tmp[i1]=range[data_I[0]];
                    break;
                else:
                    profile_tmp[i1]=range[d2];

        #normalize the profile (should be normalized by default)
        profile_O = self.normalize_profile(profile_tmp);        

        return profile_O;
    def normalize_profile(self,profile_I):
        '''Normalize a profile to the range [0,1,2,...,max) in increments of +1
        INPUT:
        profile_I = list representation of a profile w/ or w/o gaps between increments
        OUTPUT:
        profile_O = list representation of a profile w/o gaps between increments'''

        profile_O = profile_I;

        #get the distance between the min and max values
        min_val = min(profile_I);
        max_val = max(profile_I);
        #get the number of unique values
        unique_val = list(set(profile_I));
        unique_val.sort();
        unique_val_dict = {};
        for i,uv in enumerate(unique_val):
            unique_val_dict[uv]=i;
        #normalize the profile
        profile_O = [unique_val_dict[x] for x in profile_I];
        return profile_O;
    def convert_data2RankDict(self,data_I):
        '''convert list of numeric data to an integer rank
        INPUT:
        data_I = list of numerics
        OUTPUT
        rank_O = {} of data[i]:int
        '''
        rank_O = [];
        data = copy(data_I);
        data.sort(); #lowest to highest
        #make the data_dict
        rank_O = {};
        for i,d in enumerate(data):
            rank_O[d]=i;
        return rank_O;
    def calculate_lbAndUb(self,data_I,
                data_stdev_I=[],
                data_lb_I=[],data_ub_I=[],
                tol_I=1e-4,
                criteria_I='difference'):
        '''convert list of numeric data to an integer rank
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
        OUTPUT
        data_lb_O = list of lower bounds
        data_ub_O = list of upper bounds
        '''
        #calculate the lb/ub (if not given)
        data_lb_O = [];
        data_ub_O = [];
        for i,d in enumerate(data_I):
            if criteria_I == 'difference':
                data_lb_O.append(d-tol_I);
                data_ub_O.append(d+tol_I);
            elif criteria_I == 'stdev':
                data_lb_O.append(d-data_stdev_I[i]);
                data_ub_O.append(d+data_stdev_I[i]);
            elif criteria_I == 'lb/ub':
                data_lb_O.append(data_lb_I[i]);
                data_ub_O.append(data_ub_I[i]);
            else:
                print("criteria not recognized.");
        return data_lb_O,data_ub_O;
    def calculate_LBUBDifference(self,data_1_lb_I,data_1_ub_I,data_2_lb_I,data_2_ub_I):
        '''Calculate the difference between two upper and lower bounds
        INPUT:
        data_1_lb_I
        data_1_ub_I
        data_2_lb_I
        data_2_ub_I
        OUTPUT:
        significant_O = boolean, true if lb/ub do not overlab
        distance_O = distance between the lb/ub if significant
        direction_O = - if range 1 is greater than range 2
                      + if range 2 is greater than range 1
        '''
        significant_O = False;
        distance_O = 0.0;
        direction_O = None;
        if data_1_lb_I > data_2_ub_I:
            significant_O=True;
            distance_O = abs(data_2_ub_I-data_1_lb_I);
            direction_O = '-';
        elif data_1_ub_I < data_2_lb_I:
            significant_O=True;
            distance_O = abs(data_1_ub_I-data_2_lb_I);
            direction_O = '+';
        return significant_O,distance_O,direction_O;

    def convert_data2RangeDict(self,data_I,
                data_lb_I=[],data_ub_I=[]):
        '''convert list of numeric data to an integer rank
        INPUT:
        data_I = list of data of mean values
        data_lb_I = list of lower bounds
        data_ub_I = list of upper bounds
        OUTPUT
        range_O = {} of data[i]:int
        TODO:
        1. get all ranges
        2. determine if the range is different from any of the other ranges
        '''
        #identify the ranges
        ranges = [];
        ranges_lb = [];
        for d1_cnt,d1 in enumerate(data_I):
            lb = data_lb_I[d1_cnt];
            ub = data_ub_I[d1_cnt];
            for d2_cnt,d2 in enumerate(data_I):
                #check if the ranges are significantly different
                significant_O,distance_O,direction_O=self.calculate_LBUBDifference(lb,ub,data_lb_I[d2_cnt],data_ub_I[d2_cnt]);
                if significant_O: continue;
                else:
                    if data_ub_I[d2_cnt]>ub:
                        #expand the ub
                        ub = data_ub_I[d2_cnt];
                    if data_lb_I[d2_cnt]<lb:
                        #expand the lb
                        lb = data_lb_I[d2_cnt];
            if not (lb,ub) in ranges: 
                ranges.append((lb,ub));
                ranges_lb.append(lb);
        #assign data points and index to each range
        ranges_lb.sort();
        range_lb_dict = {};
        for i,lb in enumerate(ranges_lb):
            range_lb_dict[lb]=i;
        range_O  = {};
        for d in data_I:
            for i,range in enumerate(ranges):
                if d<=range[1] and d>=range[0]:
                    range_O[d]=range_lb_dict[range[0]];
                    break;            
        
        return range_O;
       
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
