from .calculate_dependencies import *
class calculate_statistics():
    def __init__(self):
        self.data=[];

    # statistical analysis
    # calculate the geometric mean and variance:
    def calculate_ave_var_geometric(self,data_I):
        """ calculate the geometric average and var of data
        with 95% confidence intervals
        """

        try:
            data_ave_O = 0.0
            # calculate the average of the sample
            for c in data_I:
                data_ave_O += numpy.log(c);
            data_ave_O = numpy.exp(data_ave_O/len(data_I));

            data_var_O = 0.0
            #calculate the variance of the sample
            for c in data_I:
                data_var_O += numpy.power(numpy.log(c/data_ave_O),2);
            data_var_O = data_var_O/(len(data_I)-1); #note: we will need to take the exp()
                                                     # to get the final variance
                                                     # but leaving it this way makes the
                                                     # downstream calculations simpler

            #calculate the 95% confidence intervals
            data_se = numpy.sqrt(data_var_O/len(data_I));
            data_lb_O = numpy.exp(numpy.log(data_ave_O) - 1.96*data_se);
            data_ub_O = numpy.exp(numpy.log(data_ave_O) + 1.96*data_se);
            
            #correct the variance for use in reporting
            data_var_O = numpy.exp(data_var_O);

            return data_ave_O, data_var_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # calculate the mean and variance:
    def calculate_ave_var(self,data_I,confidence_I = 0.95):
        """calculate the average and var of data
        with 95% confidence intervals"""

        try:
            data = numpy.array(data_I);

            data_ave_O = 0.0
            # calculate the average of the sample
            data_ave_O = numpy.mean(data);

            data_var_O = 0.0
            #calculate the variance of the sample
            data_var_O = numpy.var(data);

            #calculate the standard error of the sample
            se = scipy.stats.sem(data)

            #calculate the 95% confidence intervals
            n = len(data);
            h = se * scipy.stats.t._ppf((1+confidence_I)/2., n-1)
            data_lb_O = data_ave_O - h;
            data_ub_O = data_ave_O + h;

            return data_ave_O, data_var_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # calculate the mean and CV:
    def calculate_ave_var_cv(self,data_I,confidence_I = 0.95):
        """calculate the average, var, %cv of data
        with 95% confidence intervals"""

        try:
            data = numpy.array(data_I);

            data_ave_O = 0.0
            # calculate the average of the sample
            data_ave_O = numpy.mean(data);

            data_var_O = 0.0
            #calculate the variance of the sample
            data_var_O = numpy.var(data);

            #calculate the standard error of the sample
            se = scipy.stats.sem(data)

            #calculate the CV% of the sample
            data_cv_O = 0.0;
            if data_ave_O !=0.0:
                data_cv_O = numpy.std(data)/data_ave_O*100;

            #calculate the 95% confidence intervals
            n = len(data);
            h = se * scipy.stats.t._ppf((1+confidence_I)/2., n-1)
            data_lb_O = data_ave_O - h;
            data_ub_O = data_ave_O + h;

            return data_ave_O, data_var_O, data_cv_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # convert CV to SD:
    def convert_cv2StDev(self,data_ave_I,data_cv_I):
        """convert the %CV to standard deviation
        INPUT
        data_ave_I = float, average/mean
        data_cv_I = float, % coefficient of variation or % relative standard deviation
        OUTPUT
        data_stdev_O = float, standard deviation"""

        try:
            data_stdev_O = 0.0;
            #calculate the SD of the sample
            if not data_ave_I is None and not data_cv_I is None:
                data_stdev_O = data_ave_I*data_cv_I/100.0;

            return data_stdev_O;
        except Exception as e:
            print(e);
            exit(-1);

    # calculate the interquartiles
    def calculate_interquartiles(self,data_I,iq_range_I = [25,75]):
        '''compute the interquartiles and return the min, max, median, iq1 and iq3'''

        min_O = numpy.min(data_I);
        max_O = numpy.max(data_I);
        iq_1_O, iq_2_O = numpy.percentile(data_I, iq_range_I)
        median_O = numpy.median(data_I);

        return min_O, max_O, median_O, iq_1_O, iq_2_O;
    # calculate the fold change
    def calculate_foldChange(self,data_1_I,data_2_I,type_I = 'relative', scale_values_I = None,scale_fold_change_I = None):
        """Calculate the fold change between two data value 1 and 2
        Input:
        data_1_I = data value 1
        data_2_I = data value 2
        type_I = 'relative', 'absolute','geometric' (relative is the default);
        scale_values_I = string, factor to scale values by prior to calculating the fold change
        scale_fold_change_I = string factor to scale the fold change by
        """
        # scale the values
        if data_1_I and data_2_I:
            fold_change_O = 0.0;
        elif data_1_I==0.0:
            print("data_1_I cannot be zero!");
            return None;
        else:
            print("bad data provided!");
            return None;
        if scale_values_I and scale_values_I == "log2":
            data_1 = numpy.log2(data_1_I);
            data_2 = numpy.log2(data_2_I);
        elif scale_values_I and scale_values_I == "log10":
            data_1 = numpy.log10(data_1_I);
            data_2 = numpy.log10(data_2_I);
        elif scale_values_I and scale_values_I == "ln":
            data_1 = numpy.log(data_1_I);
            data_2 = numpy.log(data_2_I);
        elif scale_values_I and scale_values_I == "abs":
            data_1 = numpy.abs(data_1_I);
            data_2 = numpy.abs(data_2_I);
        elif scale_values_I and scale_values_I == "exp":
            data_1 = numpy.exp(data_1_I);
            data_2 = numpy.exp(data_2_I);
        elif scale_values_I and scale_values_I == "exp2":
            data_1 = numpy.exp2(data_1_I);
            data_2 = numpy.exp2(data_2_I);
        elif scale_values_I and scale_values_I == "^10":
            data_1 = numpy.power(data_1_I,10);
            data_2 = numpy.power(data_2_I,10);
        elif scale_values_I and scale_values_I == "^2":
            data_1 = numpy.power(data_1_I,2);
            data_2 = numpy.power(data_2_I,2);
        elif scale_values_I:
            print("scale_values_I not recognized.  No scaling will be applied.");
            data_1 = data_1_I;
            data_2 = data_2_I;
        else:
            data_1 = data_1_I;
            data_2 = data_2_I;
        # relative, absolute, or geometric
        if type_I == 'relative':
            fold_change_O = data_2/data_1;
        elif type_I == 'absolute':
            fold_change_O = numpy.abs(data_2/data_1);
        elif type_I == 'geometric':
            fold_change_O = numpy.log(numpy.exp(data_2)/numpy.exp(data_1));
        else:
            fold_change_O = data_2/data_1;
            print("type_I not recognized.  Relative type will be applied as default.");
        # scale the fold change
        if scale_fold_change_I and scale_fold_change_I == "log2":
            fold_change_O = numpy.log2(fold_change_O);
        elif scale_fold_change_I and scale_fold_change_I == "log10":
            fold_change_O = numpy.log10(fold_change_O);
        elif scale_fold_change_I and scale_fold_change_I == "ln":
            fold_change_O = numpy.log(fold_change_O);
        elif scale_fold_change_I and scale_fold_change_I == "abs":
            fold_change_O = numpy.abs(fold_change_O);
        elif scale_fold_change_I and scale_fold_change_I == "sqrt":
            fold_change_O = numpy.sqrt(fold_change_O);
        elif scale_fold_change_I:
            print("scale_fold_change_I not recognized.  No scaling will be applied.");
            fold_change_O = fold_change_O; 
        else:
            fold_change_O = fold_change_O;
        return fold_change_O;
    # calculate the difference
    def calculate_difference(self,data_1_I,data_2_I,type_I = 'relative',scale_values_I = None,scale_difference_I = None):
        """Calculate the differencefold change between two data value 1 and 2
        Input:
        data_1_I = data value 1
        data_2_I = data value 2
        type_I = 'relative', 'absolute','geometric' (relative is the default);
        scale_values_I = string, factor to scale values by prior to calculating the difference
        scale_difference_I = string factor to scale the fold change by
        """
        if data_1_I and data_2_I:
            difference_O = 0.0;
        else:
            print("bad data provided!");
            return None;
        if scale_values_I and scale_values_I == "log2":
            data_1 = numpy.log2(data_1_I);
            data_2 = numpy.log2(data_2_I);
        elif scale_values_I and scale_values_I == "log10":
            data_1 = numpy.log10(data_1_I);
            data_2 = numpy.log10(data_2_I);
        elif scale_values_I and scale_values_I == "ln":
            data_1 = numpy.log(data_1_I);
            data_2 = numpy.log(data_2_I);
        elif scale_values_I and scale_values_I == "abs":
            data_1 = numpy.abs(data_1_I);
            data_2 = numpy.abs(data_2_I);
        elif scale_values_I and scale_values_I == "exp":
            data_1 = numpy.exp(data_1_I);
            data_2 = numpy.exp(data_2_I);
        elif scale_values_I and scale_values_I == "exp2":
            data_1 = numpy.exp2(data_1_I);
            data_2 = numpy.exp2(data_2_I);
        elif scale_values_I and scale_values_I == "^10":
            data_1 = numpy.power(data_1_I,10);
            data_2 = numpy.power(data_2_I,10);
        elif scale_values_I and scale_values_I == "^2":
            data_1 = numpy.power(data_1_I,2);
            data_2 = numpy.power(data_2_I,2);
        elif scale_values_I:
            print("scale_values_I not recognized.  No scaling will be applied.");
            data_1 = data_1_I;
            data_2 = data_2_I;
        else:
            data_1 = data_1_I;
            data_2 = data_2_I;
        # relative, absolute, or geometric
        if type_I == 'relative':
            difference_O = data_2-data_1;
        elif type_I == 'absolute':
            difference_O = numpy.abs(data_2-data_1);
        elif type_I == 'geometric': #i.e. distance
            difference_O = numpy.log(numpy.exp(data_2)+numpy.exp(data_1));
        else:
            difference_O = data_2-data_1;
            print("type_I not recognized.  Relative type will be applied as default.");
        
        if scale_difference_I and scale_difference_I == "log2":
            difference_O = numpy.log2(difference_O);
        elif scale_difference_I and scale_difference_I == "log10":
            difference_O = numpy.log10(difference_O);
        elif scale_difference_I and scale_difference_I == "ln":
            difference_O = numpy.log(difference_O);
        elif scale_difference_I and scale_difference_I == "abs":
            difference_O = numpy.abs(difference_O);
        elif scale_difference_I and scale_difference_I == "sqrt":
            difference_O = numpy.sqrt(difference_O);
        elif scale_difference_I:
            print("scale_difference_I not recognized.  No scaling will be applied.");
            difference_O = difference_O; 
        else:
            difference_O = difference_O;
        return difference_O;
