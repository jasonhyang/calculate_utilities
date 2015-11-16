from .calculate_dependencies import *
class calculate_statisticsSampledPoints():
    def __init__(self):
        self.data=[];
    # calculate the confidence intervals
    def calculate_ciFromPoints(self,data_I, alpha=0.05):
        """Calculate the confidence intervals from sampled points"""
        data_sorted = numpy.sort(data_I)
        n = len(data_sorted)
        lb = data_sorted[int((alpha/2.0)*n)]
        ub = data_sorted[int((1-alpha/2.0)*n)]
        return lb,ub
    def bootstrap(self,data, num_samples=100000, statistic=numpy.mean, alpha=0.05):
        """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
        n = len(data)
        idx = npr.randint(0, n, (num_samples, n))
        samples = data[idx]
        stat = numpy.sort(statistic(samples, 1))
        return (stat[int((alpha/2.0)*num_samples)],
                stat[int((1-alpha/2.0)*num_samples)])
    # calculate the p-value difference
    def permutation_resampling(self,case, control, num_samples=50, statistic=numpy.mean):
        '''calculate the pvalue of two data sets using a resampling approach'''

        observed_diff = abs(statistic(case) - statistic(control))
        num_case = len(case)

        combined = numpy.concatenate([case, control])
        diffs = []
        for i in range(num_samples):
            xs = npr.permutation(combined)
            diff = numpy.mean(xs[:num_case]) - numpy.mean(xs[num_case:])
            diffs.append(diff)

        pval = (numpy.sum(diffs > observed_diff) +
                numpy.sum(diffs < -observed_diff))/float(num_samples)
        return pval, observed_diff, diffs
    def calculate_pvalue_permutation(self,data_1_I,data_2_I,n_permutations_I=10,n_resamples_I=10):    
        '''calculate the pvalue of two data by determining
        the lack of overlap between sample points using a permutation test.
        If the sample points of the data sets is not equal,
        a subset of samples of matching length is resampled from the larger data set'''
        
        data_1 = None; #sample set with fewer points
        data_2 = None; #sample set with more points
        n_resamples = 0;

        # check the length of data_1 and data_2
        if len(data_1_I)>len(data_2_I):
            data_1=numpy.array(data_2_I);
            data_2=numpy.array(data_1_I);
            n_resamples=n_resamples_I;
        elif len(data_1_I)<len(data_2_I):
            data_1=numpy.array(data_1_I);
            data_2=numpy.array(data_2_I);
            n_resamples=n_resamples_I;
        else:
            data_1=numpy.array(data_1_I);
            data_2=numpy.array(data_2_I);

        n_samples_min = len(data_1);

        vals = []
        for i in range(0,n_permutations_I):
            if n_resamples==0:
                cond1 = numpy.random.permutation(data_1)
                cond2 = numpy.random.permutation(data_2)
                z = cond1 - cond2
                x = len(z[z>0]) + 1
                y = len(z[z<0]) + 1
                k = min(x,y)
                vals.append(k)
            else:
                cond1 = numpy.random.permutation(data_1)
                cond2 = numpy.random.permutation(data_2)
                for resample in range(n_resamples):
                    cond2_int = numpy.random.randint(0,n_samples_min);
                    z = cond1 - cond2[cond2_int]
                    x = len(z[z>0]) + 1
                    y = len(z[z<0]) + 1
                    k = min(x,y)
                    vals.append(k)
        p = numpy.mean(vals)/len(data_1)*2
        return p;
