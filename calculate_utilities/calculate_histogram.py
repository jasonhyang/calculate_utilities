from .calculate_dependencies import *
class calculate_histogram():
    def __init__(self):
        self.data=[];

    # histogram and kde
    def histogram(self, data_I, n_bins_I=50, calc_bins_I=True):
        '''generate lower bound of the bins, the bin widths, and frequency of the data'''

        x_O = []; #the lower bound of the bin (inclusive)
        dx_O = []; #the width of the bin; x + dx is the upper bound (exclusive).
        y_O = []; #the count (if frequency is true), or the probability (if frequency is false).

        if calc_bins_I:
            n_bins = ceil(sqrt(len(data_I)));
        else:
            n_bins = n_bins_I;

        hist = numpy.histogram(data_I,n_bins);
        y_O = hist[0];
        edges = hist[1];

        for i in range(len(edges)-1):
            x_O.append(edges[i])
            dx_O.append(edges[i+1]-edges[i])

        return x_O,dx_O,y_O
    def pdf_kde(self,data_I,min_I=None,max_I=None,points_I=1000,bandwidth_I=None):
        '''compute the pdf from the kernal density estimate'''

        if min_I and max_I:
            min_point = min_I;
            max_point = max_I;
        else:
            min_point = min(data_I);
            max_point = max(data_I);

        x_grid = numpy.linspace(min_point, max_point, 1000)
        try:
            kde_scipy=scipy.stats.gaussian_kde(data_I, bw_method=bandwidth_I);
        except RuntimeError as e:
            print(e)
            return [0],[0];
        pdf = kde_scipy(x_grid);

        return x_grid,pdf