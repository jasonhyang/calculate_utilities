from .calculate_dependencies import *
class calculate_smoothingFunctions():
    def __init__(self):
        self.data=[];

    # smoothing functions
    def fit_trajectories(self,x_I,y_I,fit_func_I='lowess',plot_textLabels_I=None,plot_fit_I=False):
        '''fit trajectory growth rate data to a smoothing function'''
        #Input:
        #   x_I = ale_time
        #   y_I = growth_rate
        #Output:
        #   x_O = ale_time_fitted
        #   y_O = growth_rate_fitted

        #cnt = 1;
        x = [];
        y = [];
        x = x_I;
        y = y_I;
        if fit_func_I=='spline':
            #spline
            tck = splrep(x,y,k=3,s=.025) #no smoothing factor
            #tck = splrep(x,y,k=3,task=-1,t=10) #no smoothing factor
            x2 = linspace(min(x),max(x),500)
            y2_spline= splev(x2,tck)
            y2 = numpy.zeros_like(y2_spline);
            for i,y2s in enumerate(y2_spline):
                if i==0:
                    y2[i] = y2s;
                elif i!=0 and y2s<y2[i-1]:
                    y2[i] = y2[i-1];
                else:
                    y2[i] = y2s;
        elif fit_func_I=='movingWindow':
            #moving window filter
            x2 = numpy.array(x);
            y2 = smooth(numpy.array(y),window_len=10, window='hanning');
        elif fit_func_I=='legendre':
            #legendre smoothing optimization
            smooth = legendre_smooth(len(x),1,1e-4,25)
            x2 = numpy.array(x);
            y2 = smooth.fit(numpy.array(y))
        elif fit_func_I=='lowess':
            #lowess
            x2 = numpy.array(x);
            y2_lowess = lowess.lowess(x2,numpy.array(y),f=0.1,iter=100)
            y2 = numpy.zeros_like(y2_lowess);
            for i,y2s in enumerate(y2_lowess):
                if i==0:
                    y2[i] = y2s;
                elif i!=0 and y2s<y2[i-1]:
                    y2[i] = y2[i-1];
                else:
                    y2[i] = y2s;
        else:
            print("fit function not recongnized");
        if plot_fit_I:
            ##QC plot using MatPlotLib
            # Create a Figure object.
            fig = pp.figure();
            # Create an Axes object.
            ax = fig.add_subplot(1,1,1) # one row, one column, first plot
            ## Add a title.
            #ax.set_title(k['sample_label'])
            # Set the axis
            pp.axis([0,max(x),0,max(y)+0.1]);
            # Add axis labels.
            ax.set_xlabel('Time [days]')
            ax.set_ylabel('GR [hr-1]')
            ## Label data points
            #tck = splrep(x,y,k=3,s=1.); #spline fit with very high smoothing factor
            #x_days = ALEsKOs_textLabels[k['sample_name_abbreviation']]['day']
            #y_days = splev(x_days,tck)
            #for i,txt in enumerate(ALEsKOs_textLabels[k['sample_name_abbreviation']]['dataType']):
            #    ax.annotate(txt, (x_days[i],y_days[i]-.15))
            # Create the plot
            #pp.plot(x_days,y_days,'rx',x,y,'b.',x2,y2,'g')
            pp.plot(x,y,'b.',x2,y2,'g')
            #display the plot
            pp.show()
        #record
        x_O = [];
        y_O = [];
        x_O = x2;
        y_O = y2;
        #cnt += 1;
        return x_O, y_O;