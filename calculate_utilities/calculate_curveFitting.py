from .calculate_dependencies import *
class calculate_curveFitting():
    def __init__(self):
        self.data=[];
    # linear regression
    def calculate_regressionParameters(self,concentrations_I,ratios_I,dilution_factors_I,fit_I,weighting_I,use_area_I):
        '''calculate regression parameters for a given component
        NOTE: intended to be used in a loop'''
        # input:
        #       concentrations_I
        #       ratios_I
        #       dilution_factors_I
        #       fit_I
        #       weighting_I
        #       use_area_I
        # ouput:
        #       slope
        #       intercept
        #       correlation
        #       lloq
        #       uloq
        #       points

        # note need to make complimentary method to query concentrations, ratios, and dilution factors
        # for each component prior to calling this function

        #TODO

        return