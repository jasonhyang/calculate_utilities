from .calculate_dependencies import *

class calculate_biomass():
    def __init__(self):
        self.data=[];

    # calculations
    # biomass normalization
    def calculate_gdw_CVSAndCVSUnitsAndODAndConversionAndConversionUnits(self,cvs_I,cvs_units_I,od600_I,conversion_I,conversion_units_I):
        # check units
        if (cvs_units_I == 'mL' and conversion_units_I == 'gDW*L-1*OD600-1'):
            gdw_O = cvs_I*1e-3*od600_I*conversion_I;
            gdw_units_O = 'gDW';
            return gdw_O, gdw_units_O;
        else:
            print('biomass conversion units do not match!')
            exit(-1);
    def calculate_cellVolume_CVSAndCVSUnitsAndODAndConversionAndConversionUnits(self,cvs_I,cvs_units_I,od600_I,conversion_I,conversion_units_I):
        # check units
        if (cvs_units_I == 'mL' and conversion_units_I == 'uL*mL-1*OD600-1'):
            cellVolume_O = cvs_I*od600_I*conversion_I*1e-6;
            cellVolume_units_O = 'L';
            return cellVolume_O, cellVolume_units_O;
        else:
            print('cell volume conversion units do not match!')
            exit(-1);
    def calculate_conc_concAndConcUnitsAndDilAndDilUnitsAndConversionAndConversionUnits(self,conc_I,conc_units_I,dil_I,dil_units_I,conversion_I,conversion_units_I):
        # check units
        if (conc_units_I == 'uM' and conversion_units_I == 'L' and dil_units_I == 'mL'):
            conc_O = (conc_I*1e-6)*(dil_I)*(1/conversion_I);
            conc_units_O = 'mM';
            return conc_O, conc_units_O;
        elif (conc_units_I == 'uM' and conversion_units_I == 'gDW' and dil_units_I == 'mL'):
            conc_O = (conc_I*1e-3)*(dil_I)*(1/conversion_I);
            conc_units_O = 'umol*gDW-1';
            return conc_O, conc_units_O;
        elif ((conc_units_I == 'height_ratio' or conc_units_I == 'area_ratio') and (conversion_units_I == 'L' or conversion_units_I == 'gDW') and dil_units_I == 'mL'):
            conc_O = (conc_I*1e-3)*(dil_I)*(1/conversion_I);
            conc_units_O = conc_units_I;
            return conc_O, conc_units_O;
        else:
            print('concentration normalization units do not match!')
            exit(-1);
    def calculate_cultureDensity_ODAndConversionAndConversionUnits(self,od600_I,conversion_I,conversion_units_I):
        cultureDensity_O = od600_I*conversion_I;
        cultureDensity_units_O = conversion_units_I.replace('*OD600-1','');
        return cultureDensity_O, cultureDensity_units_O;
    def calculate_biomass_CVSAndCVSUnitsAndODAndConversionAndConversionUnits(self,cvs_I,cvs_units_I,od600_I,conversion_I,conversion_units_I):
        # check units
        if (cvs_units_I == 'mL' and conversion_units_I == 'gDW*L-1*OD600-1'):
            # return the biomass in gDW
            gdw_O = cvs_I*1e-3*od600_I*conversion_I;
            gdw_units_O = 'gDW';
            return gdw_O, gdw_units_O;
        elif (cvs_units_I == 'mL' and conversion_units_I == 'uL*mL-1*OD600-1'):
            # return the cell volume in L
            cellVolume_O = cvs_I*od600_I*conversion_I*1e-6;
            cellVolume_units_O = 'L';
            return cellVolume_O, cellVolume_units_O;
        else:
            print('biomass conversion units do not match!')
            exit(-1);
    def calculate_yield_growthRateAndUptakeRates(self,growthRate_ss_I,uptakeRate_ss_I,
                                                 growthRate_units_I=None,uptakeRate_units_I=None,
                                                 per_carbon_I=False):
        '''Calculate the biomass yield

        TODO: automatically calculate the output units

        Input:
        growthRate_ss_I = float, steady-state growth rate
        uptakeRate_ss_I = [float,float,...], steady-state uptake rates for (e.g., 1 or more carbon sources)
        growthRate_units_I = string, e.g., hr-1
        uptakeRate_units_I = string, e.g., mmol*gDCW*hr-1
        per_carbon_I = normalize uptake rates to a per carbon basis (default=False), not yet implemented
        
        Output:
        yield_ss_O = biomass yield 
        description of biomass yield: "The biomass yield (Y_X/S_ss) was calculated as the quotient of the growth rate and the glucose uptake rates during the exponential growth phase"
        yield_units_O = string, e.g.,mmol-1*gDCW-1
        '''

        # check input
        yield_ss_O = None;
        yield_units_O = None;
        if not growthRate_ss_I:
            print('no growth rate has been supplied')
            return yield_ss_O,yield_units_O;
        if not uptakeRate_ss_I:
            print('no uptake rates have been supplied')
            return yield_ss_O,yield_units_O;
        if None in uptakeRate_ss_I or (0.0 in uptakeRate_ss_I and len(uptakeRate_ss_I)<1):
            print('None or 0.0 values found in uptake rate')
            return yield_ss_O,yield_units_O;
        # calculate yield
        # TODO: implement on a per carbon basis
        if per_carbon_I:
            print('per carbon normaization not yet supported.  No normalization will be applied')
        yield_ss_O = growthRate_ss_I/numpy.sum(uptakeRate_ss_I);
        # determine the units
        if not growthRate_units_I or not uptakeRate_units_I:
            print('no growthRate or no uptakeRate units provided!')
            return yield_ss_O,yield_units_O;
        else:
            yield_units_O = "%s*(%s)*-1" %(growthRate_units_I,uptakeRate_units_I);
        return yield_ss_O,yield_units_O;
    
    def calculate_growthRate(self,time_I,biomass_I):
        '''calculate exponential growth'''

        x = numpy.array(time_I);
        y = numpy.log(biomass_I); #weight the biomass by the natural logarithmn

        slope, intercept, r_value, p_value, std_err = linregress(x,y);
        r2 = r_value**2; #coefficient of determination

        return slope, intercept, r2, p_value, std_err;
    def interpolate_biomass(self,time_I, slope, intercept):
        '''interpolate the biomass from an exponential fit of the growth rate'''

        biomass = [];
        for t in time_I:
            biomass.append(t*slope+intercept);
        return biomass;
    def calculate_uptakeAndSecretionRate(self,dcw_I,conc_I,gr_I):
        '''calculate uptake and secretion rates'''

        x = numpy.array(dcw_I);
        y = numpy.array(conc_I); 

        slope, intercept, r_value, p_value, std_err = linregress(x,y);
        r2 = r_value**2; #coefficient of determination

        rate = slope*gr_I;

        return slope, intercept, r2, p_value, std_err, rate;
    def interpolate_biomass(self,time_I, slope, intercept):
        '''interpolate the biomass from an exponential fit of the growth rate'''

        biomass = [];
        for t in time_I:
            biomass.append(t*slope+intercept);
        return biomass;