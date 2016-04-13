import mortdist
import mortdist_cpsh
import util
import patsy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.integrate as si
import scipy.stats as ss
import itertools

#legacy import
import DeDayMLEcore as core



#_x = np.linspace(-10, 10, 10)
#_y = np.linspace(-15, 0, 10)
#_grid2 = np.array(list(itertools.product(_x, _x)))
#_grid3 = np.array(list(itertools.product(_x, _x, _y)))
#_grid4 = np.array(list(itertools.product(_x, _x, _y, _x)))

#search grid
#Functionality to add, user defineable search grid TODO
#17**4=85312, likely a over-kill
_x = np.linspace(-10, 10, 17)
_y = np.linspace(-15, 0, 17)
_z = np.linspace(-20, 10, 17)
_grid2 = np.array(list(itertools.product(_x, _x)))
_grid3 = np.array(list(itertools.product(_x, _x, _y)))
_grid4 = np.array(list(itertools.product(_x, _x, _y, _z)))


#parameter labels
_dist_plabels ={'wb2' : ['lambda','kappa'],
                'wb3' : ['lambda','kappa','c'],
                'gp2' : ['alpha','beta'],
                'gp3' : ['alpha','beta','c'],
                'lg2' : ['mu','s'],
                'lg3' : ['mu','s','c'],
                'gl3' : ['alpha','beta','d'],
                'gl4' : ['alpha','beta','c','d'],
                'llg2': ['a','b'],
                'llg3': ['a','b','c']}

_dist_llabels ={'wb2' : [r'\lambda',r'\kappa'],
                'wb3' : [r'\lambda',r'\kappa','c'],
                'gp2' : [r'\alpha',r'\beta'],
                'gp3' : [r'\alpha',r'\beta','c'],
                'lg2' : [r'\mu','s'],
                'lg3' : [r'\mu','s','c'],
                'gl3' : [r'\alpha',r'\beta','d'],
                'gl4' : [r'\alpha',r'\beta','c','d'],
                'llg2': ['a','b'],
                'llg3': ['a','b','c']}





class survt_df(object):
    '''
    survival object:
    ----------------
    survt_df(parental_df, time_cols=['Time0', 'Time1']))
    parental_df is a Pandas.Dataframe of time-to-event data, such as:


    treatment Time0 Time1
    eat_apple 10    11    #this is a case of interval censored data
    eat_ice   5     5     #this is a case of event-at-time of exactly 5 (time unit)
    eat_dirt  10    +inf  #this is a right censored data (knowning only the left bound)
    eat_dust  0     17    #this is a left censored data (kowning only the right bound)

    Default time column names are 'Time0' and 'Time1'.

    Attributes:
    -----------
    emp_csdh:
    pandas.dataframe
    emperical CDF (cumulative distribution fucnction), SF (survival function) and HF (hazard)
 
    evt_grouped:
    pandas.dataframe
    time-to-event data grouped into each levels of treatment conbimations.

    mid_time:
    pandas.dataframe
    if the origial data is interval censored, mid_time contains (Time0+Time1)*0.5
    (also known as middle time estimates)
  
    treat_cols:
    list of strings
    The column names of treatments
    
    trt_levels
    list of tuples
    The combinations of treatment levels.    
    '''

    def __init__(self, parental_df, time_cols=['Time0', 'Time1']):
        self.treat_cols = [item for item in parental_df.columns if item not in time_cols]
        
        #Groupby treatment levels. If no treatment, creat a single level.
        if not self.treat_cols:
            p_df             = parental_df.copy()
            p_df['Group']    = 0
            self.trt_levels  = [0,]
            self.treat_cols  = ['Group']
            self.evt_grouped = p_df.groupby(self.treat_cols).agg(lambda x: tuple(x))
        else:
            self.evt_grouped = parental_df.groupby(self.treat_cols).agg(lambda x: tuple(x))
        self.mid_time = pd.DataFrame(index=self.evt_grouped.index)

        #every thing are converted to interval data format
        if self.evt_grouped.shape[1] == 2: #not interval data
            self.evt_grouped.columns = ['Time0', 'Time1']
            self.mid_time['Time']    = [tuple(item.mean(0)) for item in 
                                        map(np.array, zip(self.evt_grouped.Time0, 
                                                          self.evt_grouped.Time1))]
        elif self.evt_grouped.shape[1] == 1:
            self.evt_grouped.columns = ['Time0']
            self.evt_grouped['Time1']= self.evt_grouped['Time0']
            self.mid_time['Time']    = self.evt_grouped['Time0']

        #a dictionary contains the emprical CDF, density, hazard, etc. NAN discarded.                   
        self.mid_time['Time'] = [filter(np.isfinite, item) for item in self.mid_time['Time']]
        self.emp_csdh    = {}
        for item in self.evt_grouped.index.values:
            self.emp_csdh[item] = self._emp_cdf(self.mid_time.loc[item]['Time']).reset_index(drop=True)
        if len(self.treat_cols)==1:
            self.trt_levels = self.evt_grouped.index.tolist()
        else:
            self.trt_levels  = map(list, self.evt_grouped.index.levels)
        self._all_evt    = (np.hstack(self.evt_grouped.Time0),
                            np.hstack(self.evt_grouped.Time1))
        self._all_med    = np.nanmedian(self._all_evt)

            

    def _emp_cdf(self, data):
        '''Emprical Cdf, Sf and Hzd functions
        '''
        df        = pd.Series(data).value_counts(normalize=True).reset_index().sort('index')
        df.columns= ['time', 'pdf']
        df['cdf'] = np.cumsum(df['pdf'])
        df['sf']  = 1-df.cdf
        df['hzd'] = df.pdf/df.sf/np.hstack((df.time.diff().values[1:],
                                            np.nan))   
        return df




class distfit_df(object):
    '''
    Distribution fit object:
    ----------------
    distfit_df(survt_df, distribution)
    survt_df is a survival object.
    distribution is one of the followings:
    'wb2' : 2 parameter Weibull
    'wb3' : 3 paremeter Weibull (with Makeham term)
    'gp2' : 2 parameter Gompertz
    'gp3' : 3 parameter Gompertz (Gompertz-Makeham)
    'gl3' : 3 parameter Gompertz-logistic
    'gl4' : 4 parameter Gompertz-logistic-Makham
    'lg2' : 2 parameter logistic
    'lg3' : 3 parameter logistic (with Makeham term)
    'llg2': 2 parameter log-logistic
    'llg3': 3 parameter log-logistic (with Makeham term)


    Attributes:
    ----------------------

    distribution:
    string
    the name of the distribution model

    mdl_all_cnst:
    numpy.array
    model fit to data, with all data pooled and treatments ignored

    mdl_all_free
    pandas.dataframe
    model fit to data of each treatment level independently

    var
    pandas.dataframe
    variance for parameter estimates in mdl_all_free

    vcov_m
    numpy.ndarray
    variance covariance matrix for parameter estimates in mdl_all_free


    Methods:
    --------
    contourf_2D:
    Plotting filled contour plots of multi-dimensional confidence regions.

    logm:
    log-mean lifespan function

    logm_nc
    log-mean lifespan function, with Makeham term ignored (if the model contains Makeham term)    

    logq:
    log-quantile function

    logq_nc:
    log-quantile function, with Makeham term ignored (if the model contains Makeham term)    

    logt_var
    variance for log survival time

    logt_var_nc
    variance for log survival time, with Makeham term ignored (if the model contains Makeham term)

    plot
    line plot of model fit, with dot plot of observed data
    '''

    def __init__(self, s_df, distribution):
        if isinstance(s_df, survt_df):
            pass
        else:
            pass #TODO complain bad input
        if distribution in ['wb2','wb3','gp2','gp3','gl3','gl4',
                            'lg2','lg3','llg2','llg3']:
            pass
        else:
            pass #TODO complain not implemented
        
        self._mid_time     = s_df.mid_time
        self._emp_csdh     = s_df.emp_csdh
        self._evt_grouped  = s_df.evt_grouped
        self._all_med      = s_df._all_med
        self.distribution  = distribution
        self._dist         = eval('mortdist.%sllk_lp'%distribution) 
        self._cpsh         = eval('mortdist_cpsh.%slg_cpsh'%distribution)
        if distribution in ['wb3', 'gp3', 'lg3', 'llg3', 'gl4']:
            if '4' in distribution:
                ref_distr  = distribution.replace('4', '3')
                self._vidx = [0, 1, 3]
            else:
                ref_distr  = distribution.replace('3', '2')
                self._vidx = [0, 1,]
            self._dist_nC  = eval('mortdist.%sllk_lp'%ref_distr)
            self._cpsh_nC  = eval('mortdist_cpsh.%slg_cpsh'%ref_distr)
            self._labl_nC  = _dist_plabels[ref_distr]
            self._haveC    = True
        self._grid         = eval('_grid%s'%distribution[-1])
        self.mdl_all_free  = map(lambda x, y: self._fit_dist(distribution, 
                                                             x, 
                                                             y, 
                                                             step=s_df._all_med*0.001), 
                                 s_df.evt_grouped.Time0,
                                 s_df.evt_grouped.Time1)
        self.mdl_all_free  = pd.DataFrame(np.vstack([map(np.hstack, self.mdl_all_free)]))
        self.mdl_all_cnst  = np.hstack(self._fit_dist(distribution,
                                                      s_df._all_evt[0],
                                                      s_df._all_evt[1], 
                                                      step=s_df._all_med*0.001))
        self.mdl_all_free.index = s_df.evt_grouped.index
        self.mdl_all_free.columns =_dist_plabels[distribution] + ['logL', 'optimizer_flag']
        self._p_free = self.mdl_all_free
        
        self._p_cnst = self.mdl_all_free+np.nan
        self._p_cnst.ix[0] = self.mdl_all_cnst
        self._p_cnst.fillna(method='pad', inplace=True)
        self._p_cnst['logL'] = [self._dist(self.mdl_all_cnst[:int(distribution[-1])],
                                           np.vstack((s_df.evt_grouped.ix[i, 'Time0'], 
                                                      s_df.evt_grouped.ix[i, 'Time1']))) 
                                for i in range(len(s_df.evt_grouped))]
        self._no_pos = self.mdl_all_free.shape[1] - 2
        self._step    = s_df._all_med*0.001
        self.vcov_m, self.var = self._var()
        
        

    def _fit_dist_grid(self, data0, data1, optimizer=so.fmin, p_guess=None):
        '''Eval llk function on the search grid, find the optimal, and refine it
        using the optimizer supplied.
        If p_guess is supplied as initial guess parameter, the grid search step is skipped.
        '''
        data = np.vstack([np.array(data0), np.array(data1)])
        if p_guess is None:
            x       = mortdist.scan_m(data, self._grid, self._dist)
            idx     = np.nanargmin(x)
            p_guess = self._grid[idx]
        R = optimizer(self._dist, p_guess, args=(data,), 
                      full_output=1, disp=0)
        p      = R[0]
        v      = R[1]
        status = R[-1]
        return [p, v, status] #maybe we should get the result status code too?

    
        
    def _fit_dist(self, dist_name, data0, data1=None, step = None, optimizer=so.fmin):
        '''Turn the data to interval censored and calculate the fit.
        Refine the fit with the original data.
        (Interval-censored data is easy to fit, due to not having to calculate PDF)
        '''
        if data1 is None:
            data1=data0
        if step is None:
            step = 1.0 #maybe 1.0 is better? #TODO
        data_c = map(lambda x: x+step, data1)
        R = self._fit_dist_grid(data0, data_c, optimizer=so.fmin)
        #make an approximate interval-censor fit
        #for interval-censored data, make a fit with a modified right-limit
        #such that mixed-censored pointes will be handled correctly
        R = self._fit_dist_grid(data0, data1, optimizer=so.fmin, p_guess=R[0]) #make a final fit
        return R



    def _var(self):
        '''Calculate variance-covariance matrix as the negative inverse of Hessian
        '''
        vcov_m = []
        for i in range(self._evt_grouped.shape[0]):
            #raise an uninverseable matrix error? 
            vcov = mortdist.einverse(mortdist.derivative2(self.mdl_all_free.irow(i)[self.mdl_all_free.columns[:-2]],
                                                         np.vstack(self._evt_grouped.irow(i)),
                                                         self._dist))
            vcov_m.append(vcov)
        v_df = pd.DataFrame(np.vstack(map(np.diag, vcov_m)),
                            index=self._evt_grouped.index, 
                            columns=self.mdl_all_free.columns[:-2])
        v_df = np.sqrt(v_df)
        return dict(zip(self._evt_grouped.index.tolist(), vcov_m)), v_df



    def logq(self, key, q=0.5):
        '''log-quantile function.
        key: which treatment combination to be calculated?
        q:   quantile, Default to 0.5 (calculate the median)
        '''
        p = self.mdl_all_free.ix[key][:-2].values
        tempfuc = lambda x: (self._cpsh(p, x)[2] - q)**2
        return np.log(so.fmin(tempfuc, 0, disp=0)[0])



    def logq_nc(self, key, q=0.5):
        '''log-quantile function. Makeham term ignored.
        key: which treatment combination to be calculated?
        q:   quantile, Default to 0.5 (calculate the median)
        '''
        if self._haveC:
            p = self.mdl_all_free.ix[key, self._labl_nC].values
            tempfuc = lambda x: (self._cpsh_nC(p, x)[2] - q)**2
            return np.log(so.fmin(tempfuc, 0, disp=0)[0])
        else:
            return self.logq(key, q)



    def logm(self, key, up=None):
        '''log-mean function.
        key: which treatment combination to be calculated?
        up:  right bound for numerical intergration. Default: None. (Inferred from data)
        '''
        if not up:
            up = self._all_med*10
        p = self.mdl_all_free.ix[key][:-2].values
        return np.log(si.quad(lambda x: self._cpsh(p, x)[1]*x, 0, up)[0])



    def logm_nc(self, key, up=None):
        '''log-mean function. Makeham term ignored
        key: which treatment combination to be calculated?
        up:  right bound for numerical intergration. Default: None. (Inferred from data)
        '''
        if not up:
            up = self._all_med*10
        if self._haveC:
            p = self.mdl_all_free.ix[key, self._labl_nC].values
            return np.log(si.quad(lambda x: self._cpsh_nC(p, x)[1]*x, 0, up)[0])
        else:
            return self.logm(key, up)



    def logt_var(self, key, t):
        '''Varince of time estimates
        key: which treatment combination to be calculated?
        up:  time value, in logscale, at which variance is to be calculated

        Example: if .logm(key) returns 0.5, which means the median lifespan 
           estimate is exp(0.5), the variance of median can then be calculated as
           .logt_var(key, 0.5)
        '''
        p = self.mdl_all_free.ix[key][:-2].values
        #remeber the loglikelihood function takes two dimension x 
        mpartial = np.matrix(mortdist.derivative2(p, [t,t], self._dist)) 
        return np.sum(np.array(mpartial.T * mpartial) * self.vcov_m[key])



    def logt_var_nc(self, key, t):
        '''Varince of time estimates with Makeham term ignored
        key: which treatment combination to be calculated?
        up:  time value, in logscale, at which variance is to be calculated

        Example: if .logm(key) returns 0.5, which means the median lifespan 
           estimate is exp(0.5), the variance of median can then be calculated as
           .logt_var(key, 0.5)
        '''
        if self._haveC:
            p        = self.mdl_all_free.ix[key, self._labl_nC].values
            mpartial = np.matrix(mortdist.derivative2(p, [t,t], self._dist_nC)) 
            vcov     = self.vcov_m[key][self._vidx, :][:, self._vidx]
            return np.sum(np.array(mpartial.T * mpartial) * vcov)
        else:
            return self.logt_var(key, t)        



    def plot(self, key, x_linspace, ax, kind='sf'):
        '''Plotting method.
        Model fit in line plot. Actual data in dot plot.
        key:        which treatment combination to be calculated?
        x_linspace: time points at which the fit will be evaluated.
        ax        : a matplotlib.axes, at which the plot will be located.
        kind      : one of 'cdf', 'pdf', 'sf', 'hzd'. Kind of plot        
        '''
        cpsh = dict(zip(['cdf', 'pdf', 'sf', 'hzd'],
                        self._cpsh(self.mdl_all_free.ix[key, _dist_plabels[self.distribution]], 
                                   x_linspace)))
        if kind=='hzd':
            ax.plot(self._emp_csdh[key].time,
                    np.log(self._emp_csdh[key].__getattr__(kind)), 
                    '+')        
            ax.plot(x_linspace, np.log(cpsh[kind]))

        else:
            ax.plot(self._emp_csdh[key].time,
                    self._emp_csdh[key].__getattr__(kind), 
                    '+')        
            ax.plot(x_linspace, cpsh[kind])



    def contourf_2D(self, key, ax, cax, tuple_prange, levels=None, cmap='jet_r', **kwargs):
        '''Filled contour plot for multidimensional confident region.
        key       : which treatment combination to be calculated?
        ax        : a matplotlib.axes, at which the plot will be located.
        cax       : a matplotlib.axes, at which the colorbar will be placed
        tuple_prange: a tuple of the parameter values at which log-likelihood function will be evalulated:
        levels    : colorbar levels
        cmap      : color map
        kwargs    : additional kwargs for matplotlib.axes.contourf()

        tuple_prange example:
        for a 3 parameter distribution:

                      (np.linspace(-5, 5, 41),
                       np.linspace(0,0,1),
                       np.linspace(-4, 4, 45)))

        the tuple should be 3 element long
        assuming the MLE is [p1, p2, p3], the log-likelihood function will be evalulated at:
        np.meshgrid(p1+np.linspace(-5, 5, 41), #41 points around p1
                    p2+np.linspace(0,0,1),     #at p2 only
                    p3+np.linspace(-4, 4, 45)) #45 points around p3
        '''

        if (np.array(map(len, tuple_prange))>1).sum()==2:
            pass
        else:
            raise('Error: 2 and only 2 so-and-so')
    
        if levels is None:
            levels = self.mdl_all_free.ix[key, self._no_pos] + \
                     np.hstack((0., ss.chi2.ppf(1-0.1**(np.arange(1,11)),2)*0.5))
            cax2 = cax.twinx()
            cax2.set_ylim((0,len(levels)-1))
            cax2.yaxis.set_label_text('$p = 10^{-n}$')
            cax2.yaxis.set_label_position('left')
            cax2.yaxis.tick_left()

    
        p_grid = []
        for i in range(self._no_pos):
            p_grid.append(self.mdl_all_free.ix[key, i] + tuple_prange[i])    
        llk_rst = mortdist.scan_m(np.vstack(self._evt_grouped.ix[key]),
                                  np.array(np.meshgrid(*p_grid)).T.reshape((-1, self._no_pos)),
                                  self._dist)

        ((ix, X) , (iy, Y)) = tuple([(idx, item) for idx, item in enumerate(p_grid) if item.size>1])
        if (ix, iy) == (0,1):
            llk_rst = llk_rst.reshape(len(X), len(Y)).T
        else:
            llk_rst = llk_rst.reshape(len(Y), len(X))
    
        CF = ax.contourf(X, Y, llk_rst, levels=levels, cmap=cmap, **kwargs)
        CB = plt.colorbar(CF, cax=cax)

        cax.yaxis.set_label_text('-log likelihood')
        cax.yaxis.set_label_position('right')
        cax.yaxis.tick_right()

        ax.set_xlabel('$log(%s)$'%_dist_llabels[self.distribution][ix])
        ax.set_ylabel('$log(%s)$'%_dist_llabels[self.distribution][iy])
#        return llk_rst




class model(object):
    '''Fit GLM to a distribution fit object. (Experimental)
    Methods:
    -------
    fit_partial_model:
    Fit a GLM, zero, one or more parameter can be constrained to a single scalar.
    '''
    def __init__(self, dst_df):
        self.distribution = dst_df.distribution
        self._evt_grouped = dst_df._evt_grouped
        self._dist        = dst_df._dist
        self._no_pos      = dst_df._no_pos
        self._p_free      = dst_df._p_free.icol(range(dst_df._no_pos)).values
        self._p_cnst      = dst_df._p_cnst.icol(range(dst_df._no_pos)).values
        self._llk_free    = dst_df._p_free.logL.sum()
        self._llk_cnst    = dst_df._p_cnst.logL.sum()
        self._vstack_evt  = map(np.vstack, zip(map(list, self._evt_grouped.Time0), 
                                               map(list, self._evt_grouped.Time1)))
        self._vstack_evtL = map(np.vstack, zip(map(list, self._evt_grouped.Time0), 
                                               [map(lambda x: x+dst_df._step, item) 
                                                for item in self._evt_grouped.Time0]))
        self.model        = None
        self.dmatrix      = None


                
    def fit_partial_model(self, model, fix_pos):
        '''Fit a GLM, zero, one or more parameter can be constrained to a single scalar.
        model  : a string of model terms
        fix_pos: parameters to be constrained to a single scalar.

        examples:
        some_weibull_2p_model.fit_partial_model('Diet*Sex', [1])
        weibull model, with the 2nd parameter fixed

        some_weibull_2p_model.fit_partial_model('Diet*Sex', [])
        weibull model, with the none of the parameters fixed
           and the same model formula will be applied to both weibull parameters
        '''
        dmatrix = patsy.dmatrices('Time0 ~ '+ model, self._evt_grouped.reset_index())
        dm = dmatrix[1]
        #self.model         = None
        #self.design_matrix = None
        #self.coefficient   = None
        return self.ceoff_model(dm, fix_pos)
    
    

    def dm_upk_pfixp(self, flatpp, fix_pos):
        '''unpack partial fixed parameter 1D array to nD array
        accroding to the fixed postion and total number of position specified,
        for use in conjunction with design matrix only: As only the 1st cell
        has the corresponding value and the the other cells of that position (fix_pos)
        are 0
        '''
        result = self.unpack_pfixp(flatpp, fix_pos)
        for i in fix_pos:
            result[range(1, len(result)), i] = 0.
        return result


    
    def model_l_pf(self, p, dma, fix_p=None):
        '''for partial-fixed models, 
           (if fix_p is empty works for full model)
        '''
        if fix_p is not None:
            pL = np.dot(dma, self.dm_upk_pfixp(p, fix_p))
        else:
            pL = np.dot(dma, p.reshape(-1, len(self._vstack_evt)).T)
        LL = map(self._dist, pL, self._vstack_evt)
        return sum(LL)


    
    def fix_p_matrix_init(self, fix_pos):
        '''Creat initial guess for partially fixed model parameters
        '''
        result = self._p_free*1.
        if isinstance(fix_pos, int):
            fix_pos=[fix_pos]
        for item in fix_pos:
            result[:, item] = np.zeros(result.shape[0]) + self._p_cnst[:,item]
        return result



    def pack_pfixp(self, arr_para, fix_pos):
        '''Pack an array of parameters (nD) to 1D, for fixed position, a single value is inserted
        '''
        fix_V = arr_para[0, fix_pos]
        _temp = arr_para*1.
        _temp[:, fix_pos] = np.nan
        _temp[0, fix_pos] = fix_V
        _temp = _temp.T.ravel()
        return _temp[np.isfinite(_temp)]

    

    def unpack_pfixp(self, flat_pp, fix_pos):
        '''The reverse of pack_pfixp
        '''
        if isinstance(fix_pos, int):
            fix_pos = [fix_pos]
        dim = (len(flat_pp)-len(fix_pos))/(self._no_pos-len(fix_pos))
        result = []
        for i in range(self._no_pos):
            for j in range(dim):
                if i in fix_pos:
                    result.append(i*dim)
                else:
                    result.append(i*dim+j)
        result = map(dict(zip(set(result), range(len(result)))).get, result)
        result = [flat_pp[idx] for idx in result] #comment this line out to return index only 
        return np.array(result).reshape(self._no_pos, -1).T



    def cst_rpf_packed(self, dm, p_ndarray, fix_pos):
        '''Based on design matrix, parameter estimates, and fixed position. Return the
        initial guess parameter for optimization.
        '''
        if isinstance(fix_pos, int):
            fix_pos = [fix_pos,]
        result_p = []
        f = lambda x, i: ((np.dot(dm, x) - p_ndarray[:,i])**2).sum() 
        if np.diff(dm.shape)==0:
            for i in range(len(p_ndarray[0])):
                if i in fix_pos:
                    result_p.append(p_ndarray[0][i])
                else:
                    result_p.append(np.linalg.solve(dm, p_ndarray[:,i]))            
        else:
            for i in range(len(p_ndarray[0])):
                if i in fix_pos:
                    result_p.append(p_ndarray[0][i])
                else:
                    result_p.append(so.fmin_bfgs(f, [0,]*(dm.shape[1]),
                                                 args = (i,), disp=0))
        return np.hstack(result_p)



    def fit_pfix(self, fix_pos, optimizer=so.fmin_slsqp):
        '''Making the actual fit
        '''
        def pf_llk(p, data):
            '''data: either self._vstack_evt or self._vstack_evtL
            '''
            return sum(map(self._dist, 
                           self.unpack_pfixp(p, fix_pos), 
                           data))

        if (fix_pos is None) or (len(fix_pos)==0):
            #----No fixed position, return a full model----
            return (self._p_free,
                    self._llk_free,
                    self._p_free.ravel(), 
                    self._p_free.size, 0)
        elif len(fix_pos) == self._no_pos:
            #----All position fixed, return a fully-fixed model----
            return (self._p_cnst,
                    self._llk_cnst,
                    self._p_cnst[0], 
                    self._p_cnst.shape[1], 0)
        else:
            #----Some position fixed, some not----
            #====compare linear fit and full-fixed fit and choose the better initial guess=====
            pd_p  = self.pack_pfixp(self.fix_p_matrix_init(fix_pos), 
                                    fix_pos)
            L0    = pf_llk(pd_p, self._vstack_evt)
            L1    = self._llk_cnst#; print 'dist_name,fix_pos,L0,L1' #L1: Likelihood of the fully-fixed model
            if (np.isfinite(L0)&(L0<L1)):
                pass
            else:
                pd_p = self.pack_pfixp(self._p_cnst, fix_pos)
                L0   = L1    #Therefore: L0 always store the better fit
            #print pd_p, pf_llk(pd_p, xL0, xL1, fix_pos, no_pos, dist_name)

            #=====Optimizer=====
            _args = (self._vstack_evtL,)
            F     = lambda p,x: pf_llk(p, *_args)
            dF    = lambda p,x: mortdist.derivative(p,x,F)
            #####=====bound is determined by distribution, for SLSQP=====#####
            if self.distribution in ['gl3', 'gp3', 'gp2', 'gl4']:
                _bnds = [(-41.23, -0.257),]*len(pd_p)
            else:
                _bnd_lst = [(-5,8),(-5,8),(-41.3,1)]
                _bnds    = [[_bnd_lst[i],] if i in fix_pos else [_bnd_lst[i],]*len(self._p_cnst) 
                            for i in range(self._no_pos)]
                _bnds    = list(itertools.chain(*_bnds))
            #####=====If initial guess is good use SLSQP+SMPLEX, else SMPLEX only=====#####
            if np.isfinite(pf_llk(pd_p, *_args)):
                #Do a fit on self._vstack_evtL?
                #R     = optimizer(F, pd_p, disp=0, args=([0,0],), full_output=1, fprime=dF, bounds=_bnds); print '1st',R[1]
                _args = (self._vstack_evt,)
                if self.distribution in []: #['gl3', 'gp3', 'gp2', 'gl4']: #if run the self._vstack_evtL OPZ step, change pd_p to R[0]
                    R = optimizer(F, pd_p, disp=0, args=([0,0],), full_output=1, fprime=dF, bounds=_bnds)#; print 'SLSQP ',R[1]
                else:
    #                R = so.fmin_powell(F, pd_p, disp=0, args=([0,0],), full_output=1); print 'POWELL',R[1]
    #            if np.isfinite(R[1]) & (R[1]<=L0):
    #                pass
    #            else: #Indicating a difficult problem: use last resort.
                    R = so.fmin(F, pd_p, disp=0, args=([0,0],), full_output=1)#; print 'SMPLEX',R[1] 
            else:
                _args = (self._vstack_evt,)
                R     = so.fmin(F, pd_p, disp=0, args=([0,0],), full_output=1)#; #print '2nd',R[1]            
            return (self.unpack_pfixp(R[0], fix_pos), 
                    float(R[1]), 
                    R[0], 
                    len(R[0]), 
                    R[-1])
    #on MacOS-64bit, the 2nd round of SLSQP optimization often fail, for models other than Gompertz class.
    #not appears to be a problem in WIN64 or WIN32.



    def ceoff_model(self, dm, fix_pos):
        '''Calculate the coefficients of the linear predictors, and their SE. 
        '''
        columns  = dm.design_info.column_names
        model    = dm.design_info.describe()
        solution = self.cst_rpf_packed(dm, self.fit_pfix(fix_pos)[0], fix_pos)
        #----at this stage, `solution` is in a packed form----
        if (np.diff(dm.shape) == 0)&(set(np.unique(np.asarray(dm))) == {0., 1.}):
            #----when it is a complete model, skip the additional optimizer run----
            LL = self.model_l_pf(solution, dm, fix_pos); print 'full model'
        else: 
            #----when it is not a full model with all possible interactions, run optimizer----
            solution, LL = so.fmin(self.model_l_pf, 
                                   solution,
                                   args = (dm, fix_pos), 
                                   disp=0, full_output=1)[:2]
        p_label  = _dist_plabels[self.distribution]
        m_index  = []
        for i, v in enumerate(p_label[:self._no_pos]):
            if i in fix_pos:
                m_index.append((v, 'FIXED'))
            else:
                for v1 in columns:
                    m_index.append((v, v1))

        #----at this stage, the solution must be in its packed from----
        df_result             = pd.DataFrame({'Estimate': solution})
        try: #----I suspect this error catching is redunant---- TODO
            df_result.index       = pd.MultiIndex.from_tuples(m_index)
        except:
            print m_index

        if (len(fix_pos)==0) and (dm.shape[1]==dm.shape[0]):
            #----when it is a complete model, simply calculate SE----
            varr = np.vstack([np.diag(mortdist.inverse(mortdist.derivative2(p, x, self._dist))) 
                              for p, x in zip(self._p_free, self._vstack_evt)])
            df_result['Error'] = np.sqrt(np.dot(np.abs(np.linalg.inv(np.array(dm))).astype(int),
                                                varr)).T.flatten()
            #print 'full model'
        else:
            #----when it is not a full model with all possible interactions, run optimizer and get the Hessian----
            f = lambda p, _x: self.model_l_pf(p, dm, fix_pos)
            H = mortdist.derivative2(solution, [0,0], f)
            try:
                df_result['Error']= np.sqrt(np.diag(mortdist.inverse(H)))
            except np.linalg.LinAlgError:
                df_result['Error']= np.nan 
                #effective only when using np.linalg.inv instead of mortdist.inverse
                #in the try block.
                #the reason: in OSX, such an error may kill ipython kernel. Ugly hack!
                print "Singular Hessian Matrix, VCOV can't be evaluated."

        df_result['Std_Estm'] = df_result.Estimate/df_result.Error
        df_result['p']        = ss.norm.sf(df_result['Std_Estm'].abs())*2
        df_result['sig_code'] = map(util.prob_to_code, df_result.p)
        print LL
        print model
        return df_result
