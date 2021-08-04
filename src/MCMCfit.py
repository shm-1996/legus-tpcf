from header import *
from Galaxy import Galaxy
from Plot_Class import myPlot
import emcee
import corner
from TPCF import bootstrap_two_point,bootstrap_two_point_angular,linear_function,linear_truncation,\
smooth_function,onepowerlaw_function,piecewise_truncation

################################################################################################
#Piecewise power law functions
def model(params,data):
    #Parameters to fit
    A_1,alpha_1,alpha_2,beta = params
    #A_2 in terms of A_1 assuming continuity
    A_2 = A_1
    function = np.piecewise(data,[np.log(data)<beta],[lambda data :  A_1 + alpha_1 * np.log(data), 
        lambda data : A_2 + (alpha_1-alpha_2)*beta + alpha_2*np.log(data)])
    return function

def lnprior(params,distance,bins):
#    A_1,A_2,alpha_1,alpha_2,beta = params
    A_1,alpha_1,alpha_2,beta = params
    # beta_limits = [500.0, 2000.0]
    # beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
    # beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
    beta_limits = [np.min(bins),np.max(bins)]
    
    if(-10<A_1<10 and 
       -5<alpha_1<0 and alpha_1<alpha_2<0 and
      np.log(beta_limits[0])<beta<np.log(beta_limits[1])) :
        return 0
    else :
        return -np.inf
def lnprob(params, bins,data,data_error,distance):
    lp = lnprior(params,distance,bins)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(params, bins,data,data_error) 

def lnlike(params,bins,data,data_error):
#    data_model = model_1(params,bins)
    data_model = model(params,bins)
    log_likelihood = -1/2. * np.sum(((data-data_model)/data_error)**2)
    return log_likelihood    

################################################################################################
#Single power law functions

#1. Model
def model_singlepower(params,data):
    A_1,alpha_1 = params
    function = A_1 + np.log(data)*alpha_1
    return function

#2. Priors
def lnprior_singlepower(params):
#    A_1,A_2,alpha_1,alpha_2,beta = params
    A_1,alpha_1 = params
    if(-10<A_1<10 and 
       -5<alpha_1<0 ) :
        return 0
    else :
        return -np.inf

#3. lnprob
def lnprob_singlepower(params, bins,data,data_error):
    lp = lnprior_singlepower(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_singlepower(params, bins,data,data_error)

#4. log-likelihoods
def lnlike_singlepower(params,bins,data,data_error):
#    data_model = model_1(params,bins)
    data_model = model_singlepower(params,bins)
    log_likelihood = -1/2. * np.sum(((data-data_model)/data_error)**2)
    return log_likelihood

###############################
# Power-law with exponential truncation
def model_linear_truncation(params,data) :
    A1,alpha_1,theta_c = params 
    function = A1 + alpha_1*np.log(data) - data/theta_c
    return function

#2. Priors
def lnprior_lineartruncation(params,bins):
    A_1,alpha_1,theta_c = params
    thetac_limits = [np.min(bins),np.max(bins)*5.0]
    beta_limits = [np.min(bins),np.max(bins)]
    if(-10<A_1<10 and 
       -5<alpha_1<0 and thetac_limits[0]<theta_c<thetac_limits[1])  :
        return 0
    else :
        return -np.inf

#3. lnprob
def lnprob_lineartruncation(params, bins,data,data_error):
    lp = lnprior_lineartruncation(params,bins)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_lineartruncation(params, bins,data,data_error)

#4. log-likelihoods
def lnlike_lineartruncation(params,bins,data,data_error):
#    data_model = model_1(params,bins)
    data_model = model_linear_truncation(params,bins)
    log_likelihood = -1/2. * np.sum(((data-data_model)/data_error)**2)
    return log_likelihood
################################################################################################
################################################################################################
#Fitting routines

def fit_MCMC_galaxy(galaxy_name,method='masked_radial',function='piecewise',omega1=False,age=None):
    """
    Fit an MCMC to a galaxy with given name. 
    Parameters: 
        galaxy_name : str 
            Name of galaxy for which to fit
        method : str
            Method of Random distribution
        function: str 
            Piecewise or single PL fits to the TPCF


    """
    galaxy_name = galaxy_name.upper()
    galaxy_class = Galaxy(galaxy_name,verbose=False)
    #Save in subdirectories for below two methods
    if(method == 'masked'):
        galaxy_class.outdir += '/Masked/'
    elif(method == 'uniform'):
        galaxy_class.outdir += '/Uniform/'
    elif(method == 'masked_radial'):
        galaxy_class.outdir += '/Masked_Radial/'
    elif(method == 'exponential'):
        galaxy_class.outdir += '/Exponential/'
    else:
        raise myError("Method not recognised.")

    print("Performing MCMC for galaxy {}".format(galaxy_name))

    galaxy_class = loadObj(galaxy_class.outdir+
            '{}_summary'.format(galaxy_class.name))

    if(function == 'piecewise'):
        print("Fitting Piecewise PL with MCMC.")
        fit_MCMC(galaxy_class,save=True,function='piecewise',omega1=omega1,age=age)
    elif(function == 'singlepl') :
        print("Fitting Single PL fit in MCMC.")
        fit_MCMC(galaxy_class,save=True,function='singlepl',omega1=omega1,age=age)
    elif(function == 'singletrunc'):
        print("Fitting Single PL with exponential truncation fit in MCMC.")
        fit_MCMC(galaxy_class,save=True,function='singletrunc',omega1=omega1,age=age)
    elif(function == 'doubletrunc'):
        print("Fitting Piecewise PL with exponential truncation fit in MCMC.")
        fit_MCMC(galaxy_class,save=True,function='doubletrunc',omega1=omega1,age=age)
    

def fit_MCMC(galaxy_class,save=True,function='piecewise',omega1=True,
    age=None,plot=False,debug=False):
    """
    Fit an MCMC to the TPCF of a galaxy and create diagnostic plots after. 

    Parameters:
        galaxy_class: Class Galaxy
            Galaxy class object 
        save : Boolean
            Flag to save the figure, else just plot. 
        function: string 
            Kind of function to fit to the data 
    """
    from Galaxy import myError

    #Safety check
    if(function not in ['singlepl','piecewise','singletrunc','doubletrunc']):
        raise ValueError("Function should be either singlepl, piecewise or singletrunc")

    #Set directory where results would be saved
    MCMC_directory = os.path.abspath(galaxy_class.outdir)
    if(function == 'singlepl'):
        funcTag = "singlePL"
    elif(function == 'piecewise'):
        funcTag = "pieceWise"
    elif(function == 'singletrunc'):
        funcTag = 'singleTrunc'
    elif(function == 'doubletrunc'):
        funcTag = 'doubleTrunc'

    #Set initial guess
    if(function == 'piecewise'):
        bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        beta_limits = [np.min(bins),np.max(bins)]
        galaxy_class.fit_TPCF(function='piecewise',method='single',omega1=omega1)
        initial_guess = galaxy_class.fit_values
        #Modify initial guess to be in b/w 50 & 300 pc
        distance = galaxy_class.distance*const.Parsec*1.e6
        # # beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
        # # beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
        # initial_guess[3] = np.log((beta_limits[0]+beta_limits[1])/2.)

    elif(function == 'singlepl'):
        initial_guess = 5.0,-1.0
    elif(function == 'singletrunc'):
        galaxy_class.fit_TPCF(function='singletrunc',method='single',omega1=omega1)
        initial_guess = galaxy_class.fit_values[0],galaxy_class.fit_values[1],\
galaxy_class.fit_values[2]

    #Properties of MCMC
    nwalkers = 500
    ndim = len(initial_guess)
    step = 1.e-7

    p0 = [np.array(initial_guess) + step * np.random.randn(ndim) for i in range(nwalkers)]

    if(age == 'young'):
        corr = galaxy_class.ycorr
        dcorr = galaxy_class.ydcorr
    elif(age == 'old'):
        corr = galaxy_class.ocorr
        dcorr = galaxy_class.odcorr
    else:
        corr = galaxy_class.corr
        dcorr = galaxy_class.dcorr 

    #Set bins/x-data: i.e. the independent variable to fit with
    separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)
    indices = np.where(np.logical_and(corr>-1.0,
            np.abs(corr)>np.abs(dcorr)))
    corr_fit = corr[indices].astype(np.float)
    dcorr_fit = dcorr[indices].astype(np.float)
    separation_bins = separation_bins[indices].astype(np.float)
    #Define MCMC samplers

    #Include distance in argument for piecewise model sampler
    if(function == 'piecewise'):
        if(omega1 is True):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                    args=(separation_bins,np.log(1+corr_fit),
                                         dcorr_fit/(1+corr_fit),distance))
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                    args=(separation_bins,np.log(corr_fit),
                                         dcorr_fit/corr_fit,distance))
    elif(function == 'singlepl'):

        if(omega1 is True):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_singlepower, 
                                    args=(separation_bins,np.log(1+corr_fit),
                                         dcorr_fit/(1+corr_fit)))
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_singlepower, 
                                    args=(separation_bins,np.log(corr_fit),
                                         dcorr_fit/corr_fit))
    elif(function == 'singletrunc'):
        
        if(omega1 is True):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_lineartruncation, 
                                    args=(separation_bins,np.log(1+corr_fit),
                                         dcorr_fit/(1+corr_fit)))
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_lineartruncation, 
                                    args=(separation_bins,np.log(corr_fit),
                                         dcorr_fit/corr_fit))
    elif(function == 'doubletrunc'):

        if(omega1 is True):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_piecewisetruncation, 
                                    args=(separation_bins,np.log(1+corr_fit),
                                         dcorr_fit/(1+corr_fit),distance))
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_piecewisetruncation, 
                                    args=(separation_bins,np.log(corr_fit),
                                         dcorr_fit/corr_fit,distance))

    #Run MCMC
    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, 10000,progress=True)

    print("MCMC operation completed....")

    #Get samples
    samples = sampler.flatchain

    if(function == 'singlepl'):
        if(age == 'young'):
            galaxy_class.fit_values_singlepl_y = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singlepl_y = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singlepl_y = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singlepl_y = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singlepl_y = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singlepl_y = np.percentile(samples,84,axis=0)
        elif(age == 'old'):
            galaxy_class.fit_values_singlepl_o = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singlepl_o = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singlepl_o = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singlepl_o = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singlepl_o = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singlepl_o = np.percentile(samples,84,axis=0)
        else:
            galaxy_class.fit_values_singlepl = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singlepl = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singlepl = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singlepl = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singlepl = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singlepl = np.percentile(samples,84,axis=0)
    elif(function == 'piecewise'):
        if(age == 'young'):
            galaxy_class.fit_values_piecewise_y = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_piecewise_y = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_piecewise_y = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_piecewise_y = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_piecewise_y = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_piecewise_y = np.percentile(samples,84,axis=0)

            #Store distribution of beta as it is used elsewhere
            galaxy_class.beta_distribution = samples[:,3]
        elif(age == 'old'):
            galaxy_class.fit_values_piecewise_o = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_piecewise_o = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_piecewise_o = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_piecewise_o = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_piecewise_o = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_piecewise_o = np.percentile(samples,84,axis=0)
        else:
            galaxy_class.fit_values_piecewise = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_piecewise = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_piecewise = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_piecewise = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_piecewise = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_piecewise = np.percentile(samples,84,axis=0)
        
    elif(function == 'singletrunc'):
        if(age == 'young'):
            galaxy_class.fit_values_singletrunc_y = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singletrunc_y = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singletrunc_y = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singletrunc_y = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singletrunc_y = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singletrunc_y = np.percentile(samples,84,axis=0)

            #Take theta distribution from young TPCF, since NGC 5253 has no old TPCF
            if(galaxy_class.name == 'NGC_5253'):
                galaxy_class.theta_distribution = samples[:,2]
        elif(age == 'old'):
            galaxy_class.fit_values_singletrunc_o = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singletrunc_o = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singletrunc_o = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singletrunc_o = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singletrunc_o = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singletrunc_o = np.percentile(samples,84,axis=0)

            # Theta distribution used to infer galaxy rc
            galaxy_class.theta_distribution = samples[:,2]
        else:
            galaxy_class.fit_values_singletrunc = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors_singletrunc = sampler.flatchain.std(axis=0)
            galaxy_class.maxlnprob_singletrunc = np.max(sampler.flatlnprobability)
            galaxy_class.median_fit_singletrunc = np.percentile(samples,50,axis=0)
            galaxy_class.low_error_singletrunc = np.percentile(samples,16,axis=0)
            galaxy_class.high_error_singletrunc = np.percentile(samples,84,axis=0)
    else:
        raise myError("Function type not recognised")

    if(plot):
        #Diagnostic plots
        print("Making diagnostic plots...")

        #Diagnostic Plot 1 : Corner Plot
        if(function == 'piecewise'):
            labels = [r'$A_1$',r'$\alpha_1$',r'$\alpha_2$',r'$\beta$']        
        elif(function == 'singlepl'):
            labels = [r'$A_1$',r'$\alpha_1$']
        elif(function == 'singletrunc'):
            labels = [r'$A_1$',r'$\alpha_1$',r'$\theta_c$']
        elif(function == 'doubletrunc'):
            labels = [r'$A_1$',r'$\alpha_1$',r'$\alpha_2$',r'$\beta$',r'$\theta_c$']

        fig = corner.corner(samples,show_titles=True,labels=labels,
            plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
        if(save):
            fig.savefig(MCMC_directory+'/{}_{}_Posterior_Dist.pdf'.format(galaxy_class.name,funcTag),bbox_inches='tight')
            plt.close()
        else:
            fig.show()

        #Diagnostic Plot 4: Best-Fit Version
        galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
        galaxy_class.fit_errors = sampler.flatchain.std(axis=0)
        pl = myPlot(galaxy_class)
           
        pl.plot_TPCF(save=save,filename=MCMC_directory+'/{}_{}_Fit.pdf'.format(galaxy_class.name,funcTag),
            function=function,omega1=omega1)

    if(debug):

    #Diagnostic Plot 2: Convergence Check
        chain = sampler.get_chain(discard=200,thin=40)
        fig, axes = plt.subplots(ndim, figsize=(8, 7), sharex=True)
        for i in range(ndim):
            ax = axes[i]
            ax.plot(chain[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(chain))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("step number")
        if(save):
            fig.savefig(MCMC_directory+'/{}_{}_Convergence.pdf'.format(galaxy_class.name,funcTag),bbox_inches='tight')
            plt.close()
        else:
            fig.show()

        #Diagnostic Plot 3: Overplot sample fits on data
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(sep_to_pc,pc_to_sep))

        separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
        separation_bins*=(1./arcsec_to_degree)
        chain = sampler.get_chain(discard=200,thin=40,flat=True)
        inds = np.random.randint(len(chain), size=50)
        
        for ind in inds:
            sample = chain[ind]
            
            if(function == 'piecewise'):
                axs.plot(separation_bins,np.exp(model(sample,separation_bins)),
                 alpha=0.1)
            elif(function == 'singlepl'):
                axs.plot(separation_bins,np.exp(model_singlepower(sample,separation_bins)),
                    alpha=0.1)
            elif(function == 'singletrunc'):
                axs.plot(separation_bins,np.exp(model_linear_truncation(sample,separation_bins)),
                 alpha=0.1)
            elif(function == 'doubletrunc'):
                axs.plot(separation_bins,np.exp(model_piecewise_truncation(sample,separation_bins)),
                 alpha=0.1)            

        if(omega1 is True):
            axs.errorbar(separation_bins,1+galaxy_class.corr,yerr=galaxy_class.dcorr,
                    fmt='.-')
        else:
            axs.errorbar(separation_bins,galaxy_class.corr,yerr=galaxy_class.dcorr,
                    fmt='.-')

        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.callbacks.connect("xlim_changed", axs_to_parsec)
        axs.legend()
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

        if(save):
            plt.savefig(MCMC_directory+'/{}_{}_Samplefits'.format(galaxy_class.name,funcTag),bbox_inches='tight')
            plt.close()
        else:
            plt.show()
        #Save sampler if in debug mode
        if(save):
            saveObj(sampler,MCMC_directory+'/{}_{}_MCMCsampler'.format(galaxy_class.name,funcTag))

    return



def sep_to_pc(sep) :

    """Angular separation in arcsec to parsec

    Parameters
    ----------
    sep : float
        separation in angular 
    Returns
    -------
    pc : float
         The separation in parsec for given angular separation
    """

    # Distance to NGC 628 = 9.9 Mpc : Olivares et al 2010
    distance_628 = 9.9*const.Parsec*1.e6  
    arcsec_to_pc = u.arcsec.to(u.radian)*distance_628/(const.Parsec)
    return sep*arcsec_to_pc

def pc_to_sep(pc) :
    """Linear separation in parsec to arcsec

    Parameters
    ----------
    pc : float
        separation in linear
    Returns
    -------
    pc : float
         The separation in parsec for given angular separation
    """
    # Distance to NGC 628 = 9.9 Mpc : Olivares et al 2010
    distance_628 = 9.9*const.Parsec*1.e6  
    pc_to_radian = pc*const.Parsec/distance_628
    return pc_to_radian*u.radian.to(u.arcsec)

def axs_to_parsec(axs) :
    """Draw twin axes corresponding to angular limits

    Parameters
    ----------
    axs : Matplotlib axis instance
          Matplotlib axis on which to plot
    Returns
    -------
    None
    """
    xmin,xmax = axs.get_xlim()
    ax2.set_ylim(sep_to_pc(xmin),sep_to_pc(xmax))
    ax2.figure.canvas.draw()


