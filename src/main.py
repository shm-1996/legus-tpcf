from header import *

from Galaxy import Galaxy
# from Galaxy import Galaxy
from Plot_Class import myPlot


def pipeline_galaxy(galaxy_name,method='masked',
    outdir=None,overwrite=False,plot=False):
    """
    Perform all required operations for one galaxy.
    Parameters
        ----------
        galaxy_name : string
            Name of galaxy
        method: string
            Method to prepare random catalog.
        outdir: string
            Output directory to store summary/plots
        overwrite: Boolean
            Flag to overwrite summary file if it already exists  
        plot : Boolean 
            Flag to plot TPCF individually + other diagnostic plots      
        Returns
        -------
        galaxy_class : Class Galaxy
            Instance of class galaxy is returned
    """
    import itertools
    print("Performing TPCF pipeline for galaxy {}.".format(galaxy_name))
    
    #Initialise galaxy class
    galaxy_class = Galaxy(galaxy_name,verbose=True)

    #Change output directory if passed by user
    if(outdir is not None):
        galaxy_class.outdir = outdir
        
    if(method not in ['masked','uniform','masked_radial','exponential']):
        raise myError("Method not recognised.")
    
    #Make sure path exists else create
    if(not os.path.exists(galaxy_class.outdir)):
        print("Path for output {} do not exist. Creating now"
            .format(galaxy_class.outdir))
        os.makedirs(galaxy_class.outdir)

     #First check if pickle file of class already exists here
    if(overwrite is False):
        if(os.path.isfile(galaxy_class.outdir+
            '{}_summary.pkl'.format(galaxy_class.name))):
            print("Summary file exists and overwrite=False. Reading existing file.")
            galaxy_class = loadObj(galaxy_class.outdir+
            '{}_summary'.format(galaxy_class.name))
            return galaxy_class
        else :
            print("TPCF not saved in output directory . Computing TPCF now..")    
    else :
        print("Overwrite flag is provided. Computing and overwriting TPCF even if it exists.")

    

    #Pipeline for each galaxy starts now - 

    #TODO: Implement this

    ###################################################################################
    # 1. Compute TPCF for -
    if(method == 'exponential'):
        #Extract old TPCF r_c fitted
        sampler = loadObj(masked_directory+'MCMC_sampler')
        samples = sampler.flatchain
        PF_median_fit = np.percentile(samples,50,axis=0)
        r_c = PF_median_fit[2]
        galaxy_class.rc = r_c

    # a. all clusters
    galaxy_class.Compute_TPCF(random_method=method,age=None)
    # b. young clusters (T<10 Myr)
    galaxy_class.Compute_TPCF(random_method=method,age='young',age_cut=1.e7)
    # c. old clusters (T>10 Myr)
    galaxy_class.Compute_TPCF(random_method=method,age='old',age_cut=1.e7)

    #Plot TPCF if required
    if(plot):
        pl = myPlot(galaxy_class)
        try:
            if(age == 'young'):
                pl.plot_TPCF(save=save,function=None,omega1=True,
                    filename=pl.galaxy.outdir+'{}_Young_TPCF.pdf'.format(galaxy_name))
            elif(age == 'old'):
                pl.plot_TPCF(save=save,function=None,omega1=True,
                    filename=pl.galaxy.outdir+'{}_Old_TPCF.pdf'.format(galaxy_name))
            else:
                pl.plot_TPCF(save=save,function=None,omega1=True,
                    filename=pl.galaxy.outdir+'{}_All_TPCF.pdf'.format(galaxy_name))
        except:
            print("TPCF could not be plotted.")
    saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'.format(galaxy_class.name))
    ###################################################################################
    # 3. Fit TPCF with MCMC - both young and old all 3 functional forms
    for func,age in itertools.product(['singlepl','piecewise','singletrunc'],
        ['young','old']):
        galaxy_class.fit_TPCF(method='mcmc',function=func,age=age,plot=plot)
    saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'.format(galaxy_class.name))

    ###################################################################################
    # 4. Compare AIC for 3 models
    galaxy_class.compare_AIC(age='young')
    galaxy_class.compare_AIC(age='old')

    ###################################################################################
    # 5. Compute Inferred Physical properties: lcorr, D2 and ExpDiskrc
    galaxy_class.getPhysicalProps()

    ###################################################################################
    # 6. Save galaxy class as pickle file
    print("Saving class object of {} as pickle file.".format(galaxy_class.name))
    saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'.format(galaxy_class.name))
    print("\n\n")
    

    return galaxy_class

def plots_galaxy(pl,method='masked_radial',function='piecewise',outdir=None,save=False,age=None):
    """
    Plot all the required output plots for each galaxy.
    Parameters
        ----------
        pl : class myPlot
            Instance of class myPlot.
        method: string
            Method to prepare random catalog.
        Returns
        -------
        None
    """

    if(outdir is not None):
        galaxy_class.outdir = outdir

    print("Plotting summary plots for Galaxy {}".format(pl.galaxy.name))
    #pl.galaxy_image(save=True)
    #pl.plot_clusters(save=save)
    pl.plot_random(save=save,random_method=method)
    #pl.class_distribution(save=save)
    # try:
    #     pl.mass_histogram(save=save)
    # except :
    #     print("Issue with mass histogram. Cannot create image.")
    # try:
    #     pl.age_histogram(save=save)
    # except :
    #     print("Issue with age histogram. Cannot create image.")
    # try:
    #     pl.bin_distribution(save=save)
    # except:
    #     print("Some bins were empty so did not analyse bin dist for this galaxy.")
    try:
        if(age == 'young'):
            pl.plot_TPCF(save=save,function=None,omega1=True,filename=pl.galaxy.outdir+'Young_TPCF.pdf')
        elif(age == 'old'):
            pl.plot_TPCF(save=save,function=None,omega1=True,filename=pl.galaxy.outdir+'Old_TPCF.pdf')
        else:
            pl.plot_TPCF(save=save,function=None,omega1=True)
    except:
        print("TPCF could not be plotted.")
    # try:
    #     pl.plot_TPCF_allclass(random_method=method,save=save,verbose=True)
    # except: 
    #     print("Some classes not available for this galaxy. Check class hist.")


def tpcf_allgalaxies(method,function,overwrite=False,save=False,age=None) :
    """
    Plot all the required output plots for all galaxies.
    Parameters
        ----------
        method : string
            Method to prepare random catalog.
        
        Returns
        -------
        None
    """
    print("Computing TPCF for all galaxies in LEGUS Survey.")
    for galaxy_name in list_of_galaxies:
        print("Calculating TPCF for galaxy {}".format(galaxy_name))
        galaxy_class = tpcf_galaxy(galaxy_name,method=method,function=function,overwrite=overwrite,
            age=age)
        plot_class = myPlot(galaxy_class)
        plots_galaxy(plot_class,method=method,function=function,save=save,age=age)
        #Save galaxy class info as pickle file
        path_exists = os.path.isfile(galaxy_class.outdir+
            '{}_summary.pkl'.format(galaxy_class.name))
        if(overwrite is True or path_exists is False):
            print("Saving class object of {} as pickle file.".format(galaxy_class.name))
            if(age is None):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'
                    .format(galaxy_class.name))            
            elif(age == 'young'):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_young_summary'
                    .format(galaxy_class.name))

            elif(age == 'old'):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_old_summary'
                    .format(galaxy_class.name))
        print("##############################################################")
        print("\n\n")


if __name__ == "__main__":

    #Parsing Arguments
    ############################################################################
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for tpcf-starclusters. All inputs optional. ')
    ap.add_argument('-method',action='store',type=str,default='masked_radial',
        help='Method to prepare the random catalog: "Uniform","Masked"' +
        '" Masked_radial (default)" ')
    ap.add_argument('-galaxy',action='store',type=str,default=None,
        help = 'Galaxy for which tpcf to be computed. By default done for all.')
    ap.add_argument('-outdir',action='store',type=str,default=None,
        help = 'Alternate output directory for plots and files.')
    ap.add_argument('-overwrite',action='store_true',
        help='overwrite computation of TPCF.Default: False.')
    ap.add_argument('-function',action='store',type=str,default='piecewise',
        help='Functional form to fit to TPCF: "piecewise" and "smooth" ')
    ap.add_argument('-show_only',action='store_true',
        help='Invoke this flag to only show the plots and not save them.')
    ap.add_argument('-age',action='store',type=str,default=None,
        help = 'Use for age cuts on clusters. Accepted values are young and old.')
    args = vars(ap.parse_args())

    method = args['method'].lower()
    if(method not in ['uniform','masked','masked_radial','exponential']):
        raise ValueError("This method does not exist. Allowed values are "+
            "'Uniform', 'Masked', and 'Masked_Radial'.")

    function = args['function'].lower()
    if(function not in ['piecewise','smooth']):
        raise ValueError("This functional form does not exist.")
    galaxy_name = args['galaxy'].upper()
    if(args['outdir'] is not None):
        output_directory = os.path.abspath(args['outdir']+'/')
        if(not os.path.exists(output_directory)):
            os.makedirs(output_directory)

    else:
        output_directory = None
    #Arguments parsed
    ############################################################################

    #Check if plots to be saved
    if(args['show_only'] is True):
        save=False
    else:
        save=True

    if(args['age'] is not None):
        
        if(args['age'].lower() == 'young'):
            age = 'young'
        elif(args['age'].lower() == 'old'):
            age = 'old'
        else:
            raise ValueError("Age flag should have value old or young. Please check.")
    else:
        age = None

    #Only for one galaxy
    if(galaxy_name is not None):
         #for all galaxies
        if(galaxy_name == 'ALL'):
            tpcf_allgalaxies(method,function=function,overwrite=args['overwrite'],save=save,
                age=age)
        elif(galaxy_name in list_of_galaxies):

            print("Running tpcf-starclusters for {}.".format(galaxy_name))
            galaxy_class = tpcf_galaxy(galaxy_name,method=method,outdir=output_directory,function=function,
                overwrite=args['overwrite'],age=age)
            plot_class = myPlot(galaxy_class)
            plots_galaxy(plot_class,method=method,outdir=output_directory,function=function,save=save,age=age)
            print("Saving class object of {} as pickle file.".format(galaxy_class.name))
            if(age is None):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'
                    .format(galaxy_class.name))            
            elif(age == 'young'):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_young_summary'
                    .format(galaxy_class.name))

            elif(age == 'old'):
                saveObj(galaxy_class,galaxy_class.outdir+'{}_old_summary'
                    .format(galaxy_class.name))
            #Save galaxy class info as pickle file



        else:
            raise myError("The provided galaxy {} is not in the list of galaxies".format(galaxy_name)+
                " for which cluster catalogs are available with LEGUS.")

   

    else :
        raise myError("Please provide galaxy name or ALL for all galaxies.")

        #TODO: Prepare galaxy info table
        #galaxy_info(output_directory)

        #TODO: Prepare fit summary table
        #galaxy_fit(output_directory)






    
        
        



