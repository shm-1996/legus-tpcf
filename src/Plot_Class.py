from header import *
from TPCF import linear_function,linear_truncation,\
smooth_function,onepowerlaw_function,piecewise_truncation
import regions

# from TPCF import bootstrap_two_point,bootstrap_two_point_angular,linear_function,linear_truncation,\
# smooth_function,onepowerlaw_function,piecewise_truncation
#from Galaxy import Galaxy
from sklearn.neighbors import KDTree


class myPlot():
    """
    Class that provides functionality for different plots for a galaxy
    Parameters
        galaxy: class Galaxy
            Instance of class galaxy that contains its properties
    """
    def __init__(self,galaxy):
        self.galaxy = galaxy
        
        
        
    def plot_TPCF(self,age='both',function=None,save=False,filename=None,omega1=True,
        axs = None,sec_axis=True,**kwargs):
        """
        Plot TPCF of a galaxy

        Parameters
        ----------
        age : string
            Age subset to plot, i.e. young or old. Can be None to plot combined, or both.
        function : string
            Fitted functional form to overplot, can be None, 'singlepl','piecewise',
            'singletrunc' or 'best', for plotting the best-fit among the three
        save : Boolean
            Flag to save the plot as a PDF. 
        filename: string
            Filename to save plot as. 
        omega1 : Boolean
            Plot 1+Omega, if false plots Omega
        axs : Matplotlib axes instance
            Plot in this axis. Creates a new axis if not provided. 
        sec_axis: Boolean
            Flag to plot the secondary axis for linear distance
        **kwargs: 
            Other optional parameters for plotting that go directly into formatting.
        
        Returns
        -------
        None
        
        """
        if(axs is None):
            fig,axs = plt.subplots(ncols=1)
        if(sec_axis is True):
            ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        separation_bins = self.galaxy.bin_centres*(1./arcsec_to_degree)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        if(function not in ['piecewise','best','singlepl','singletrunc',None]):
            raise ValueError("This funtional form does not exist.")
        if(age not in ['young','old','both',None]):
            raise ValueError("This age group cannot be plotted.")

        lw = kwargs.pop('lw',0.5)
        if(age in ['young','both']):
            lc = kwargs.pop('c','#F56B5C')
        elif(age == 'old'):
            lc = kwargs.pop('c','#4591F5')
        else:
            lc = kwargs.pop('c','k')
        fmt = kwargs.pop('fmt','.')
        
            
        if(age is None):
            separation_bins,corr_fit,dcorr_fit = filter_bins(self.galaxy.bin_centres,
                self.galaxy.corr,self.galaxy.dcorr)
            axs.errorbar(separation_bins,1+corr_fit,yerr=dcorr_fit,
            fmt=fmt,lw=lw,c=lc,**kwargs)
        if(age in ['both','young']):
            separation_bins,corr_fit,dcorr_fit = filter_bins(self.galaxy.bin_centres,
                self.galaxy.ycorr,self.galaxy.ydcorr)
            axs.errorbar(separation_bins,1+corr_fit,yerr=dcorr_fit,
            fmt=fmt,lw=lw,c=lc,label=r'$T <= 10 \, \mathrm{Myr}$')
        if(age in ['both','old']):
            if(age == 'both'):
                lc = kwargs.pop('c','#4591F5')
            separation_bins,corr_fit,dcorr_fit = filter_bins(self.galaxy.bin_centres,
                self.galaxy.ocorr,self.galaxy.odcorr)
            axs.errorbar(separation_bins,1+corr_fit,yerr=dcorr_fit,
            fmt=fmt,lw=lw,c=lc,label=r'$T > 10 \, \mathrm{Myr}$')

        if(age in ['both','young']):

            ls = '--'
            lw = 2.0
            lc = '#F56B5C'

            if(function == 'best'):
                best_model = self.galaxy.best_model_y
                if(best_model == 'piecewise'):
                    axs.plot(plot_points,np.exp(linear_function(plot_points,self.galaxy.best_fit_y[0],
                        self.galaxy.best_fit_y[1],self.galaxy.best_fit_y[2],self.galaxy.best_fit_y[3])),
                        ls = ls,lc=lc,lw=lw)
                elif(best_model == 'singlepl'):
                    axs.plot(plot_points,np.exp(onepowerlaw_function(plot_points,self.galaxy.best_fit_y[0],
                        self.galaxy.best_fit_y[1]),ls=ls,color=lc,lw=lw))
                elif(best_model == 'singletrunc'):
                    axs.plot(plot_points,np.exp(linear_truncation(plot_points,self.galaxy.best_fit_y[0],
                        self.galaxy.best_fit_y[1],self.galaxy.best_fit_y[2]),ls=ls,color=lc,lw=lw))
            
            elif(function == 'piecewise'):
                axs.plot(plot_points,np.exp(linear_function(plot_points,
                    self.galaxy.fit_values_piecewise_y[0],self.galaxy.fit_values_piecewise_y[1],
                    self.galaxy.fit_values_piecewise_y[2],self.galaxy.fit_values_piecewise_y[3])),
                    ls=ls,lc=lc,lw=lw)
                break_theta = np.exp(self.galaxy.fit_values_piecewise_y[3])
                break_theta_error = np.exp(self.galaxy.fit_values_piecewise_y[3])
                
            elif(function == 'singlepl'):
                axs.plot(plot_points,np.exp(onepowerlaw_function(plot_points,
                    self.galaxy.fit_values_singlepl_y[0],self.galaxy.fit_values_singlepl_y[1])),
                    ls=ls,lc=lc,lw=lw)
            elif(function == 'singletrunc'):
                axs.plot(plot_points,np.exp(linear_truncation(plot_points,
                    self.fit_values_singletrunc_y[0],self.fit_values_singletrunc_y[1],
                    self.fit_values_singletrunc_y[2])),
                    ls=ls,lc=lc,lw=lw)
                theta_c = self.galaxy.fit_values[2]
                theta_c_error = self.galaxy.fit_errors[2]


        elif(age in ['both','old']):
            ls = '--'
            lw = 2.0
            lc = '#4591F5'

            if(function == 'best'):
                best_model = self.galaxy.best_model_o
                if(best_model == 'piecewise'):
                    axs.plot(plot_points,np.exp(linear_function(plot_points,self.galaxy.best_fit_o[0],
                        self.galaxy.best_fit_o[1],self.galaxy.best_fit_o[2],self.galaxy.best_fit_o[3])),
                        ls = ls,lc=lc,lw=lw)
                elif(best_model == 'singlepl'):
                    axs.plot(plot_points,np.exp(onepowerlaw_function(plot_points,self.galaxy.best_fit_o[0],
                        self.galaxy.best_fit_o[1]),ls=ls,color=lc,lw=lw))
                elif(best_model == 'singletrunc'):
                    axs.plot(plot_points,np.exp(linear_truncation(plot_points,self.galaxy.best_fit_o[0],
                        self.galaxy.best_fit_o[1],self.galaxy.best_fit_o[2]),ls=ls,color=lc,lw=lw))
            
            elif(function == 'piecewise'):
                axs.plot(plot_points,np.exp(linear_function(plot_points,
                    self.galaxy.fit_values_piecewise_o[0],self.galaxy.fit_values_piecewise_o[1],
                    self.galaxy.fit_values_piecewise_o[2],self.galaxy.fit_values_piecewise_o[3])),
                    ls=ls,lc=lc,lw=lw)
                
            elif(function == 'singlepl'):
                axs.plot(plot_points,np.exp(onepowerlaw_function(plot_points,
                    self.galaxy.fit_values_singlepl_o[0],self.galaxy.fit_values_singlepl_o[1])),
                    ls=ls,lc=lc,lw=lw)
            elif(function == 'singletrunc'):
                axs.plot(plot_points,np.exp(linear_truncation(plot_points,
                    self.fit_values_singletrunc_o[0],self.fit_values_singletrunc_o[1],
                    self.fit_values_singletrunc_o[2])),
                    ls=ls,lc=lc,lw=lw)
                

        
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        if(omega1 is True):
            axs.set_ylabel(r"$1+\omega\left(\theta \right)$")
        else:
            axs.set_ylabel(r"$\omega\left(\theta \right)$")
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.legend()
        if(sec_axis is True):
            #axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
            ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_TPCF.pdf'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            return


    def plot_TPCF_allclass(self,random_method = 'masked_radial',
        save=False,filename=None,verbose=False):
        """
        Plot TPCF of a galaxy for all classes comparing between them

       Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        random_method: string
            random method to use
        filename : string
            File to save to, else default filename 
        verbose : string
            Whether to print verbose
        
        Returns
        -------
        None
        
        """
        if(verbose):
            print("Plotting comparison of TPCF for different classes.")
        #Initialise figure
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        

        #Compute TPCF for each class
        for i in range(1,4):
            if(verbose):
                print("Computing TPCF for class {} clusters".format(i))
            self.galaxy.Compute_TPCF(cluster_class=i,random_method=random_method,
                verbose=verbose)
            separation_bins = self.galaxy.bin_centres*(1./arcsec_to_degree)        
            indices = np.where(self.galaxy.corr>0.0)
            corr_fit = self.galaxy.corr[indices].astype(np.float)
            dcorr_fit = self.galaxy.dcorr[indices].astype(np.float)
            separation_bins = separation_bins[indices].astype(np.float)
            axs.errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
                fmt='.-',label='Class {}'.format(i))
        # Combined
        if(verbose):
            print("Computing TPCF for combined class clusters")

        #TODO: Figure out how to do this. Currently the properties for 
        self.galaxy.Compute_TPCF(cluster_class=-1,random_method=random_method,
            verbose=verbose)
        separation_bins = self.galaxy.bin_centres*(1./arcsec_to_degree)
        indices = np.where(self.galaxy.corr>0.0)
        corr_fit = self.galaxy.corr[indices].astype(np.float)
        dcorr_fit = self.galaxy.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        axs.errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',label='Class 1+2+3')


        #Rest of plot
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs.set_xscale('log')
        axs.set_yscale('log')

        #Secondary axis
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        
        axs.legend()
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_Classes_TPCF.pdf'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def plot_clusters(self,save=False,filename=None):
        """
        Plot spatial distribution of clusters in the galaxy

       Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 
        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)

        ages = self.galaxy.get_cluster_ages()
        cmap = cmr.waterlily

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        x = self.galaxy.ra*3600.0
        y = self.galaxy.dec*3600.0

        #Get central pixel values
        xcen = (np.min(x)+np.max(x))/2.
        ycen = (np.min(y)+np.max(y))/2.

        #Scale offset around bounding box to ~ 5% of axes width
        offset_x = (np.max(x)-np.min(x))*0.05
        offset_y = (np.max(y)-np.min(y))*0.05

        xmin,xmax = np.min(x)-offset_x-xcen, np.max(x)+offset_x-xcen
        ymin,ymax = np.min(y)-offset_y-ycen, np.max(y)+offset_y-ycen
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)

        im1 = ax1.scatter(x-xcen,y-ycen,c=np.log10(ages),alpha=0.5,cmap=cmap)
        cbar = fig.colorbar(im1,ax = ax1,use_gridspec=False,
                            orientation='vertical',pad=0.00,aspect=30)

        # #Draw 100 arcsec scale bar
        import matplotlib.lines as lines
        #No of pixels in axes
        total_pixels = np.int(np.floor(ax1.transData.transform((xmax,ymax))[0]) - \
        np.floor(ax1.transData.transform((xmin,ymin))[0]))

        length_per_pixel = (xmax-xmin)/(total_pixels)
        #Convert to parsec 
        length_per_pixel = self.sep_to_pc(length_per_pixel)
        #Scale bar of 50 arcsec
        length = self.sep_to_pc(50)
        no_pixels = np.abs(length/length_per_pixel)
        no_pixels = no_pixels/total_pixels

        scale = lines.Line2D([0.8,0.8+no_pixels],[0.1],
                                         lw=1,color='black',
                                        transform=ax1.transAxes)

        ax1.add_line(scale)
        ax1.annotate(r'$50^{\prime \prime} = %d \, \mathrm{pc}$'%length,(0.65,0.15),xycoords='axes fraction',
                            fontsize=12)

        ax1.set_xlabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)
        ax1.set_ylabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_clusters'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()




    def plot_random(self,save=False,random_method='masked_radial',filename=None):
        """
        Plot spatial distribution of random catalog
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        random_method : string
            Random method in use
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 
        from TPCF import uniform_sphere,masked_random_sample,\
        masked_random_sample,exponential_sample
        #Obtain one instance of the random sample
        if(random_method == 'uniform') :
            #Draw uniform random sample with N points
            ra_R, dec_R = uniform_sphere(self.galaxy,2 * len(self.galaxy.ra))
                                      
        elif(random_method == 'masked') :
            #Draw a masked random sample with N points
            ra_R, dec_R = masked_random_sample(self.galaxy,10*len(self.galaxy.ra))

        elif(random_method == 'masked_radial') :
            #Draw a random sample 
            ra_R, dec_R = masked_radial_random_sample(self.galaxy,1000*len(self.galaxy.ra))
        elif(random_method == 'exponential'):
            r_c = self.galaxy.rc
            ra_R, dec_R = exponential_sample(self.galaxy,r_c=r_c,len_random=len(self.galaxy.ra))

        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)

        fig = plt.figure(constrained_layout=True)
        ax1 = fig.add_subplot(111)

        x = ra_R*3600.0
        y = dec_R*3600.0

        #Get central pixel values
        xcen = (np.min(x)+np.max(x))/2.
        ycen = (np.min(y)+np.max(y))/2.

        #Scale offset around bounding box to ~ 5% of axes width
        offset_x = (np.max(x)-np.min(x))*0.05
        offset_y = (np.max(y)-np.min(y))*0.05

        xmin,xmax = np.min(x)-offset_x-xcen, np.max(x)+offset_x-xcen
        ymin,ymax = np.min(y)-offset_y-ycen, np.max(y)+offset_y-ycen
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)

        im1 = ax1.scatter(x-xcen,y-ycen)

        ax1.set_xlabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)
        ax1.set_ylabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_random'.format(self.galaxy.name) + \
                '_{}'.format(random_method)
            plt.savefig(filename,bbox_inches='tight')
        else :
            plt.show()

    def class_distribution(self,save=False,filename=None):
        """
        Plot the distribution of clusters lying in different cluster classes. 
        Parameters:
        ----------
            save: Boolean
                Flag to save the plot
            filename: String
                File name of the plot
        Returns
        -------
        None
        
        """
        #Read file for distribution of classes
        file = np.loadtxt(self.galaxy.catalog_file)
        N0 = np.size(np.where(file[:,33]==0))
        N1 = np.size(np.where(file[:,33]==1))
        N2 = np.size(np.where(file[:,33]==2))
        N3 = np.size(np.where(file[:,33]==3))
        N4 = np.size(np.where(file[:,33]==4))

        #Plot now
        fig,axs = plt.subplots(ncols=1)
        label = ['Class 0', 'Class 1', 'Class 2', 'Class 3','Class 4']
        axs.bar([0,1,2,3,4],[N0,N1,N2,N3,N4],color='#F59005',tick_label=label)
        axs.set_xlabel('Cluster Class')
        axs.set_ylabel('Number')
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_ClassDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def mass_histogram(self,save=False,filename=None):
        """
        Plot distribution of masses in each class of clusters.
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 

        #Read file for distribution of masses in each class
        file = np.loadtxt(self.galaxy.catalog_file)
        M1 = file[np.where(file[:,33]==1)][:,19]
        M2 = file[np.where(file[:,33]==2)][:,19]
        M3 = file[np.where(file[:,33]==3)][:,19]
        
        #Plot now
        fig,axs = plt.subplots(ncols=1)
        axs.hist(np.log10(M1),bins='auto',histtype='step',log=True,
            label='Class 1',color='#F51557')
        axs.hist(np.log10(M2),bins='auto',histtype='step',log=True,
            label='Class 2',color='#33A7F4')
        axs.hist(np.log10(M3),bins='auto',histtype='step',log=True,
            label='Class 3',color='#28F56E')

        axs.set_xlabel(r'$\log_{10} \, \mathrm{Mass} \, (M_{\odot})$')
        axs.set_ylabel(r'$\mathrm{Number}$')
        axs.legend()
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_MassDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def age_histogram(self,save=False,filename=None):
        """
        Plot distribution of ages in each class of clusters.
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 

        #Read file for distribution of ages in each class
        file = np.loadtxt(self.galaxy.catalog_file)
        A1 = file[np.where(file[:,33]==1)][:,16]
        A2 = file[np.where(file[:,33]==2)][:,16]
        A3 = file[np.where(file[:,33]==3)][:,16]

        # # Ages are in yr. Convert to Myr.
        # A1 /= 1.e6
        # A2 /= 1.e6
        # A3 /= 1.e6

        
        #Plot now
        fig,axs = plt.subplots(ncols=1)
        axs.hist(np.log10(A1),bins='auto',histtype='step',log=True,
            label='Class 1',color='#F51557')
        axs.hist(np.log10(A2),bins='auto',histtype='step',log=True,
            label='Class 2',color='#33A7F4')
        axs.hist(np.log10(A3),bins='auto',histtype='step',log=True,
            label='Class 3',color='#28F56E')

        axs.set_xlabel(r'$\log_{10} \, \mathrm{Age} \, (\mathrm{Myr})$')
        axs.set_ylabel(r'$\mathrm{Number}$')
        axs.legend()
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_AgeDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()

    def bin_distribution(self,save=False,filename=None):
        """
        Plot distribution of pairs in each TPCF bin. 
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        Returns
        -------
        None
        
        """         

        #Get no of pairs for each class using KDTree two point correlation
        counts_DD = np.zeros((3,self.galaxy.no_bins))
        bins_transform = angular_dist_to_euclidean_dist(self.galaxy.bins)
        for i in range(1,4):
            self.galaxy.get_ra_dec(cluster_class=i)
            xyz_clusters = np.asarray(ra_dec_to_xyz(self.galaxy.ra, 
                                            self.galaxy.dec), order='F').T
            KDT_D = KDTree(xyz_clusters)
            counts = KDT_D.two_point_correlation(xyz_clusters, bins_transform)
            counts_DD[i-1] = np.diff(counts)
        #Reset RA/DEC
        self.galaxy.get_ra_dec(cluster_class=-1)    

        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[0],
            align='edge',width=np.diff(self.galaxy.bins_arcsec),
            color='#F51557',label='Class 1')
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[1],
                    align='edge',width=np.diff(self.galaxy.bins_arcsec),
                    bottom=counts_DD[0],color='#4983FC',label='Class 2')
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[2],
                    align='edge',width=np.diff(self.galaxy.bins_arcsec),
                    bottom=counts_DD[1],color='#FAC90E',label='Class 3')
            
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r'$\mathrm{N_{\mathrm{pairs}}}$')
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.legend(loc='upper left')
        
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_BinDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def hst_image(self,save=False,filename=None):

        """
        Plot optical HST image of galaxy. 
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """         

        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)
        fig = plt.figure(constrained_layout=True)
        ax1 = fig.add_subplot(111,projection=wcs)

        with np.errstate(divide='ignore', invalid='ignore'):
            im = ax1.imshow(np.log10(hdu.data),vmin=-2.0)
        cbar = fig.colorbar(im,ax = ax1,orientation='vertical')
        cbar.ax.set_ylabel(r"$\log_{10} \, I$",rotation=90,labelpad=5,fontsize=20)
        ax1.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        ax1.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)

        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_HSTImage'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()

    def massage_image(self,save=False,filename=None):
        """
        Overplot clusters on optical HST image of galaxy, colored by age and mass. 
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """  

        hdu = fits.open(self.galaxy.fits_file)[0]
        hdu.data[np.where(hdu.data <= 1.e-1)] = np.nan

        #Read ages and masses
        ages = self.galaxy.get_cluster_ages()
        masses = self.galaxy.get_cluster_masses()

        #Colormap for scatter
        cmap = cmr.waterlily

        #Convert ra/dec to pixel coordinates
        wcs_galaxy = WCS(hdu.header)
        xpix,ypix = wcs_galaxy.all_world2pix(self.galaxy.ra_raw,self.galaxy.dec_raw,0)

       
        fig = plt.figure(tight_layout=False,figsize=(8,4))
        
        ax1 = fig.add_subplot(121,projection=wcs_galaxy)

        #Plot HST image
        with np.errstate(divide='ignore', invalid='ignore'):
            im = ax1.imshow(np.log10(hdu.data),vmin=-2.0,alpha=0.4)
        

        #Plot scatter
        
        im1 = ax1.scatter(xpix,ypix,c=np.log10(ages),alpha=0.6,cmap=cmap,lw=0.3)
        cbar = fig.colorbar(im1,ax = ax1,use_gridspec=False,
                            orientation='vertical',pad=0.00,aspect=30)
        cbar.ax.set_ylabel(r"$\log_{10} \, \mathrm{Age} \, (\mathrm{yr})$",rotation=90,labelpad=5,fontsize=20)
        ax1.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=16)
        ax1.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-1.2,
                       fontsize=16)
        #Find extents of clusters
        ra_min,ra_max = np.min(self.galaxy.ra_raw),np.max(self.galaxy.ra_raw)
        dec_min,dec_max = np.min(self.galaxy.dec_raw),np.max(self.galaxy.dec_raw)
        imax,jmin = np.floor(wcs_galaxy.all_world2pix(ra_min,dec_min,0)).astype(int)
        imin,jmax = np.floor(wcs_galaxy.all_world2pix(ra_max,dec_max,0)).astype(int)


        #Set axis limits to these limits
        ax1.set_xlim(imin-500,imax+500)
        ax1.set_ylim(jmin-500,jmax+500)
        scale = self.scale_to_plot(0.1,0.1,ax1)
        ax1.add_line(scale)
        ax1.annotate(r'$1 \, \mathrm{kpc}$',(0.1,0.05),xycoords='axes fraction',
                    fontsize=12)   


        #####################################################################
        ax2 = fig.add_subplot(122,projection=wcs_galaxy)

        with np.errstate(divide='ignore', invalid='ignore'):
            im = ax2.imshow(np.log10(hdu.data),vmin=-2.0,alpha=0.4)
        im1 = ax2.scatter(xpix,ypix,c=np.log10(masses),alpha=0.6,cmap=cmap,lw=0.3)
        cbar = fig.colorbar(im1,ax = ax2,use_gridspec=False,
                            orientation='vertical',pad=0.00,aspect=30)
        cbar.ax.set_ylabel(r"$\log_{10} \, \mathrm{Mass} \, (M_{\odot})$",rotation=90,labelpad=5,fontsize=20)
        dec = ax2.coords['DEC']
        dec.set_ticklabel_visible(False)
        ax2.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=16)

        #Set axis limits to these limits
        ax2.set_xlim(imin-500,imax+500)
        ax2.set_ylim(jmin-500,jmax+500)
        scale = self.scale_to_plot(0.1,0.1,ax2)
        ax2.add_line(scale)
        ax2.annotate(r'$1 \, \mathrm{kpc}$',(0.1,0.05),xycoords='axes fraction',
                    fontsize=12)    


        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_galaxyImage'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def scale_to_plot(self,xcen,ycen,axs,length=1.e3) :
        """
        Return scale bar to add on plot.    

        Returns
        -------
        scale_bar: lines.line2D object
        """
        import matplotlib.lines as lines    
        
        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs_galaxy = WCS(hdu.header)

        #Find length per pixel
        xmin,xmax = axs.get_xlim()
        ymin,ymax = axs.get_ylim()
        total_pixels = xmax-xmin

        ra_min = wcs_galaxy.all_pix2world(xmin,ymin,0)[0]
        ra_max = wcs_galaxy.all_pix2world(xmax,ymin,0)[0]
        length_per_pixel = np.abs(ra_max-ra_min)*3600./(total_pixels)
        length_per_pixel = self.sep_to_pc(length_per_pixel)
        
        no_of_pixels = np.abs(length/length_per_pixel)
        no_of_pixels = no_of_pixels/total_pixels
        scale_bar = lines.Line2D([xcen,xcen+no_of_pixels],[ycen],
                                 lw=1,color='black',
                                transform=axs.transAxes) 
        return scale_bar


    def sep_to_pc(self,sep):
        """
        Angular separation in arcsec to parsec

        Parameters
        ----------
        sep : float
            angular separation in arcsec
        Returns
        -------
        pc : float
             The separation in parsec for given angular separation
        """

        distance = self.galaxy.distance*const.Parsec*1.e6
        arcsec_to_pc = u.arcsec.to(u.radian)*distance/(const.Parsec)
        return sep*arcsec_to_pc
    def pc_to_sep(self,pc) :
        """
        Linear separation in parsec to arcsec

        Parameters
        ----------
        pc : float
             The separation in parsec 
        Returns
        -------
        sep : float
            angular separation in arcsec
        """
        distance = self.galaxy.distance*const.Parsec*1.e6
        pc_to_radian = pc*const.Parsec/distance
        return pc_to_radian*u.radian.to(u.arcsec)
    def axs_to_parsec(axs) :
        """
        Draw twin axes corresponding to angular limits
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

#Some useful utility functions for plotting

def filter_bins(bins,corr,dcorr):
    """
    Filter non-physical bins where the 1+omega is negative or the
    error is higher than the correlation. The returned quantities
    can be directly plotted. 
    Parameters
    ----------
    bins : ndarray
        Array of bins where TPCF is calculated
    corr: ndarray
        Array of TPCF values at these bins
    dcorr: ndarray
        Array of error on the above TPCF values
    Returns
    -------
    separation_bins : ndarray
        Filtered bins
    corr : ndarray
        corr_fit : ndarray
        Filtered corr function values
    dcorr : ndarray
        Filtered error values
    """

    #Isolate non-zero correlation points
    indices_err = np.abs(corr)>np.abs(dcorr)
    indices_neg = corr>-1
    indices = np.where(np.logical_and(indices_err,indices_neg))
    corr_fit = corr[indices].astype(np.float)
    dcorr_fit = dcorr[indices].astype(np.float)
    separation_bins = bins*(1./arcsec_to_degree)
    separation_bins = separation_bins.astype(np.float)
    separation_bins = separation_bins[indices].astype(np.float)
    separation_bins = separation_bins.astype(np.float)
    return separation_bins, corr_fit, dcorr_fit

def bbox(img):
    """
    Filter out pixels in an image outside the bounding box of non-zero pixel values. 

    Parameters
    ----------
    image : 2D ndarray
        Array containing the pixel values of the image
    Returns
    -------
    rmin : integer
        Minimum row index
    rmax : integer
        Maximum row index
    cmin : integer
        Minimum column index
    cmin : integer
        Maximum column index
    """
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return rmin, rmax, cmin, cmax

def deproject_region_centre(region,i,xpix_c,ypix_c,galaxy_class):
    """
    Deproject the region file of the galaxy about the galaxy centre. 
    Parameters:
    -------
        region: ds9 region Instance array
            Array of Region instance of the galaxy to deproject
        i : integer
            The index containing the specific region to deproject
        xpix_c : float
            Centre of the galaxy in pixel coordinates
        ypix_c : float
            Centre of the galaxy in pixel coordinates
        galaxy_class: Instance of Class Galaxy
            The instance of the galaxy class
    Returns:
    -------
    None
    """
    #rotate clockwise by angle PA
    region_rotated = region[i].rotate(regions.PixCoord(xpix_c,ypix_c),-galaxy_class.pa*u.deg)
    try:
        size = np.size(region_rotated.vertices)
    except:
        return None
    x = np.zeros(size)
    y = np.zeros(size)
    for i in range(0,size):
        x[i] = region_rotated.vertices[i].x/np.cos(np.deg2rad(galaxy_class.inclination)) -\
        xpix_c/np.cos(np.deg2rad(galaxy_class.inclination)) + xpix_c
        y[i] = region_rotated.vertices[i].y
    regions_dep = regions.PolygonPixelRegion(vertices=regions.PixCoord(x=x,y=y))    
    return regions_dep

def get_separations(sides,pl,parsec=True):
    """
    Get the length of the sides of the polygo FoV observed in the galaxy. 
    Parameters:
    -------
        sides: Tuple
            Tuple of the vertices of the FoV polygon 
        pl : Class Instance of Plot_Class
            The plot_class instance for this galaxy
        parsec : Boolean
            Flag to return sides in parsec units, otherwise in angular separation
    Returns:
    -------
        sizes : List
            List of the length of sides of the FoV polygon
    """
    i = 0
    sizes = []
    while i<np.size(sides):
        if(i == np.size(sides)-1):
            s = sides[0].separation(sides[i]).arcsec
        else:
            s = sides[i+1].separation(sides[i]).arcsec
        if(parsec):
            s = pl.sep_to_pc(s)
        sizes.append(s)
        i = i+1
    return sizes


