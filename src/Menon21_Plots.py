from header import *
from Plot_Class import *
import matplotlib.lines as lines
from regions import read_ds9
from astropy.coordinates import SkyCoord
import regions
from TPCF import linear_function,linear_truncation,\
        smooth_function,onepowerlaw_function,piecewise_truncation,\
        bootstrap_two_point,bootstrap_two_point_angular

#Axes limits in parsec
global_axes_limits = [8,1.e4]
#from PaperPlots import get_separations, deproject_region_centre
from matplotlib.patches import Patch



def Figure1(save=True,outdir='../results/',indir='../results/'):
    """
    Creates Figure 1 of Menon et al 2021b. 
    Parameters:
        save : Boolean
            Flag to save the figure as a pdf
        outdir: string
            Directory to save figure in
        indir: string
            Directory where galaxy summary pickle files lie
    """

    #Create figure and axs instance
    #fig,axs = plt.subplots(nrows=4,ncols=3,figsize=(12,16))
    print("Creating Figure 1")
    fig = plt.figure(figsize=(14,16),constrained_layout=True)

    if(indir == None):
        indir = os.path.abspath('../results/')+'/'
        
    #Directories
    indir = os.path.abspath(indir)+'/'
    outir = os.path.abspath(outdir)+'/'

    i,j = 0,0
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        
        galaxy_dir = indir+'/'            
        
        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        pl = myPlot(galaxy_class)

        image = fits.getdata(pl.galaxy.fits_file)
        #Find extents of clusters
        

        if(galaxy_name != 'NGC_5457'):
            image[np.where(image <= 1.e-2)] = np.nan

        wcs_galaxy = WCS(fits.getheader(pl.galaxy.fits_file))
        xpix,ypix = wcs_galaxy.all_world2pix(pl.galaxy.ra_raw,pl.galaxy.dec_raw,0)
        ages = galaxy_class.get_cluster_ages()

        axs = plt.subplot2grid((4,3),(i,j),projection=wcs_galaxy,fig=fig)       

        cmap = cmr.fusion_r

        with np.errstate(divide='ignore', invalid='ignore'):
            im = axs.imshow(np.log10(image),vmin=-2.0,alpha=0.4)
        im1 = axs.scatter(xpix,ypix,c=np.log10(ages),
            alpha=0.6,cmap=cmap,lw=0.3,s=8.0,edgecolors='none')
        cbar = fig.colorbar(im1,ax = axs,use_gridspec=False,
                    orientation='vertical',pad=0.00,aspect=30)

        imin,imax,jmin,jmax = bbox(~np.isnan(image))
        xextent,yextent = (imax-imin),(jmax-jmin)
        imin,imax = imin - xextent*0.1,imax + xextent*0.1
        jmin,jmax = jmin - yextent*0.1,jmax + yextent*0.1
        axs.set_ylim(imin,imax)
        axs.set_xlim(jmin,jmax)

        ra = axs.coords['RA']
        dec = axs.coords['DEC']
        ra.set_ticklabel(exclude_overlapping=True)
        dec.set_ticklabel(exclude_overlapping=True)
        ra.set_auto_axislabel(False)
        dec.set_auto_axislabel(False)
        

        # #Draw 50 arcsec scale bar
        #No of pixels in axes
        total_pixels = imax-imin
        xmin,ymin = wcs_galaxy.all_pix2world(imin,jmin,0)
        xmax,ymax = wcs_galaxy.all_pix2world(imax,jmax,0)

        #This is in degrees
        length_per_pixel = (xmax-xmin)/(total_pixels)
        #Convert to arcsec
        length_per_pixel = length_per_pixel*3600.
        #Convert to parsec 
        length_per_pixel = pl.sep_to_pc(length_per_pixel)
        length = pl.sep_to_pc(50)
        no_pixels = np.abs(length/length_per_pixel)
        no_pixels = no_pixels/total_pixels

        if(galaxy_name in ['NGC_5194','NGC_7793','NGC_5253']):
            scale = lines.Line2D([0.17,0.17+no_pixels],[0.07],
                                         lw=1,color='black',
                                        transform=axs.transAxes)
            axs.add_line(scale)
            axs.annotate(r'$50^{\prime \prime} = %d \, \mathrm{pc}$'%length,(0.07,0.15),
                xycoords='axes fraction',
                                fontsize=12)
        else:
            scale = lines.Line2D([0.7,0.7+no_pixels],[0.1],
                                         lw=1,color='black',
                                        transform=axs.transAxes)
            axs.add_line(scale)
            axs.annotate(r'$50^{\prime \prime} = %d \, \mathrm{pc}$'%length,(0.6,0.15),
                xycoords='axes fraction',
                                fontsize=12)

        if(j==2):
            cbar.ax.set_ylabel(r"$\log_{10} \, \mathrm{Age} \, (\mathrm{yr})$",
                rotation=90,labelpad=5,fontsize=20)
        if(j==0):
            axs.set_ylabel(r"$\mathrm{Declination \; (J2000)}$")
        else:
            dec.set_axislabel('')

        if(i == 3):
            axs.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$")
        else:
            ra.set_axislabel('')            

        axs.text(0.07,0.9,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs.transAxes)

        #Get position of subplot
        j +=1
        if(j==3):
            j = 0
            i +=1


    if(save):    
        filename = outdir+'Fig1.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


    return



def Figure2(save=False,outdir='../results/',indir='../results/',omega_neg=True):
    """
    Plot Figure 2 of Menon et al 2021b

    Parameters:
        save: Boolean
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
        omega_neg: Boolean
            Flag to allow TPCF values with omega<0 to be plotted
    Returns:
        None


    """
    from matplotlib.lines import Line2D
    print("Creating Figure 2")

    #Create figure and axs instance
    fig,axs = plt.subplots(nrows=4,ncols=3,figsize=(12,14))

    if(indir == None):
        indir = os.path.abspath('../results/')+'/'        

    #Directories
    indir = os.path.abspath(indir)+'/'
    outir = os.path.abspath(outdir)+'/'

    fit_labels = [r'$\mathrm{Model \; S}$',r'$\mathrm{Model \; PW}$',r'$\mathrm{Model \; PF}$']

    i,j = 0,0   
    
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        galaxy_dir = indir

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        plot_class = myPlot(galaxy_class)

        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        #Set axes
        axs[i,j].set_yscale('log')            
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        

        ########################################
        corr = galaxy_class.ycorr
        dcorr = galaxy_class.ydcorr    

        #Isolate non-zero correlation points
        indices_err = np.abs(corr)>np.abs(dcorr)
        if(omega_neg):
            indices_neg = corr>-1
        else:
            indices_neg = corr>0
        indices = np.where(np.logical_and(indices_err,indices_neg))
        corr_fit = corr[indices].astype(np.float)
        dcorr_fit = dcorr[indices].astype(np.float)
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        separation_bins = separation_bins.astype(np.float)

        galaxy_function = galaxy_class.best_model_y
        galaxy_AIC = [galaxy_class.aic_singlepl_y,galaxy_class.aic_piecewise_y,
            galaxy_class.aic_singletrunc_y]
        fit_label_young = fit_labels[np.argmin(galaxy_AIC)]
        
        
        axs[i,j].errorbar(separation_bins,1+corr_fit,yerr=dcorr_fit,
            fmt='.',lw=0.5,color='#4591F5',label=r'$T < 10 \, \mathrm{Myr}$')

        #Plot fits

        ls = '--'
        lw = 2.0
        lc = '#4591F5'
        if(galaxy_function == 'piecewise'):
            ls = '--'
            axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.best_fit_y[0],
                galaxy_class.best_fit_y[1],galaxy_class.best_fit_y[2],galaxy_class.best_fit_y[3])),
                ls=ls,color=lc,lw=lw)        
        elif(galaxy_function == 'singlepl'):
            ls = 'dotted'
            axs[i,j].plot(plot_points,np.exp(onepowerlaw_function(plot_points,galaxy_class.best_fit_y[0],
            galaxy_class.best_fit_y[1])),
            ls=ls,lw=lw,color=lc)        
        elif(galaxy_function == 'singletrunc'):
            ls = '-.'
            axs[i,j].plot(plot_points,np.exp(linear_truncation(plot_points,galaxy_class.best_fit_y[0],
                galaxy_class.best_fit_y[1],galaxy_class.best_fit_y[2])),
                ls=ls,lw=lw,color=lc)            
            #axs[i,j].axvline(theta_c,ls=':',label=r'$\theta_c = {:2.1f} \pm {:2.1f}$'.format(theta_c,
            #    theta_c_error))

        ls_y = ls

        ########################################
        if(galaxy_name not in ["NGC_3344","NGC_5253","NGC_7793"]):
            
            corr = galaxy_class.ocorr
            dcorr = galaxy_class.odcorr    

            #Isolate non-zero correlation points
            indices_err = np.abs(corr)>np.abs(dcorr)
            if(omega_neg):
                indices_neg = corr>-1
            else:
                indices_neg = corr>0
            indices = np.where(np.logical_and(indices_err,indices_neg))
            corr_fit = corr[indices].astype(np.float)
            dcorr_fit = dcorr[indices].astype(np.float)
            separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
            separation_bins = separation_bins.astype(np.float)
            separation_bins = separation_bins[indices].astype(np.float)
            separation_bins = separation_bins.astype(np.float)

            galaxy_function = galaxy_class.best_model_o
            galaxy_AIC = [galaxy_class.aic_singlepl_o,galaxy_class.aic_piecewise_o,
                galaxy_class.aic_singletrunc_o]
            fit_label_old = fit_labels[np.argmin(galaxy_AIC)]

            axs[i,j].errorbar(separation_bins,1+corr_fit,yerr=dcorr_fit,
                fmt='.',lw=0.2,color='#F56B5C',label=r'$T > 10 \, \mathrm{Myr}$')
            
            #Plot fits
            ls = '--'
            lw = 2.0
            lc = '#F56B5C'
            if(galaxy_function == 'piecewise'):
                ls = '--'
                axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.best_fit_o[0],
                    galaxy_class.best_fit_o[1],galaxy_class.best_fit_o[2],galaxy_class.best_fit_o[3])),
                    ls=ls,color=lc,lw=lw)        
            elif(galaxy_function == 'singlepl'):
                ls = 'dotted'
                axs[i,j].plot(plot_points,np.exp(onepowerlaw_function(plot_points,galaxy_class.best_fit_o[0],
                galaxy_class.best_fit_o[1])),
                ls=ls,lw=lw,color=lc)        
            elif(galaxy_function == 'singletrunc'):
                ls = '-.'
                axs[i,j].plot(plot_points,np.exp(linear_truncation(plot_points,galaxy_class.best_fit_o[0],
                    galaxy_class.best_fit_o[1],galaxy_class.best_fit_o[2])),
                    ls=ls,lw=lw,color=lc)            
                #axs[i,j].axvline(theta_c,ls=':',label=r'$\theta_c = {:2.1f} \pm {:2.1f}$'.format(theta_c,
                #    theta_c_error))
            ls_o = ls
        else:
            fit_label_old = None
            ls_o = None

        
        #Set X-Axis Limits
        distance = galaxy_class.distance*const.Parsec*1.e6
        axs_limits = [0.0,0.0]
        axs_limits[0] =  global_axes_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs_limits[1] =  global_axes_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs[i,j].set_xlim(axs_limits[0],axs_limits[1])

        axs[i,j].set_ylim(0.1,1.2e2)

        #Figure out edge effect region
        region = read_ds9(galaxy_class.region_file)
        hdu = fits.open(galaxy_class.fits_file)[0]
        wcs = WCS(hdu.header)
        name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
        ra_dec = SkyCoord.from_name(name)
        xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)

        for k in range(0,np.size(region)):
            region_dep = deproject_region_centre(region,k,xpix_c,ypix_c,galaxy_class)
            if(region_dep is not None):
                sides = region_dep.to_sky(wcs).vertices
                sizes = get_separations(sides,plot_class)

        ########################################
        #Get probable boundary scale
        boundary_scale = np.max(sizes)/5.0
        boundary_scale = plot_class.pc_to_sep(boundary_scale)

        #Shade region beyond which edge effects play role
        axs[i,j].axvspan(boundary_scale,axs_limits[1],alpha=0.3,color='#8D717490')

        # General Plot stuff
        #X-labels only on bottom row
        if(i==3):
            axs[i,j].set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        #Y-labels only on left column
        if(j == 0):
            axs[i,j].set_ylabel(r"$1+\omega\left(\theta \right)$")
        axs[i,j].set_xscale('log')
        axs[i,j].callbacks.connect("xlim_changed", plot_class.axs_to_parsec)
        
        axs[i,j].tick_params(axis='x',which='both',bottom=True,top=False)

        if(i+j == 0):
            handles, labels = axs[i,j].get_legend_handles_labels()
            handles = [h[0] for h in handles]
            line = [Line2D([0],[0],color='k',lw=lw,ls=ls)]
            label = ['fit']
            #handles = handles+line
            #labels = labels+label

            fit_handle = [Line2D([0],[0],color='#4591F5',lw=lw,ls=ls_y)]
            fit_label = [fit_label_young]
            if(ls_o is not None):
                fit_handle = fit_handle + [Line2D([0],[0],color='#F56B5C',lw=lw,ls=ls_o)]
                fit_label = fit_label + [fit_label_old]

            handles = handles + fit_handle
            labels = labels + fit_label

            age_legend = axs[i,j].legend(handles,labels,loc='upper right',numpoints=1,markerscale=2,
                frameon=True,fancybox=True,framealpha=0.3,facecolor='grey',handletextpad=0.4,fontsize=12)
            axs[i,j].add_artist(age_legend)
        
        handles = [Line2D([0],[0],color='#4591F5',lw=lw,ls=ls_y)]
        labels = [fit_label_young]
        if(ls_o is not None):
            handles = handles + [Line2D([0],[0],color='#F56B5C',lw=lw,ls=ls_o)]
            labels = labels + [fit_label_old]
        if(galaxy_name not in ['NGC_0628']):            
            axs[i,j].legend(handles,labels,loc='upper right',
                frameon=True,fancybox=True,framealpha=0.3,facecolor='grey',handletextpad=0.4,fontsize=12)


        #Secondary axis label only for top row
        if(i==0):
            ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

        
        axs[i,j].text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs[i,j].transAxes)

        if(galaxy_name == "NGC_0628"):
              axs[i,j].text(boundary_scale+10,0.2,"Edge",rotation=90)
              axs[i,j].text(boundary_scale+50,0.2,"Effects",rotation=90)           


        #Get position of subplot
        j +=1
        if(j==3):
            j = 0
            i +=1


    if(save):
        filename = outdir+'Fig2.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()        


def Figure3(save=False,outdir='../results/',indir='../results/'):
    """
    Plot Figure 3 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 3")
    galaxy = "NGC_5194"
    galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
    group_no = 1
    labels = [r'$T<2 \, \mathrm{Myr}$',r'$2<T<10 \, \mathrm{Myr}$',
             r'$10<T<100 \, \mathrm{Myr}$',r'$T>100 \, \mathrm{Myr}$']
    colors = ["#DB0E04","#F5A100","#9D23DB","#324BDB"]

    fig,axs = plt.subplots(ncols=1)
    for group_no in range(1,5):
        i = group_no -1 
        if(group_no == 1):
            corr,dcorr = galaxy_class.g1_corr,galaxy_class.g1_dcorr
            AIC_S,AIC_PF = galaxy_class.g1_aics,galaxy_class.g1_aicpf
            fit_values,fit_errors = galaxy_class.g1_fit_values,galaxy_class.g1_fit_errors
        elif(group_no == 2):
            corr,dcorr = galaxy_class.g2_corr,galaxy_class.g2_dcorr
            AIC_S,AIC_PF = galaxy_class.g2_aics,galaxy_class.g2_aicpf
            fit_values,fit_errors = galaxy_class.g2_fit_values,galaxy_class.g2_fit_errors
        if(group_no == 3):
            corr,dcorr = galaxy_class.g3_corr,galaxy_class.g3_dcorr
            AIC_S,AIC_PF = galaxy_class.g3_aics,galaxy_class.g3_aicpf
            fit_values,fit_errors = galaxy_class.g3_fit_values,galaxy_class.g3_fit_errors
        if(group_no == 4):
            corr,dcorr = galaxy_class.g4_corr,galaxy_class.g4_dcorr
            AIC_S,AIC_PF = galaxy_class.g4_aics,galaxy_class.g4_aicpf
            fit_values,fit_errors = galaxy_class.g4_fit_values,galaxy_class.g4_fit_errors

        galaxy_class.set_bins(set_no_bins=15)
        galaxy_class.corr = corr
        galaxy_class.dcorr = dcorr
        pl = myPlot(galaxy_class)
        separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
        separation_bins*=(1./arcsec_to_degree)
        if(group_no == 1):
            pl.plot_TPCF(function=None,omega1=True,
                         axs=axs,sec_axis=True,fmt='-',lw=0.8,c=colors[i],
                        label=labels[i],age=None)
        else:
            pl.plot_TPCF(function=None,omega1=True,
                         axs=axs,sec_axis=False,fmt='-',lw=0.8,c=colors[i],
                        label=labels[i],age=None)

        plot_points = np.linspace(np.min(separation_bins),
                                  np.max(separation_bins),1000)
        if(AIC_S<AIC_PF):
            axs.plot(plot_points,np.exp(onepowerlaw_function(plot_points,fit_values[0],
                            fit_values[1])),
                            ls='--',c=colors[i],lw=2.0)

        else:
            axs.plot(plot_points,np.exp(linear_truncation(plot_points,fit_values[0],
                            fit_values[1],fit_values[2])),
                            ls='--',c=colors[i],lw=2.0)

    axs.set_ylabel(r'$1+\omega(\theta)$')
    axs.text(0.7,0.55,r'$\alpha_1=-0.55 \, (\mathrm{S})$',transform=axs.transAxes,
             c=colors[0],fontsize=12)
    axs.text(0.7,0.47,r'$\alpha_1=-0.38 \, (\mathrm{S})$',transform=axs.transAxes,
             c=colors[1],fontsize=12)
    axs.text(0.7,0.39,r'$\alpha_1=-0.28 \, (\mathrm{PF})$',transform=axs.transAxes,
             c=colors[2],fontsize=12)
    axs.text(0.7,0.31,r'$\alpha_1=-0.05 \, (\mathrm{PF})$',transform=axs.transAxes,
             c=colors[3],fontsize=12)
    if(save):
        filename = outdir+'Fig3.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()        



def Figure4(save=True,outdir='../results/',indir='../results/'):
    """
    Plot Figure 3 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 4")
    fig,axs = plt.subplots(ncols=1)

    colors = ['#F58300','#F500A8','#006B4A']

    #Draw Model S
    bins = np.logspace(np.log10(1.e-3),np.log10(1.0),100)
    separation_bins = (bins[1:]+bins[:-1])/2.
    #Model S
    axs.plot(separation_bins,7*separation_bins**(-0.5),ls='-',lw=4.0,
             c=colors[0],alpha=0.7)

    #Model PW
    axs.plot(separation_bins,np.exp(linear_function(separation_bins,np.log(1.e-1),
                                                    -1.0,-0.1,np.log(0.05))),
            lw=4.0,c=colors[1],alpha=0.7)

    #Model PF
    axs.plot(separation_bins,np.exp(linear_truncation(separation_bins,0.03,
                                                      -0.1,0.1)),
            lw=4.0,c=colors[2],alpha=0.7)



    #Annotate Model S
    x1,x2 = 0.02,0.1
    y1,y2 = 12*x1**(-0.5),12*x2**(-0.5)
    axs.plot([x1,x2],[y1,y2],ls='--',lw=4.0,c=colors[0],alpha=0.7)
    angle = np.arctan((y2-y1)/(x2-x1))
    axs.text(0.02,1.5e2,r'$\alpha_1 = D_2-2$',c=colors[0],
             size=15.0,rotation=-12.0,rotation_mode='anchor')
    xarrow = np.max(separation_bins)
    yarrow = 12*xarrow**(-0.5)
    axs.arrow(xarrow,yarrow+12,dx=0.0,dy=-12.0,color=colors[0],
                  width=0.04,length_includes_head=False,
                  head_width=0.2,head_length=4.0,lw=1.0,alpha=0.5)
    axs.text(xarrow-0.55,yarrow+17.0,
             r'$\theta_{\mathrm{max}} < l_{\mathrm{corr}}$',
             c=colors[0],size=15.0)

    #Annotate Model PW
    x1,x2 = 2.e-3,1.e-2
    y1,y2 = 0.05*x1**(-1.0),0.05*x2**(-1.0)
    axs.plot([x1,x2],[y1,y2],ls='--',lw=4.0,c=colors[1],alpha=0.7)
    angle = np.arctan((y2-y1)/(x2-x1))
    axs.text(1.5e-3,10.0,r'$\alpha_1 = D_2-2$',c=colors[1],size=15.0,rotation=-24.0,rotation_mode='anchor')

    x1,x2 = 0.2,0.5
    y1,y2 = 2.*x1**(-0.1),2.*x2**(-0.1)
    axs.plot([x1,x2],[y1,y2],ls='--',lw=4.0,c=colors[1],alpha=0.7)
    axs.text(0.19,4.0,r'$\alpha_2 \approx 0$',size=15.0,rotation=-2.0,
             rotation_mode='anchor',c=colors[1])

    axs.arrow(0.05,8.0,dx=0.0,dy=-4.0,color=colors[1],
                  width=0.002,length_includes_head=False,
                  head_width=0.01,head_length=1.3,lw=1.0,alpha=0.5)
    axs.text(0.03,9.7,r'$\beta \approx l_{\mathrm{corr}}$',c=colors[1],size=18.0)

    #Annotate Model PF
    x1,x2 = 2.e-3,9.e-3
    y1,y2 = 0.5*x1**(-0.15),0.5*x2**(-0.15)
    axs.plot([x1,x2],[y1,y2],ls='--',lw=4.0,c=colors[2],alpha=0.7)
    angle = np.arctan((y2-y1)/(x2-x1))
    axs.text(1.5e-3,0.6,r'$\alpha_1 = D_2-2$',c=colors[2],
             size=15.0,rotation=-3.0,rotation_mode='anchor')

    axs.arrow(0.1,0.1,dx=0.0,dy=0.15,color=colors[2],
                  width=0.005,length_includes_head=False,
                  head_width=0.02,head_length=0.12,lw=1.0,alpha=0.5)
    axs.text(0.07,0.05,r'$\theta_c \approx r_c$',size=15.0,c=colors[2])


    #Set Legend
    legend_elements = [Patch(facecolor=colors[0], edgecolor=None,
                             label=r'$\mathrm{Model \; S}$',alpha=0.7),
                       Patch(facecolor=colors[1], edgecolor=None,
                             label=r'$\mathrm{Model \; PW}$',alpha=0.7),
                      Patch(facecolor=colors[2], edgecolor=None,
                             label=r'$\mathrm{Model \; PF}$',alpha=0.7)]
    axs.legend(handles=legend_elements,loc='lower left')

    axs.set_xlim(7.e-4,2.5)
    axs.set_ylim(0.02,5.e2)

    axs.tick_params(axis='both',which='both',length=0.0,labelsize=0.0)
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.set_xlabel(r"$\theta \, (\mathrm{arbitrary \, units})$")
    axs.set_ylabel(r"$1+\omega \, (\mathrm{arbitrary \, units})$")

    if(save):
        filename = outdir+'Fig4.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    return

def Figure5(save=True,outdir='../results/',indir='../results/'):
    """
    Plot Figure 5 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 5")
    all_galaxies = ['NGC_0628','NGC_1313','NGC_1566','NGC_3344','NGC_3627','NGC_4449',                       
                    'NGC_5194','NGC_5253','NGC_5457','NGC_6503','NGC_7793']

    #Read in values required for plotting trends
    l_corr = np.zeros(np.size(all_galaxies))
    dl_corr = np.zeros((2,np.size(all_galaxies)))
    mstar = np.zeros(np.size(all_galaxies))
    sfr_uv = np.zeros(np.size(all_galaxies))
    sigma_sfr = np.zeros(np.size(all_galaxies))
    sigma_star = np.zeros(np.size(all_galaxies))
    R_25 = np.zeros(np.size(all_galaxies))
    Tvalue = np.zeros(np.size(all_galaxies))
    i = 0

    #For galaxies whose FoV is smaller than the extent, we use the Sigma_SFR in the region probed
    #This is stored in the file ../data/galaxyinfo/SFR_Table_extra_short.txt and was 
    # provided by Daniela Calzetti and would be published in Adamo et al 2021, in prep
    newsfr_galaxies = ['NGC_0628','NGC_3627','NGC_5194','NGC_5457']
    file = '../data/galaxyinfo/SFR_Table_extra_short.txt'
    sfr_data = np.loadtxt(file,usecols=10)
    for galaxy in all_galaxies:
        galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
        galaxy_class.getPhysicalProps()

        l_corr[i] = galaxy_class.lcorr[0]
        dl_corr[0,i] = galaxy_class.lcorr[1]
        dl_corr[1,i] = galaxy_class.lcorr[2]
        
        
        galaxy_class.read_galaxyprops()
        mstar[i] = galaxy_class.mstar
        sfr_uv[i] = galaxy_class.sfr
        if(galaxy in newsfr_galaxies):
            index = np.where(galaxy == np.array(newsfr_galaxies))
            sigma_sfr[i] = sfr_data[index]
        else:
            sigma_sfr[i] = galaxy_class.sigma_sfr
        R_25[i] = galaxy_class.r25
        sigma_star[i] = galaxy_class.mstar/(np.pi*galaxy_class.r25**2)
        Tvalue[i] = galaxy_class.T_value
        i = i+1

    #Plot specific settings for markers/colours etc
    lowlims = []
    uplims = []
    for galaxy in all_galaxies:
        if(galaxy in ['NGC_4449','NGC_5253']):
            uplims.append(True)
            lowlims.append(False)
        elif(galaxy in ['NGC_1313','NGC_1566','NGC_3627','NGC_5194']):
            uplims.append(False)
            lowlims.append(True)
        else:
            uplims.append(False)
            lowlims.append(False)
    lowlims = np.array(lowlims)
    uplims = np.array(uplims)

    indices_dwarf = np.where(np.logical_or(np.array(all_galaxies) == 'NGC_4449',
                                 np.array(all_galaxies) == 'NGC_5253'))[0]
    indices_singlepl = np.where(lowlims)
    indices_pwpl = np.where(~np.logical_or(lowlims,uplims))
    indices_spiral = np.delete(np.arange(0,np.size(all_galaxies),1),indices_dwarf)

    mfc = ['#6301DB']*np.size(all_galaxies)
    mec = mfc
    marker = ['s']*np.size(all_galaxies)
    marker = np.array(marker)
    mfc = np.array(mfc)
    mec = np.array(mec)
    mfc[indices_singlepl] = 'white'
    mfc[indices_dwarf] = 'white'
    marker[indices_dwarf] = 'o'
    alpha = [1.0]*np.size(all_galaxies)
    alpha = np.array(alpha)
    alpha[indices_pwpl] = 0.7


    #Get correlation coefficients for each pair of variables
    corr = np.zeros(6)
    pvalues = np.zeros(6)
    corr[0],pvalues[0] = scipy.stats.pearsonr(R_25,l_corr)
    corr[1],pvalues[1] = scipy.stats.pearsonr(mstar,l_corr)
    corr[2],pvalues[2] = scipy.stats.pearsonr(sfr_uv,l_corr)
    corr[3],pvalues[3] = scipy.stats.pearsonr(Tvalue,l_corr)
    corr[4],pvalues[4] = scipy.stats.pearsonr(sigma_star,l_corr)
    corr[5],pvalues[5] = scipy.stats.pearsonr(sigma_sfr,l_corr)

    #Plot trends
    fig,axs = plt.subplots(ncols=3,nrows=2,figsize=(14,8),sharey='row')
    for i in range(np.size(all_galaxies)):
        
        #For lower/upper limits set representative arrow size
        yerr = np.array(dl_corr[:,i]).reshape(2,1)
        if(mfc[i] == 'white'):
            yerr = 0.3*l_corr[i]
            
        axs[0,0].errorbar(R_25[i],l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])    
        
        axs[0,1].errorbar(mstar[i]/1.e10,l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        axs[0,2].errorbar(sfr_uv[i],l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        axs[1,0].errorbar(Tvalue[i],l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])    
        
        
        axs[1,1].errorbar(sigma_star[i]/1.e7,l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        axs[1,2].errorbar(sigma_sfr[i]*1.e3,l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        
    axs[0,0].set_xlabel(r'$R_{25} \, (\mathrm{kpc})$')    
    axs[0,0].set_ylabel(r'$l_{\mathrm{corr}} \, \left(\mathrm{pc} \right)$')
    axs[0,1].set_xlabel(r'$M_* \, \left( 10^{10} M_{\odot} \right)$')
    axs[0,2].set_xlabel(r'$\mathrm{SFR}_{\mathrm{UV}} \, \left(\mathrm M_\odot \, \mathrm{yr}^{-1} \right)$')
    axs[1,0].set_xlabel(r'$\mathrm{Morphological} \; T$')
    axs[1,2].set_xlabel(r'$\Sigma_\mathrm{SFR} \, \left(10^{-3} \, \mathrm M_\odot \, \mathrm{yr}^{-1} \, \mathrm{kpc}^{-2} \, \right)$')
    axs[1,0].set_ylabel(r'$l_{\mathrm{corr}} \, \left(\mathrm{pc} \right)$')
    axs[1,1].set_xlabel(r'$\Sigma_{*} \, \left(10^7 \mathrm M_\odot \, \mathrm{kpc}^{-2} \, \right)$')
    
    #Remove ticks on right side to avoid confusion
    axs[0,0].tick_params(axis='y',which='both',left=True,right=False)
    axs[0,1].tick_params(axis='y',which='both',left=True,right=False)
    axs[1,0].tick_params(axis='y',which='both',left=True,right=False)
    axs[1,1].tick_params(axis='y',which='both',left=True,right=False)


    for i in range(0,2):
        for j in range(0,3):
            axs[i,j].set_yscale('log')
            xlim = axs[i,j].get_xlim()
            xlen = (xlim[1]-xlim[0])
            axs[i,j].set_xlim(xlim[0]-xlen*0.1,xlim[1]+xlen*0.1)
            
            ylim = axs[i,j].get_ylim()
            ylen = (ylim[1]-ylim[0])
            axs[i,j].set_ylim(ylim[0]-ylen*0.1,ylim[1]+ylen*0.1)
            index = int(3*i+j)
            axs[i,j].text(0.1,0.9,
                 r'$\rho = \,{:2.2f} \; (p = {:2.2f})$'.format(corr[index],pvalues[index]),
                 transform=axs[i,j].transAxes,fontsize=16)
    import matplotlib.lines as mlines
    #Make Legend
    fill_square = mlines.Line2D([],[],mfc=mfc[0],marker='s',
                                label='Spirals (PW)',ls='',mec=mec[0],
                               mew=1.0,alpha=0.7)
    open_square = mlines.Line2D([],[],mfc='white',marker='s',
                                label='Spirals (S)', ls='',mec=mec[0],
                               mew=1.0)
    open_circle = mlines.Line2D([],[],mfc='white',marker='o',
                                label='Dwarfs (PF)', ls='',mec=mec[0],
                               mew=1.0)
    axs[0,0].legend(handles=[fill_square,open_square,open_circle],
                    loc='lower right',numpoints=1,markerscale=1.0,
                   frameon=True,fancybox=True,framealpha=0.15,facecolor='grey')
            
    plt.subplots_adjust(wspace=0.0,hspace=0.25)

    if(save):
        filename = outdir+'Fig5.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

#Functions required for Figure 6
def toomre_length(sigma_g,v_rot,r):
    omega = v_rot/r
    kappa = np.sqrt(2)*omega
    value = 4*np.pi**2*const.G_GravityConstant*sigma_g/(kappa**2)
    return value

def get_distances(galaxy_class,age='young'):
    from TPCF import ra_dec_to_xyz
    name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]

    #Center coordinates
    ra_dec = SkyCoord.from_name(name)
    ra_centre,dec_centre = ra_dec.ra.value,ra_dec.dec.value
    #Deproject FOV to prepare mask
    hdu = fits.open(galaxy_class.fits_file)[0]
    wcs = WCS(hdu.header)
    xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)

    galaxy_class.get_ra_dec(age=age)
    galaxy_ra,galaxy_dec = galaxy_class.ra,galaxy_class.dec
    
    xyz_centre = np.array(ra_dec_to_xyz(ra_centre,dec_centre))
    xyz_data = np.array(ra_dec_to_xyz(galaxy_ra,galaxy_dec))
    
    # Radial distances from centre in deg
    distance_data = np.sqrt((xyz_data[0]-xyz_centre[0])**2 + (xyz_data[1]-\
            xyz_centre[1])**2 + (xyz_data[2]-xyz_centre[2])**2)
    distance = distance_data*galaxy_class.distance*1.e6*const.Parsec/(1.e3*const.Parsec)
    return distance
    

def Figure6(save=True,outdir='../results/',indir='../results/'):
    """
    Plot Figure 6 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 6")
    spiral_galaxies = ['NGC_0628','NGC_1313','NGC_1566','NGC_3344','NGC_3627', 'NGC_5194',                      
        'NGC_5457','NGC_6503','NGC_7793']
    v_rot = [217.0,220.0,161.0,158.0,192.0,219.0,180.0,110.0,125.0]
    sigma_hi = np.zeros(np.size(spiral_galaxies))
    l_toomre = np.zeros(np.size(spiral_galaxies))
    l_corr = np.zeros(np.size(spiral_galaxies))
    dl_corr = np.zeros((2,np.size(spiral_galaxies)))
    sigma_gas = loadObj('../data/galaxyinfo/AveragedSigmaH2')
    sigma_gas = np.delete(sigma_gas,[5,6,8])

    #Radius to calculate
    r = np.zeros(np.size(spiral_galaxies))

    i = 0
    for galaxy in spiral_galaxies:
        galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
        galaxy_class.getPhysicalProps()

        l_corr[i] = galaxy_class.lcorr[0]
        dl_corr[0,i] = galaxy_class.lcorr[1]
        dl_corr[1,i] = galaxy_class.lcorr[2]

        galaxy_class.read_galaxyprops()
        sigma_hi[i] = galaxy_class.sigma_hi
        
        #Radius calculation
        distances = get_distances(galaxy_class)
        r[i] = np.median(distances)
        
        
        i = i+1
    sigma_gas[1] = sigma_hi[1]
    l_toomre = toomre_length(sigma_gas * const.SolarMass/(const.Parsec**2),
                                   np.array(v_rot)*1.e5,r = r*1000.*const.Parsec)/const.Parsec

    #Settings for next plot
    lowlims = []
    uplims = []
    for galaxy in spiral_galaxies:
        if(galaxy in ['NGC_4449','NGC_5253']):
            uplims.append(True)
            lowlims.append(False)
        elif(galaxy in ['NGC_1313','NGC_1566','NGC_3627','NGC_5194']):
            uplims.append(False)
            lowlims.append(True)
        else:
            uplims.append(False)
            lowlims.append(False)
    lowlims = np.array(lowlims) 



    #Settings for the following two plots
    indices_dwarf = np.where(np.logical_or(np.array(spiral_galaxies) == 'NGC_4449',
                                 np.array(spiral_galaxies) == 'NGC_5253'))[0]
    indices_singlepl = np.where(lowlims)
    indices_pwpl = np.where(~np.logical_or(lowlims,uplims))
    indices_spiral = np.delete(np.arange(0,np.size(spiral_galaxies),1),indices_dwarf)


    #mfc = ['#5C5850']*np.size(spiral_galaxies)
    mfc = ['#6301DB']*np.size(spiral_galaxies)

    mec = mfc
    marker = ['s']*np.size(spiral_galaxies)
    marker = np.array(marker)
    mfc = np.array(mfc)
    mec = np.array(mec)
    mfc[indices_singlepl] = 'white'
    mfc[indices_dwarf] = 'white'
    #mfc[indices_dwarf] = '#6301DB'
    #mec[indices_dwarf] = '#6301DB'
    marker[indices_dwarf] = 'o'
    alpha = [1.0]*np.size(spiral_galaxies)
    alpha = np.array(alpha)
    alpha[indices_pwpl] = 0.7

    fig,axs = plt.subplots(ncols=1)
    for i in range(np.size(spiral_galaxies)):
        
        yerr = np.array(dl_corr[:,i]).reshape(2,1)
        if(mfc[i] == 'white'):
            yerr = 0.3*l_corr[i]
        
        axs.errorbar(l_toomre[i],l_corr[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
    corr,pvalues = scipy.stats.spearmanr(l_corr,l_toomre).correlation, scipy.stats.spearmanr(l_corr,l_toomre).pvalue     
    axs.text(0.1,0.9,
                 r'$\rho = \,{:2.2f} \; (p = {:2.2f})$'.format(corr,pvalues),
                 transform=axs.transAxes,fontsize=16)

    axs.set_xlabel(r'$l_{\mathrm{T}} \, \left(\mathrm{pc} \right)$')
    axs.set_ylabel(r'$l_{\mathrm{corr}} \, \left(\mathrm{pc} \right)$')

    axs.set_xscale('log')
    axs.set_yscale('log')

    xlim = axs.get_xlim()
    xlen = (xlim[1]-xlim[0])
    axs.set_xlim(xlim[0]-xlen*0.1,xlim[1]+xlen*0.7)

    ylim = axs.get_ylim()
    ylen = (ylim[1]-ylim[0])
    axs.set_ylim(ylim[0]-ylen*0.1,ylim[1]+ylen*0.1)



    import matplotlib.lines as mlines
    #Make Legend
    fill_square = mlines.Line2D([],[],mfc=mfc[0],marker='s',
                                label='Spirals (PW)',ls='',mec=mec[0],
                               mew=1.0,alpha=0.7)
    open_square = mlines.Line2D([],[],mfc='white',marker='s',
                                label='Spirals (S)', ls='',mec=mec[0],
                               mew=1.0)

    axs.legend(handles=[fill_square,open_square],
                    loc='lower right',numpoints=1,markerscale=1.0,
                   frameon=True,fancybox=True,framealpha=0.15,facecolor='grey')

    if(save):
        filename = outdir+'Fig6.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()    
        

def Figure7(save=True,outdir='../results/',indir='../results/'):
    """
    Plot Figure 7 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 7")
    all_galaxies = ['NGC_0628','NGC_1313','NGC_1566','NGC_3344','NGC_3627','NGC_4449',                       
                    'NGC_5194','NGC_5253','NGC_5457','NGC_6503','NGC_7793']

    #Read in values required for plotting trends
    D_2 = np.zeros(np.size(all_galaxies))
    D_2_error = np.zeros((2,np.size(all_galaxies)))
    mstar = np.zeros(np.size(all_galaxies))
    sfr_uv = np.zeros(np.size(all_galaxies))
    sigma_sfr = np.zeros(np.size(all_galaxies))
    sigma_star = np.zeros(np.size(all_galaxies))
    R_25 = np.zeros(np.size(all_galaxies))
    Tvalue = np.zeros(np.size(all_galaxies))
    i = 0
    newsfr_galaxies = ['NGC_0628','NGC_3627','NGC_5194','NGC_5457']
    file = '../data/galaxyinfo/SFR_Table_extra_short.txt'
    sfr_data = np.loadtxt(file,usecols=10)
    for galaxy in all_galaxies:
        galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
        galaxy_class.getPhysicalProps()
        D_2[i] = galaxy_class.D2[0]
        D_2_error[1,i] = galaxy_class.D2[1]
        D_2_error[0,i] = galaxy_class.D2[2]
        
        galaxy_class.read_galaxyprops()
        mstar[i] = galaxy_class.mstar
        sfr_uv[i] = galaxy_class.sfr
        if(galaxy in newsfr_galaxies):
            index = np.where(galaxy == np.array(newsfr_galaxies))
            sigma_sfr[i] = sfr_data[index]
        else:
            sigma_sfr[i] = galaxy_class.sigma_sfr
        R_25[i] = galaxy_class.r25
        sigma_star[i] = galaxy_class.mstar/(np.pi*galaxy_class.r25**2)
        Tvalue[i] = galaxy_class.T_value
        i = i+1

    #Plot specific settings for markers/colours etc
    lowlims = []
    uplims = []
    for galaxy in all_galaxies:
        if(galaxy in ['NGC_4449','NGC_5253']):
            uplims.append(True)
            lowlims.append(False)
        elif(galaxy in ['NGC_1313','NGC_1566','NGC_3627','NGC_5194']):
            uplims.append(False)
            lowlims.append(True)
        else:
            uplims.append(False)
            lowlims.append(False)
    lowlims = np.array(lowlims)
    uplims = np.array(uplims)

    indices_dwarf = np.where(np.logical_or(np.array(all_galaxies) == 'NGC_4449',
                                 np.array(all_galaxies) == 'NGC_5253'))[0]
    indices_singlepl = np.where(lowlims)
    indices_pwpl = np.where(~np.logical_or(lowlims,uplims))
    indices_spiral = np.delete(np.arange(0,np.size(all_galaxies),1),indices_dwarf)

    mfc = ['#6301DB']*np.size(all_galaxies)
    mec = mfc
    marker = ['s']*np.size(all_galaxies)
    marker = np.array(marker)
    mfc = np.array(mfc)
    mec = np.array(mec)
    mfc[indices_singlepl] = 'white'
    mfc[indices_dwarf] = 'white'
    marker[indices_dwarf] = 'o'
    alpha = [1.0]*np.size(all_galaxies)
    alpha = np.array(alpha)
    alpha[indices_pwpl] = 0.7


    #Get correlation coefficients for each pair of variables
    corr = np.zeros(6)
    pvalues = np.zeros(6)
    corr[0],pvalues[0] = scipy.stats.pearsonr(R_25,D_2)
    corr[1],pvalues[1] = scipy.stats.pearsonr(mstar,D_2)
    corr[2],pvalues[2] = scipy.stats.pearsonr(sfr_uv,D_2)
    corr[3],pvalues[3] = scipy.stats.pearsonr(Tvalue,D_2)
    corr[4],pvalues[4] = scipy.stats.pearsonr(sigma_star,D_2)
    corr[5],pvalues[5] = scipy.stats.pearsonr(sigma_sfr,D_2)

    #Plot trends
    fig,axs = plt.subplots(ncols=3,nrows=2,figsize=(14,8),sharey='row')
    for i in range(np.size(all_galaxies)):
        
        #For lower/upper limits set representative arrow size
        yerr = np.array(D_2_error[:,i]).reshape(2,1)
            
        axs[0,0].errorbar(R_25[i],D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])    
        
        axs[0,1].errorbar(mstar[i]/1.e10,D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        axs[0,2].errorbar(sfr_uv[i],D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        axs[1,0].errorbar(Tvalue[i],D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])    
        
        
        axs[1,1].errorbar(sigma_star[i]/1.e7,D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        axs[1,2].errorbar(sigma_sfr[i]*1.e3,D_2[i],
                      marker=marker[i],mfc=mfc[i],mec=mec[i],ls='',
                      yerr=yerr,
                      capsize=2.0,
                      lolims=lowlims[i],uplims=uplims[i],ms=6.0,ecolor='k',
                         mew=1.0,alpha=alpha[i])
        
        
    axs[0,0].set_xlabel(r'$R_{25} \, (\mathrm{kpc})$')    
    axs[0,0].set_ylabel(r'$D_2$')
    axs[0,1].set_xlabel(r'$M_* \, \left( 10^{10} M_{\odot} \right)$')
    axs[0,2].set_xlabel(r'$\mathrm{SFR}_{\mathrm{UV}} \, \left(\mathrm M_\odot \, \mathrm{yr}^{-1} \right)$')
    axs[1,0].set_xlabel(r'$\mathrm{Morphological} \; T$')
    axs[1,2].set_xlabel(r'$\Sigma_\mathrm{SFR} \, \left(10^{-3} \, \mathrm M_\odot \, \mathrm{yr}^{-1} \, \mathrm{kpc}^{-2} \, \right)$')
    axs[1,0].set_ylabel(r'$D_2$')
    axs[1,1].set_xlabel(r'$\Sigma_{*} \, \left(10^7 \mathrm M_\odot \, \mathrm{kpc}^{-2} \, \right)$')
    
    #Remove ticks on right side to avoid confusion
    axs[0,0].tick_params(axis='y',which='both',left=True,right=False)
    axs[0,1].tick_params(axis='y',which='both',left=True,right=False)
    axs[1,0].tick_params(axis='y',which='both',left=True,right=False)
    axs[1,1].tick_params(axis='y',which='both',left=True,right=False)


    for i in range(0,2):
        for j in range(0,3):
            xlim = axs[i,j].get_xlim()
            xlen = (xlim[1]-xlim[0])
            axs[i,j].set_xlim(xlim[0]-xlen*0.1,xlim[1]+xlen*0.1)
            
            ylim = axs[i,j].get_ylim()
            ylen = (ylim[1]-ylim[0])
            axs[i,j].set_ylim(ylim[0]-ylen*0.1,ylim[1]+ylen*0.1)
            index = int(3*i+j)
            axs[i,j].text(0.1,0.9,
                 r'$\rho = \,{:2.2f} \; (p = {:2.2f})$'.format(corr[index],pvalues[index]),
                 transform=axs[i,j].transAxes,fontsize=16)
    import matplotlib.lines as mlines
    #Make Legend
    fill_square = mlines.Line2D([],[],mfc=mfc[0],marker='s',
                                label='Spirals (PW)',ls='',mec=mec[0],
                               mew=1.0,alpha=0.7)
    open_square = mlines.Line2D([],[],mfc='white',marker='s',
                                label='Spirals (S)', ls='',mec=mec[0],
                               mew=1.0)
    open_circle = mlines.Line2D([],[],mfc='white',marker='o',
                                label='Dwarfs (PF)', ls='',mec=mec[0],
                               mew=1.0)
    axs[0,0].legend(handles=[fill_square,open_square,open_circle],
                    loc='lower right',numpoints=1,markerscale=1.0,
                   frameon=True,fancybox=True,framealpha=0.15,facecolor='grey')
            
    plt.subplots_adjust(wspace=0.0,hspace=0.25)

    if(save):
        filename = outdir+'Fig7.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def Figure8(save=True,outdir='../results/',indir='../results/'):
    """
    Plot Figure 8 of Menon et al 2021b
    Parameters
        save: Boolean   
            Flag to save the plot
        outdir: string
            Output directory in which to store plot. Default is results directory.
        indir : string
            Input directory from which to read the results. Default is results directory.
    Returns:
        None 
    """
    print("Creating Figure 8")
    r_c = np.zeros(np.size(list_of_galaxies))
    r_c_e = np.zeros((2,np.size(list_of_galaxies)))
    i = 0
    for galaxy in list_of_galaxies:
        galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
        galaxy_class.getPhysicalProps()
        rc = galaxy_class.rc_disk
        if(rc is not None):
            r_c[i], r_c_e[0,i], r_c_e[1,i] = rc  
        else:
            r_c[i], r_c_e[0,i], r_c_e[1,i] = 0,0,0
        i = i+1

    #Read Spitzer survey scale radii        
    file = '../data/galaxyinfo/Scale_Radii.dat'
    data = np.loadtxt(file,usecols=6)
    rc_spitzer = np.zeros(np.size(list_of_galaxies))
    rc_spitzererr = np.zeros((2,np.size(list_of_galaxies)))
    i = 0
    for galaxy in list_of_galaxies:
        galaxy_class = loadObj(indir+'{}_summary'.format(galaxy))
        pl = myPlot(galaxy_class)

        distance = galaxy_class.distance*const.Parsec*1.e6
        distance_err = galaxy_class.errordist*const.Parsec*1.e6
        distance_samples = np.random.normal(distance,distance_err,size=10000)
        arcsec_to_pc = u.arcsec.to(u.radian)*distance_samples/(const.Parsec)
        
        rc_spitzer[i] = np.median(data[i]*arcsec_to_pc)
        rc_spitzererr[0,i] = rc_spitzer[i] - np.percentile(data[i]*arcsec_to_pc,16)
        rc_spitzererr[1,i] = np.percentile(data[i]*arcsec_to_pc,84) - rc_spitzer[i]
        i = i+1

    #Filter out NGC 3344 and NGC 7793 since we don't have r_c estimates for them
    indices = np.where(~np.logical_or(np.array(list_of_galaxies) == 'NGC_3344', 
                                 np.array(list_of_galaxies) == 'NGC_7793'))

    fig,axs = plt.subplots(ncols=1)
    axs.errorbar(rc_spitzer[indices]/1000.,
                 r_c[indices]/1000.,marker='s',ls='',yerr=r_c_e[:,indices].reshape(2,10)/1000.,
                c='#D95B55',capsize=5.0,xerr=rc_spitzererr[:,indices].reshape(2,10)/1000.)
    axs.axline((0,0),slope=1.0,ls='--',c='#7238FF',lw=3.0,alpha=0.7)
    axs.set_ylim(-4.0,25.0)
    axs.set_ylabel(r'$r_c \, \left( \mathrm{kpc} \right)$')
    axs.set_xlabel(r'$r_{\mathrm{spitzer}} \, \left( \mathrm{kpc} \right)$')
    for i,txt in enumerate(list_of_galaxies):
        text = txt.split('_')[1]
        if(r_c[i]>0.0):
            if(txt == 'NGC_5253'):
                axs.annotate(text,(rc_spitzer[i]/1000.,r_c[i]/1000.),xytext=(-12,-12),
                    textcoords='offset pixels')
            elif(txt == 'NGC_3738'):
                axs.annotate(text,(rc_spitzer[i]/1000.,r_c[i]/1000.),xytext=(-5,10),
                    textcoords='offset pixels')
            else:
                axs.annotate(text,(rc_spitzer[i]/1000.,r_c[i]/1000.),xytext=(5,5),
                    textcoords='offset pixels')

    if(save):
        filename = outdir+'Fig8.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()


#Appendix Plots
def calzetti_pure(r,A,alpha):
    return A*(r**(alpha)) -  1.0


def FigA1(save=True,outdir='../Menon21_Plots/'):

    fig,axs = plt.subplots(ncols=1)
    fractal_dim = 1.5
    sizes = [0.05,0.1,0.2]
    directory = '/Users/shm/Desktop/WORK/Cluster_2ptCorrelation/Toy_Models/Fractals/\
    Gadi/Monte_Carlo_Boundaries/Fractal_{}'.format(fractal_dim)
    directory_full = '/Users/shm/Desktop/WORK/Cluster_2ptCorrelation/Toy_Models/Fractals/\
    Gadi/'
    directory = os.path.abspath(directory) + '/'
    directory_full = os.path.abspath(directory_full) + '/'
    colors = ['#FB005A','#F29A0C','#6A00DB']
    i = 0
    for size in sizes:
        bins,corr,dcorr = loadObj(directory+'Size_{}'.format(size))
        separation_bins = (bins[1:]+bins[:-1])/2.
        indices = np.where(corr>0.0)
        dcorr_lz = dcorr[indices]
        separation_bins = separation_bins[indices]
        corr_lz = corr[indices]
        bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
        plot_bins = np.linspace(np.min(separation_bins),np.max(separation_bins),200)

        popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/corr_lz)
        perr= np.sqrt(np.diag(pcov))

        r_max = size/2.
        axs.errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                     fmt='.-',lw=0.2,c=colors[i],label=r'$R_{\mathrm{max}} =$'+r'${}$'.format(r_max))
        
        f = calzetti_pure(plot_bins,np.exp(popt[0]),popt[1])
        axs.plot(plot_bins,f,'-',c='k')

        
        r_edge = r_max/5.

        f = calzetti_pure(r_edge,np.exp(popt[0]),popt[1])

        i = i+1

    r_edge= sizes[0]/10.
    axs.arrow(r_edge,4.e-2,dx=0.0,dy=0.13,color=colors[0],
                  width=0.0000005,length_includes_head=True,
                  head_width=0.0002,head_length=0.01,lw=4.0,alpha=0.5)

    r_edge= sizes[1]/10.
    axs.arrow(r_edge,4.e-2,dx=0.0,dy=0.12,color=colors[1],
                  width=0.000001,length_includes_head=True,
                  head_width=0.0004,head_length=0.01,lw=4.0,alpha=0.5)
        
        
    r_edge = sizes[2]/10.
    axs.arrow(r_edge,4.e-2,dx=0.0,dy=0.12,color=colors[2],
                  width=0.000001,length_includes_head=True,
                  head_width=0.001,head_length=0.01,lw=4.0,alpha=0.5,
                      )


        
    axs.set_xscale('log')
    axs.set_yscale('log')    
    axs.set_xlim(9.e-5,0.1)
    axs.set_ylim(1.e-3,12.0)
    axs.legend(loc='lower left')
    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$\omega$")

    if(save):
        filename = outdir+'FigA1.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def FigC1(save=True,outdir='../Menon21_Plots/'):
    fig,axs = plt.subplots(ncols=2,figsize=(12,4))

    directory = '../Toy_Models/Fractals/Gadi/'
    fractal_dims = [0.8,1.0,1.3,1.5,1.7,2.0]
    i = 0
    colors=["#D9AB38","#9D23DB","#324BDB","#7BDB4F","#DB0E04","k"]
    for D_in in fractal_dims:
        color = colors[i]
        label = "{}".format(D_in)
        if(D_in == 2.0):
            xf,yf = np.random.uniform(0,1,size=5000),np.random.uniform(0,1,size=5000)
            bins,corr,dcorr = compute_tpcf(xf,yf)
            dcorr = dcorr/5.
        else:
            xf,yf,bins,corr,dcorr = loadObj(directory+'D_{}'.format(int(D_in*10.0)))
            xf,yf = cull_points(xf,yf,5000)
            bins,corr,dcorr = compute_tpcf(xf,yf)
        separation_bins,corr_lz,dcorr_lz = filter_stuff(bins,corr,dcorr)
        axs[0].errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                     fmt='-',lw=3.0,label=label,color=color,elinewidth=1.0,errorevery=2)
        i +=1
        
    axs[0].set_xlim(9.e-4,1.0)
    axs_limits = axs[0].get_xlim()
    axs[0].axvspan(1./8,axs_limits[1],alpha=0.2,color='#7BB9E3')
    axs[0].set_xlabel(r"$\Delta x$")
    axs[0].set_ylabel(r"$1+\omega$")
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].legend(loc='upper right')


    fractal_dims = [0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.0]
    D_out = np.zeros(np.size(fractal_dims))
    D_out_err = np.zeros(np.size(fractal_dims))

    i = 0
    for D_in in fractal_dims:
        if(D_in == 2.0):
            xf,yf = np.random.uniform(0,1,size=5000),np.random.uniform(0,1,size=5000)
            D_out[i] = 0.0
            D_out_err[i] = 0.01
            continue
        else:
            xf,yf,bins,corr,dcorr = loadObj(directory+'D_{}'.format(int(D_in*10.0)))
            xf,yf = cull_points(xf,yf,5000)
        bins,corr,dcorr = compute_tpcf(xf,yf)
        popt,perr = fit_tpcf(bins,corr,dcorr,fit='piecewise')
        D_out[i]= popt[1]
        D_out_err[i]= perr[1]
        i = i+1
        
    D_out[-1] = 0.0
    axs[1].errorbar(fractal_dims,2+D_out,yerr=D_out_err,fmt='o',ms=8.0)
    axs[1].plot(fractal_dims,fractal_dims,ls='--',lw=3.0)
    axs[1].set_xlim(0.6,2.1)
    axs[1].set_ylim(0.6,2.1)
    axs[1].set_xlabel(r'$D_{2}$')
    axs[1].set_ylabel(r'$\alpha +2$')

    if(save):
        filename = outdir+'FigC1.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def FigC2(save=True,outdir='../Menon21_Plots'):
    fig,axs = plt.subplots(ncols=1)
    lmaxs = [64,32,16,8]
    colors=["#D9AB38","#9D23DB","#324BDB","#DB0E04"]
    labels = [r'$L_{\mathrm{max}} = 1/64$',
              r'$L_{\mathrm{max}} = 1/32$',r'$L_{\mathrm{max}} = 1/16$',
             r'$L_{\mathrm{max}} = 1/8$']
    i = 0
    for lmax in lmaxs:
        x,y = loadObj('../Toy_Models/Fractals/Fractal_Lmax{}'.format(int(lmax)))
        x,y = cull_points(x,y)
        bins,corr,dcorr = compute_tpcf(x,y,no_bins=51)
        separation_bins,corr_lz,dcorr_lz = filter_stuff(bins,corr,dcorr)
        axs.errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                         fmt='-',lw=3.0,label=labels[i],color=colors[i],
                         elinewidth=1.0,errorevery=2,alpha=0.6)
        popt,perr = fit_tpcf(bins,corr,dcorr,fit='piecewise')
        axs.axvline(np.exp(popt[3]),ls='--',lw=2.0,color=colors[i],alpha=0.5)
        i = i+1
        
    #Arrows for max hierarchy length scale
    axs.arrow(1./64.,0.3,dx=0.0,dy=0.6,color=colors[0],
                  width=0.001,length_includes_head=True,
                  head_width=0.003,head_length=0.15,lw=2.0,alpha=0.7)

    axs.arrow(1./32.,0.3,dx=0.0,dy=0.6,color=colors[1],
                  width=0.002,length_includes_head=True,
                  head_width=0.006,head_length=0.15,lw=2.0,alpha=0.7)

    axs.arrow(1./16.,0.3,dx=0.0,dy=0.6,color=colors[2],
                  width=0.005,length_includes_head=True,
                  head_width=0.015,head_length=0.15,lw=2.0,alpha=0.7)

    axs.arrow(1./8.,0.3,dx=0.0,dy=0.6,color=colors[3],
                  width=0.008,length_includes_head=True,
                  head_width=0.024,head_length=0.15,lw=2.0,alpha=0.7)

    axs.set_xlim(1.e-3,1.2)    
    axs.set_ylim(1.e-1,1.e3)
    axs.set_xscale('log')
    axs.set_yscale('log')
    handles,legends = axs.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    line = [Line2D([0],[0],color='k',lw=2.0,ls='--')]
    label = [r'$\beta \, (\mathrm{Model \, PW})$']
    handles = handles+line
    labels = labels+label
    axs.legend(handles,labels,loc='best',fontsize=10)
    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$1+ \omega \left(\theta \right)$")

    if(save):
        filename = outdir+'FigC2.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def FigC3(save=True,outdir='../Menon21_Plots/'):
    fig,axs = plt.subplots(ncols=2,figsize=(12,4))
    array_rc = [0.05,0.1,0.2,0.4]

    R_max = 1.0
    r_min = 0.01
    size = 5000
    inclin = 0.
    z_h = 0.05

    #Exp disk with spiral arms
    file = '../Toy_Models/Bruce_Spirals/tpcf_spiral.xyzr'
    data =  np.loadtxt(file)
    x,y,z = data[:,0],data[:,1],data[:,2]
    x_m,y_m = x[np.where(np.sqrt(x**2+y**2)>0.01)],y[np.where(np.sqrt(x**2+y**2)>0.01)]
    bins_s,corr_s,dcorr_s = compute_tpcf(x_m,y_m,R_min=1.e-2,no_bins=31)

    colors=["#D9AB38","#025AF6","#ED00F5","#DB0E04"]
    for i in tqdm.trange(np.size(array_rc)):
        color=colors[i]
        r_c = array_rc[i]
        label = r"$r_c={}$".format(r_c)
        x,y = Create_Galaxy(r_c=r_c,z_h=z_h,i=inclin,r_min=r_min)
        bins,corr,dcorr = compute_tpcf(x,y,R_min=1.e-2,R_max=min(4.0*r_c,1.0))
        separation_bins,corr_lz,dcorr_lz = filter_stuff(bins,corr,dcorr)
        axs[0].errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                     fmt='.-',lw=1.0,label=label,color=color,elinewidth=1.0,errorevery=5,alpha=0.8)
    separation_bins,corr_s,dcorr_s = filter_stuff(bins_s,corr_s,dcorr_s)
    axs[0].plot(separation_bins,1+corr_s,
                     ls='--',lw=3.0,label=r'$r_c=0.2 \, (\mathrm{Spiral)}$',color='#ED00F5',zorder=10)
    axs[0].set_xlim(5.e-3,1.0)
    axs[0].set_ylabel(r"$1+\omega$")
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].legend(loc='lower left')
    axs[0].set_xlabel(r"$\Delta x$")

    r_cs = np.logspace(np.log10(0.01),np.log10(0.5),10)
    R_max = 1.0
    size = 8000
    z_h = 0.05
    inclin = 0.0
    r_min = 0.0

    rc_fits = np.zeros(np.size(r_cs))
    rc_fitserr = np.zeros(np.size(r_cs))
    i = 0
    for r_c in r_cs:
        x,y = Create_Galaxy(r_c=r_c,z_h=z_h,i=inclin,R_max=R_max,r_min=r_min)
        bins,corr,dcorr = compute_tpcf(x,y,R_max=r_c)
        popt,perr = fit_tpcf(bins,corr,dcorr,fit='singletrunc')
        rc_fits[i],rc_fitserr[i] = popt[2],perr[2]
        i = i+1
        
    axs[1].errorbar(r_cs,rc_fits,yerr=rc_fitserr,fmt='o',ms=8.0,zorder=3)
    axs[1].axline((r_cs[0],r_cs[0]),(r_cs[-1],r_cs[-1]),ls='--',lw=3.0,c='green')
    #axs.plot(r_cs,r_cs,ls='--',lw=3.0)
    axs[1].set_xlim(0.007,1.0)
    #axs.set_ylim(0.6,2.1)
    axs[1].set_xlabel(r'$r_c$')
    axs[1].set_ylabel(r'$\theta_c$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_ylim(4.e-3,1.0)
    axs[1].set_xlim(4.e-3,1.0)

    if(save):
        filename = outdir+'FigC3.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def FigC4(save=True,outdir='../Menon21_Plots/') :

    r_c = 0.2
    R_max = 1.0
    i = 30.0
    r_min = 0.05

    fig,axs = plt.subplots(ncols=1)


    x,y = Create_Galaxy(r_c=r_c,R_max=2.0*r_c,i=i,r_min=r_min)
    bins,corr,dcorr = compute_tpcf(x,y,R_max=r_c,R_min=1.e-2,no_bins=30)
    plot_tpcf(bins,corr,dcorr,axs=axs,fit=None,
              label=r'$R_{\mathrm{max}}=2r_c$',c='#FB005A',fmt='.-',lw=1.0,
              elinewidth=1.0,errorevery=5)
    bins_best ,corr_best ,dcorr_best = bins,corr,dcorr
    popt,perr = fit_tpcf(bins,corr,dcorr,fit='singletrunc')

    fit_best,fiterr_best = popt,perr
    axs.axvline(popt[2],ls='--',c='#FB005A',lw=2.0)

    x,y = Create_Galaxy(r_c=r_c,R_max=3.0*r_c,i=i,r_min=r_min)
    bins,corr,dcorr = compute_tpcf(x,y,R_max=r_c,R_min=1.e-2,no_bins=40)
    plot_tpcf(bins,corr,dcorr,axs=axs,fit=None,
              label=r'$R_{\mathrm{max}}=3r_c$',fmt='.-',lw=0.2,c='#F29A0C',
             errorevery=5)
    popt,perr = fit_tpcf(bins,corr,dcorr,fit='singletrunc')
    axs.axvline(popt[2],ls='--',c='#F29A0C',lw=2.0)

    x,y = Create_Galaxy(r_c=r_c,R_max=2.0,i=i,r_min=r_min)
    bins,corr,dcorr = compute_tpcf(x,y,R_max=1.0,R_min=1.e-2)
    plot_tpcf(bins,corr,dcorr,axs=axs,fit=None,
              label=r'$R_{\mathrm{max}}=5r_c$',fmt='.-',lw=0.2,c='#6A00DB',
             errorevery=5)
    popt,perr = fit_tpcf(bins,corr,dcorr,fit='singletrunc')
    axs.axvline(popt[2],ls='--',c='#6A00DB',lw=2.0)

    axs.arrow(0.2,0.05,dx=0.0,dy=-0.038,color='#5C5850',
                  width=0.005,length_includes_head=True,
                  head_width=0.04,head_length=0.007,lw=2.0,alpha=0.5)
    axs.text(0.18,0.07,r'$r_c$',fontsize=18,c='#5C5850')

    axs.legend(loc='lower left')
    axs.set_ylabel(r"$\omega$")

    if(save):
        filename = outdir+'FigC4.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    #Parsing Arguments
    ############################################################################
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for tpcf-starclusters. All inputs optional. ')

    ap.add_argument('-outdir',action='store',type=str,default='../results/',
        help = 'Alternate output directory for plots.')

    ap.add_argument('-indir',action='store',type=str,default='../results/Menonetal2021/',
        help = 'Alternate input directory where summary pickle files stored.')
    args = vars(ap.parse_args())
    

    print("Preparing plots of Menon et al 2021b.")
    print("Using summary pickle files from the directory: {}".format(args['indir']))
    print("Resulting plots would be saved in {}".format(args['outdir']))
    Figure1(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure2(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure3(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure4(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure5(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure6(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure7(save=True,outdir=args['outdir'],indir=args['indir'])
    Figure8(save=True,outdir=args['outdir'],indir=args['indir'])

    print("Appendix plots not added yet. Work in progress...")
