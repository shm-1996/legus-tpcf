from header import *
from astropy.coordinates import SkyCoord
from TPCF import bootstrap_two_point,bootstrap_two_point_angular
import copy

# Some defaults : 20 bins in rang 10-5000 pc
default_bin_limits = 10,5000
default_no_of_bins = 20

class Galaxy(object):
    """
    A class describing the properties of the galaxy, and
    providing methods to perform calculations using those properties.

    Parameters
        galaxy_name : string
          Name of galaxy 
        filename : string
          name of file from which to read galaxy description
        verbose : Boolean
              print out information about the galaxy as we read it
    Class attributes - ordered by category
    ########################################################
    I. Basic
        name : string
            Name of the galaxy
        distance : float
            Distance to the galaxy in Mpc 
        errordist : float
            Error in distance to the galaxy in Mpc
        inclination : float
            Inclination of the galaxy in the line-of-sight
        pa : float
            Position angle of the galaxy
        inclination_flag : bool
            Flag to perform correction for inclination in the galaxies
        inclination_done : bool
            Flag to inform that inclination correction already done
        centre : tuple
            ra,dec of galaxy centre in fk5 notation
        catalog_file : string
            filename of the cluster catalog file
        fits_file : string
            filename of HST image
        region_file : string 
            filename of the ds9 region file 
        morph_type: string 
            Morphological type of galaxy 
        T_value : float
            Morphological T-type of galaxy
        sfr : float
            Star formation rate of galaxy
        mstar : float
            Stellar mass in the galaxy
        mhi : float
            HI mass in the galaxy
        r25: float
            R25 radius
        sigma_sfr: float
            SFR surface density
        sigma_hi: float
            HI mass surface density
        ra : ndarray
            RA of star clusters
        dec: ndarray
            DEC of star clusters
        ra_raw: ndarray
            RA of star clusters before inclination correction
        dec_raw: ndarray
            DEC of star clusters before inclination correction
        ########################################################

        II. TPCF related
        no_bins : integer
            Number of bins to compute TPCF
        bin_limits : tuple
            limits of the TPCF bins in arcsecs
        bin_centres : ndarray
            Centres of the bins of TPCF
        corr : ndarray
            TPCF values for all clusters
        dcorr: ndarray 
            TPCF value errors for all clusters
        ycorr : ndarray
            TPCF values for young clusters
        ydcorr: ndarray 
            TPCF value errors for young clusters
        ocorr : ndarray
            TPCF values for old clusters
        odcorr: ndarray 
            TPCF value errors for old clusters
        ########################################################
        III. MCMC fit related - 
        a. replace y with o for young clusters
        b. singlepl with piecewise oor singletrunc for Model PW or PF

        fit_values_singlepl_y : tuple
            Best-fit values for Model S in young cluster TPCF
        fit_errors_singlepl_y : tuple
            Errors on best-fit values for Model S in young cluster TPCF
        median_fit_singlepl_y : tuple
            Median fit values for Model S
        low_error_singlepl_y : tuple
            16th percentile error for Model S
        high_error_singlepl_y : tuple
            84th percentile error for Model S
        beta_distribution: ndarray
            Full posterior distribution of beta parameter in Model PW
        theta_distribution : ndarray
            Full posterior distribution of theta parameter in Model PF
        aic_singlepl_y : float
            AIC value for Model S in young cluster TPCF
        best_model_y: string
            Best fit model obtained with AIC criterion for young clusters
        best_fit_y : tuple
            The best-fit parameters corresponding to the best fit model
        best_fit_error_y: tuple
            Error on best-fit parameters corresponding to best model
        #############################################################
        IV. Derived physical parameters
        lcorr: tuple
            Correlation length and its 16th and 84th percentile errors
        alpha: tuple
            Slope of TPCF and its 16th and 84th percentile errors
        D2: tuple 
            Fractal dimension of young star clusters and its 16th and 84th percentile errors
        rc_disk: tuple
            Derived exponential disk scale radius and its 16th and 84th percentile errors
    """


    ####################################################################
    # Method to initialize
    ####################################################################
    def __init__(self, galaxy_name=None, filename=None,verbose=False):
        """
        Parameters
            galaxy_name : string
              Name of galaxy
            filename : string
              name of file from which to read galaxy properties
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None

        """
        self.name= None
        self.distance = 0.0
        self.inclination = 0.0
        self.deproject_galaxy = True
        self.inclination_done = False
        self.centre = 0,0
        self.no_bins = 0
        self.bin_limits = [0,0]
        self.catalog_file = None
        self.outdir = None
        self.fits_file = None

        # Read from info file if provided
        if(filename != None):
            self.readFile(filename,verbose=verbose)

        #If not provided see if name of galaxy provided
        if(filename is None):
            if(galaxy_name):
                self.name = galaxy_name.strip()
                self.defaultGalaxyInfo(verbose=verbose)
            else:
                raise myError("No provided galaxy or information file.")

        self.get_data()
        # Obtain RA/DEC of all clusters
        self.get_ra_dec()
        # Set bins based on bin limits & no of bins
        self.set_bins()
        # Read galaxy properties
        self.read_galaxyprops()


    ####################################################################
    # Method to set default galaxy properties
    ####################################################################

    def defaultGalaxyInfo(self,verbose=False):
        """
        Sets default galaxy information by reading its info file, or creating it if
        it does not exist. This info is required for computing the TPCF.
        Parameters
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None
        """
        
        if self.name in list_of_galaxies:
            if(verbose):
                print("Starting analysis for "+self.name)

            filename = os.path.abspath('../data/galaxyinfo/'+
                self.name+'.info')
            #Check catalog file exists else throw error
            if(not os.path.exists(filename)) :
                #Create file if it does not exist
                if(verbose):
                    print("Creating galaxy info file as it does not exist.")
                self.create_file(verbose=verbose)
            else:
                self.readFile(filename=filename,verbose=verbose)

            self.outdir = os.path.abspath('../results/')
            if(verbose):
                print("Setting output directory to {}".format(self.outdir))                        
        else :
            raise myError("The galaxy information+catalog is not available."+
                " Create information file or manually input info.")



    ####################################################################
    # Method to create galaxy info file
    ####################################################################

    def create_file(self,verbose=False):
        """
        Creates info file for galaxy, by reading information from Table I
        of Calzetti et al 2015 (LEGUS). 
        
        Parameters
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None

        """
        info_directory = os.path.abspath('../data/galaxyinfo')
        
        #Read galaxy info from Table 1 Calzetti et al 2015
        if(verbose):
            print("Reading galaxy information from Table 1 of Calzetti et al 15.")
        Legus_Table = info_directory+'/Calzetti_Table.txt'
        NewDistances_Table = info_directory+'/GalaxyProps_New.dat'
        list_gals = np.array(list_of_galaxies)
        index_listgals = np.where(self.name == list_gals)

        # Read columns from table
        galaxy_names = np.loadtxt(Legus_Table,usecols=0,delimiter='\t',dtype=str)
        #Convert self galaxy name to table name format
        galaxy_formatted = self.name.split("_")[1]
        galaxy_formatted = '%04d'%int(galaxy_formatted)
        galaxy_formatted = self.name.split("_")[0] + ' ' + galaxy_formatted
        index = np.where(galaxy_formatted == galaxy_names)
        
        
        # Read in distance, inclination and position angles from new table
        galaxy_distances = np.loadtxt(NewDistances_Table,usecols=0)
        galaxy_errordist = np.loadtxt(NewDistances_Table,usecols=1)
        galaxy_inclinations = np.loadtxt(NewDistances_Table,usecols=2)
        galaxy_pa = np.loadtxt(NewDistances_Table,usecols=3)
        distance = galaxy_distances[index_listgals][0]
        errordist = galaxy_errordist[index_listgals][0]
        inclination = galaxy_inclinations[index_listgals][0]
        pa = galaxy_pa[index_listgals][0]

        #Setting default bin info, defined in header.py

        #Treat some exceptions in bin limits
        if(self.name == 'NGC_3344'):
            bin_limits = 10,2000 
            no_bins = 15
        elif(self.name == 'NGC_3738'):
            bin_limits = 10,2000
            no_bins = 20
        elif(self.name == 'NGC_5194'):
            bin_limits = 20,5000
            no_bins = 20
        elif(self.name == 'NGC_5253'):
            bin_limits = 50,1000
            no_bins = 20
        else:
            bin_limits = default_bin_limits
            no_bins = default_no_of_bins
        if(verbose):
            print("Setting default bin limits of {}".format(bin_limits))            
            print("Setting default no of bins = {}".format(no_bins))
        catalog_file = os.path.abspath('../data/' + 
            self.name+'.tab')
        if(verbose):
            print("Setting Catalog file = {}".format(catalog_file))
        fits_file = os.path.abspath('../data/' + 
            self.name+'.fits')
        if(verbose):
            print("Setting fits file = {}".format(fits_file))
        region_file = os.path.abspath('../data/') + '/' +\
                self.name + '.reg'
        if(verbose):
            print("Setting region file = {}".format(region_file))  
        #Write to info file
        info_file = info_directory + '/' + self.name + '.info'
        #Write in required format
        try: 
            fp = open(info_file,'w')
        except IOError:
            raise myError("Cannot write to file "+info_file)
        

        self.distance = distance 
        self.errordist = errordist
        self.inclination = inclination
        self.no_bins = no_bins
        self.bin_limits = bin_limits 
        self.catalog_file = catalog_file
        self.fits_file = fits_file
        self.pa = pa
        self.region_file = region_file

        if(verbose):
            print("Writing info file for {}".format(self.name))

        
        # Write header
        fp.write("## Info file for {} with relavant info \
for computing TPCF.\n".format(self.name))
        fp.write('################################################# \n\n')

        # Name of galaxy
        fp.write("# Name of Galaxy \n")
        fp.write("Name = {} \n\n".format(self.name))

        # Distance
        fp.write("# Distance in Mpc \n")
        fp.write("Distance = {} \n\n".format(self.distance))

         # Error in Distance
        fp.write("# Error in Distance in Mpc \n")
        fp.write("Error_Distance = {} \n\n".format(self.errordist))

        # Inclination
        fp.write("# Inclination in degrees \n")
        fp.write("Inclination = {} \n\n".format(self.inclination))

        # Position Angle
        fp.write("# Position Angle in degrees \n")
        fp.write("Position_Angle = {} \n\n".format(self.pa))

        # Number of bins
        fp.write("# Number of bins \n")
        fp.write("No_Bins = {} \n\n".format(self.no_bins))

        #Bin Limits
        fp.write("# Bin Limits in parsec\n")
        fp.write("Bin_Low = {} \n".format(self.bin_limits[0]))
        fp.write("Bin_High = {} \n\n".format(self.bin_limits[1]))

        # Catalog file
        fp.write("# Catalog file \n")
        fp.write("Catalog_File = {} \n\n".format(os.path.abspath(self.catalog_file)))

        # Fits file
        fp.write("# Fits file \n")
        fp.write("Fits_File = {} \n\n".format(os.path.abspath(self.fits_file)))

        # Region file
        fp.write("# Region file \n")
        fp.write("Region_File = {} \n\n".format(os.path.abspath(self.region_file)))

        #Close file
        fp.close()


    ####################################################################
    # Method to read galaxy metadata require to compute TPCF from a file
    ####################################################################
    def readFile(self,filename,verbose=False): 
        #Try reading file
        try :
            fp = open(filename,'r')
        except IOError :
            raise myError("cannot open file "+filename)
        for line in fp : 
            # Skip empty and comment lines
            if line=='\n':
                continue
            if line.strip()[0] == "#":
                continue
             # Break line up based on equal sign
            linesplit = line.split("=")
            if len(linesplit) < 2:
                raise myError("Error parsing input line: "+line)
            if linesplit[1] == '':
                raise myError("Error parsing input line: "+line)

            # Trim trailing comments from portion after equal sign
            linesplit2 = linesplit[1].split('#')

            # Read stuff based on the token that precedes the equal sign
            if linesplit[0].upper().strip() == 'NAME':
                self.name = str(linesplit2[0]).strip()
                if(verbose) :
                    print("Setting galaxy = "+str(self.name))
            elif linesplit[0].upper().strip() == 'DISTANCE':
                self.distance = float(linesplit2[0])
                if(verbose):
                    print("Setting distance = "+str(self.distance) + " Mpc")
            elif linesplit[0].upper().strip() == 'ERROR_DISTANCE':
                self.errordist = float(linesplit2[0])
                if(verbose):
                    print("Setting error in distance = "+str(self.errordist) + " Mpc")
            elif linesplit[0].upper().strip() == 'INCLINATION':
                self.inclination = float(linesplit2[0])
                if(verbose):
                    print("Setting inclination = "+str(self.inclination) + " degrees")
            elif linesplit[0].upper().strip() == 'POSITION_ANGLE':
                self.pa = float(linesplit2[0])
                if(verbose):
                    print("Setting position angle = "+str(self.pa) + " degrees")
            elif linesplit[0].upper().strip() == 'NO_BINS':
                self.no_bins = int(linesplit2[0])
                if(verbose):
                    print("Setting number of bins = "+str(self.no_bins))
            elif linesplit[0].upper().strip() == 'BIN_LOW':
                self.bin_limits[0] = float(linesplit2[0])
                if(verbose):
                    print("Setting lower bin limit for TPCF = "+str(self.bin_limits[0]) +
                         " parsec")
            elif linesplit[0].upper().strip() == 'BIN_HIGH':
                self.bin_limits[1] = float(linesplit2[0])
                if(verbose):
                    print("Setting upper bin limit for TPCF = "+str(self.bin_limits[1]) +
                        " parsec")
            elif linesplit[0].upper().strip() == 'CATALOG_FILE':
                self.catalog_file = os.path.abspath(str(linesplit2[0]).strip())
                if(verbose):
                    print("Setting cluster catalog file = "+str(self.catalog_file))

            elif linesplit[0].upper().strip() == 'FITS_FILE':
                self.fits_file = os.path.abspath(str(linesplit2[0]).strip())

                if(verbose):
                    print("Setting HST fits file = "+str(self.fits_file))

            elif linesplit[0].upper().strip() == 'REGION_FILE':
                self.region_file = os.path.abspath(str(linesplit2[0]).strip())

                if(verbose):
                    print("Setting HST region file = "+str(self.region_file))

            else:
                # Line does not correspond to any known keyword, so
                # throw an error
                raise myError("unrecognized token " +
                    linesplit[0].strip() + " in file " + filename)

        # Close file
        fp.close()

        #Compute ra/dec of centre of galaxy
        name = self.name.split('_')[0] + ' ' + self.name.split('_')[1]
        ra_dec = SkyCoord.from_name(name)
        ra = ra_dec.ra.value
        dec = ra_dec.dec.value
        self.centre = ra,dec

        # Some safety checks below
        #Make sure lower and upper limits of bins are consistent
        if(self.bin_limits[0] >= self.bin_limits[1]) :
            raise myError("Provided bin input for lower limit greater than higher limit.")


    def read_galaxyprops(self):
        """
        Reads galaxy properties from other ancillary data/references. This is not
        required directly for computing the TPCF.

        """
        info_directory = os.path.abspath('../data/galaxyinfo')
        Legus_Table = info_directory+'/Calzetti_Table.txt'

        #Convert self galaxy name to table name format
        galaxy_formatted = self.name.split("_")[1]
        galaxy_formatted = '%04d'%int(galaxy_formatted)
        galaxy_formatted = self.name.split("_")[0] + ' ' + galaxy_formatted

        galaxy_names = np.loadtxt(Legus_Table,usecols=0,delimiter='\t',dtype=str)
        index = np.where(galaxy_formatted == galaxy_names)


        morph_type = np.loadtxt(Legus_Table,usecols=2,delimiter='\t',dtype=str)
        self.morph_type = morph_type[index][0]

        T_value = np.loadtxt(Legus_Table,usecols=3,delimiter='\t',dtype=str)
        self.T_value = float(T_value[index][0].split('(')[0])

        #Account for different formats in table for different type of galaxies
        if(self.name in ['NGC_1313',]):
            sfr = np.loadtxt(Legus_Table,usecols=10,delimiter='\t',dtype=str)
            mstar = np.loadtxt(Legus_Table,usecols=11,delimiter='\t',dtype=str)
            mhi = np.loadtxt(Legus_Table,usecols=12,delimiter='\t',dtype=str)
        elif(self.name in ['NGC_3738','NGC_4449','NGC_4656','NGC_6503','NGC_5457','NGC_5253']):
            sfr = np.loadtxt(Legus_Table,usecols=10,delimiter='\t',dtype=str)
            mstar = np.loadtxt(Legus_Table,usecols=11,delimiter='\t',dtype=str)
            mhi = np.loadtxt(Legus_Table,usecols=12,delimiter='\t',dtype=str)
        else :
            sfr = np.loadtxt(Legus_Table,usecols=11,delimiter='\t',dtype=str)
            mstar = np.loadtxt(Legus_Table,usecols=12,delimiter='\t',dtype=str)
            mhi = np.loadtxt(Legus_Table,usecols=13,delimiter='\t',dtype=str)

        self.sfr = float(sfr[index][0])
        self.mstar = float(mstar[index][0])
        self.mhi = float(mhi[index][0])

        R25File = info_directory+'/R25.txt'
        galaxy_names = np.loadtxt(R25File,usecols=0,delimiter='\t',dtype=str)
        index = np.where(galaxy_formatted == galaxy_names)
        galaxy_r25 = np.loadtxt(R25File,usecols=1,delimiter='\t',dtype=float)[index][0]

        #Convert from arcmin to degree
        galaxy_r25 = galaxy_r25*60*u.arcsec.to(u.radian)
        galaxy_r25 = galaxy_r25 * self.distance*const.Parsec*1.e6
        #Convert from degrees to arcsec

        # Store in kpc
        self.r25 = galaxy_r25/(const.Parsec*1.e3)

        # Store in Msun per year per kpc^2
        self.sigma_sfr = self.sfr/(np.pi * self.r25**2)

        #Store sigma_hi in Msun per pc^2
        self.sigma_hi = self.mhi/(np.pi * (self.r25*1.e3)**2)

        #TODO: Add sigma_h2 and vrot here


    ####################################################################
    # Method to obtain ra dec of clusters
    ####################################################################
    def get_ra_dec(self,cluster_class=-1,age=None,age_cut=1.e7,verbose=False):
        """
        Method to read the RA/DEC coordinates for the clusters in the catalog
        Parameters:
            cluster_class : integer
                The cluster class to read. Use -1 to read Class 1+2+3. 
            age: string 
                Age group to compute the TPCF on : 'young', or 'old'
            age_cut: float
                The cutoff age below/above which young/old is defined
            verbose: boolean
                print out what is being done
        Returns:
            None

        """
        file = np.loadtxt(self.catalog_file)
        Class0_sources = np.where(file[:,33]==0)

        # Pick out Machine Learning Catalogs for NGC 5194 for Class 1 & 2
        if(self.name == 'NGC_5194'):
            Class1_sources = np.where((file[:,33]==1) | (file[:,34]==1))
            Class2_sources = np.where((file[:,33]==2) | (file[:,34]==2))
        else :
            Class1_sources = np.where(file[:,33]==1)
            Class2_sources = np.where(file[:,33]==2)
        
        Class3_sources = np.where(file[:,33]==3)
        Class4_sources = np.where(file[:,33]==4)
        Cluster_sources = np.append(Class1_sources,Class2_sources)
        if(self.name == 'NGC_5194'):
            print("Class 3 Sources ignored for NGC 5194.")
        else:
            Cluster_sources = np.append(Class3_sources,Cluster_sources)


        
        # Compute TPCF for a subset of clusters if passed in argument
        if(cluster_class == 1) :
            RA = file[Class1_sources][:,3]
            DEC = file[Class1_sources][:,4]
        elif(cluster_class == 2) :
            RA = file[Class2_sources][:,3]
            DEC = file[Class2_sources][:,4]
        elif(cluster_class == 3) :
            RA = file[Class3_sources][:,3]
            DEC = file[Class3_sources][:,4]
        elif(cluster_class == -1) :
            RA = file[Cluster_sources][:,3]
            DEC = file[Cluster_sources][:,4]

        else :
            raise myError("Invalid cluster class passed.")

        if(age is not None):
            ages = self.get_cluster_ages(cluster_class=cluster_class)
            if(age == 'young'):
                RA = RA[np.where(ages <= age_cut)]
                DEC = DEC[np.where(ages <= age_cut)] 
            elif(age == 'old'):
                RA = RA[np.where(ages > age_cut)]
                DEC = DEC[np.where(ages > age_cut)]
            
        self.ra = RA 
        self.dec = DEC

        self.ra_raw = RA 
        self.dec_raw = DEC

        # Deproject positions if deproject_galaxy is set to true. 
        # Note this is true by default. 
        if(self.deproject_galaxy == True):
            self.correct_inclination(force=True,verbose=verbose)

    ####################################################################
    # Method to obtain ages of clusters
    ####################################################################

    def get_cluster_ages(self,cluster_class=-1,verbose=False):
        """
        Get the ages for the clusters
        Parameters:
            cluster_class : integer
                The cluster class to read. Use -1 to read Class 1+2+3.           
            verbose: boolean
                print out what is being done
        Returns:
            ages: ndarray
                Ages of the clusters as an ndarray

        """

        file = np.loadtxt(self.catalog_file)
        Class0_sources = np.where(file[:,33]==0)

        # Pick out Machine Learning Catalogs for NGC 5194 for Class 1 & 2
        if(self.name == 'NGC_5194'):
            Class1_sources = np.where((file[:,33]==1) | (file[:,34]==1))
            Class2_sources = np.where((file[:,33]==2) | (file[:,34]==2))
        else :
            Class1_sources = np.where(file[:,33]==1)
            Class2_sources = np.where(file[:,33]==2)
        
        Class3_sources = np.where(file[:,33]==3)
        Class4_sources = np.where(file[:,33]==4)
        Cluster_sources = np.append(Class1_sources,Class2_sources)
        if(self.name == "NGC_5194"):
            print("Class 3 sources ignored in ages.")
        else:
            Cluster_sources = np.append(Class3_sources,Cluster_sources)
        
        # Compute TPCF for a subset of clusters if required
        if(cluster_class == 1) :
            ages = file[Class1_sources][:,16]            
        elif(cluster_class == 2) :
            ages = file[Class2_sources][:,16]
        elif(cluster_class == 3) :
            ages = file[Class3_sources][:,16]
        elif(cluster_class == -1) :
            ages = file[Cluster_sources][:,16]
        else :
            raise myError("Invalid cluster class passed.")

        return ages


    ####################################################################
    # Method to obtain masses of clusters
    ####################################################################

    def get_cluster_masses(self,cluster_class=-1,verbose=False):
        """
        Get the masses for the clusters
        Parameters:
            cluster_class : integer
                The cluster class to read. Use -1 to read Class 1+2+3.           
            verbose: boolean
                print out what is being done
        Returns:
            masses: ndarray
                Masses of the clusters as an ndarray

        """
        file = np.loadtxt(self.catalog_file)
        Class0_sources = np.where(file[:,33]==0)

        # Pick out Machine Learning Catalogs for NGC 5194 for Class 1 & 2
        if(self.name == 'NGC_5194'):
            Class1_sources = np.where((file[:,33]==1) | (file[:,34]==1))
            Class2_sources = np.where((file[:,33]==2) | (file[:,34]==2))
        else :
            Class1_sources = np.where(file[:,33]==1)
            Class2_sources = np.where(file[:,33]==2)
        
        Class3_sources = np.where(file[:,33]==3)
        Class4_sources = np.where(file[:,33]==4)
        Cluster_sources = np.append(Class1_sources,Class2_sources)
        if(self.name == "NGC_5194"):
            print("Class 3 sources ignored in ages.")
        else:
            Cluster_sources = np.append(Class3_sources,Cluster_sources)
        
        # Compute TPCF for a subset of clusters if required
        if(cluster_class == 1) :
            masses = file[Class1_sources][:,19]            
        elif(cluster_class == 2) :
            masses = file[Class2_sources][:,19]
        elif(cluster_class == 3) :
            masses = file[Class3_sources][:,19]
        elif(cluster_class == -1) :
            masses = file[Cluster_sources][:,19]
        else :
            raise myError("Invalid cluster class passed.")

        return masses
        

    ####################################################################
    # Method to obtain set bins.
    ####################################################################
    def set_bins(self,set_no_bins=None,set_bin_limits=None):
        """
        Method to set the number of bins, allowing user to change no of bins
        and bin limits, in which case bins are recomputed. 
        Parameters
            set_no_bins : integer
              Number of bins to set
            set_bin_limits : tuple
              Bin limits to set provided as a tuple in units of parsecs
        Returns
            None

        """
        
        bin_limits = [0,0]
        if(set_no_bins == None and set_bin_limits == None):
            #Convert bin limits in parsec to arcsec
            distance = self.distance*const.Parsec*1.e6
            bin_limits[0]= self.bin_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
            bin_limits[1] = self.bin_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

            bin_min,bin_max = np.log10(bin_limits[0]*u.arcsec.to(u.deg)),\
            np.log10(bin_limits[1]*u.arcsec.to(u.deg))
            bins = np.logspace(bin_min,bin_max,self.no_bins+1)
            self.bins = bins
            self.bin_centres = (self.bins[1:]+self.bins[:-1])/2
            self.bins_arcsec = bins*(1./arcsec_to_degree)

        else :
            if(set_no_bins is not None):
                no_bins = int(set_no_bins)
                self.no_bins = no_bins
                print("Changing number of bins to {}".format(set_no_bins))
            if(set_bin_limits is not None):
                set_bin_limits = list(set_bin_limits)
                self.bin_limits = list(self.bin_limits)
                self.bin_limits[0] = float(set_bin_limits[0])
                self.bin_limits[1] = float(set_bin_limits[1])
                print("Changing bin limits to {}".format(set_bin_limits))

            print("Recomputing bins")
            #Convert bin limits in parsec to arcsec
            distance = self.distance*const.Parsec*1.e6
            bin_limits[0]= self.bin_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
            bin_limits[1] = self.bin_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

            bin_min,bin_max = np.log10(bin_limits[0]*u.arcsec.to(u.deg)),\
            np.log10(bin_limits[1]*u.arcsec.to(u.deg))
            bins = np.logspace(bin_min,bin_max,self.no_bins+1)
            self.bins = bins
            self.bin_centres = (self.bins[1:]+self.bins[:-1])/2
            self.bins_arcsec = bins*(1./arcsec_to_degree)



    ####################################################################
    # Method to correct for inclination.
    ####################################################################
    def correct_inclination(self,force=False,verbose=True):
        """
        Method to correct for the inclination. See Equation 5 & 6 of Menon et al 2021b
        Parameters:
            force: boolean
                Force perform the inclination irrespective of whether correction done
                or not. Use this only if the raw RA/DEC positions of clusters are changed,
                as otherwise inclination would be repeated, and the final positions wrong.
            verbose: boolean
                print out what is being done
        Returns:
            None

        """
        if(force == False):
            if(self.inclination_done == True):
                raise myError("Inclination correction already performed." +
                    " Use flag force =True if you want to force deproject.")
        if(verbose):
            print("Correcting for inclination of galaxy. This galaxy has PA = {}"
                .format(self.pa) + " and inclination = {}".format(self.inclination))


        hdu = fits.open(self.fits_file)[0]
        wcs = WCS(hdu.header)

        #Convert RA/DEC to appropriate pixels in the WCS 
        xpix,ypix = wcs.all_world2pix(self.ra_raw,self.dec_raw,0)

        #Coordinates of centre of galaxy
        name_galaxy = self.name.split('_')[0] + ' ' + self.name.split('_')[1]
        ra_dec = SkyCoord.from_name(name_galaxy)
        xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)

        # See Eq 1 Grasha et al 2017
        #Basically rotating all points by PA in clockwise direction to align major axis with north. 
        ra_dep = (xpix-xpix_c)*np.cos(np.deg2rad(self.pa)) + (ypix-ypix_c)*np.sin(np.deg2rad(self.pa))
        dec_dep = -1.*(xpix-xpix_c)*np.sin(np.deg2rad(self.pa)) + (ypix-ypix_c)*np.cos(np.deg2rad(self.pa))

        #Correct for inclination : separation only along minor axis changes
        ra_dep = ra_dep/(np.cos(np.deg2rad(self.inclination)))

        ra_dep = ra_dep+xpix_c
        dec_dep = dec_dep+ypix_c

        #Convert from pixel back to RA/DEC
        ra_dep,dec_dep = wcs.all_pix2world(ra_dep,dec_dep,0)

        self.ra = ra_dep 
        self.dec = dec_dep

        self.inclination_done = True

    ####################################################################
    # Method to create ds9 region file for FOV footprint
    ####################################################################

    def create_region_file(self,verbose=False):
        """
        Creates ds9 region file for galaxy, used to prepare the random catalog.
        This is a wrapper to call the footprintfinder.py tool developed by StSCI. 
        Parameters
            verbose: Boolean
                Print out information

        Returns
            None

        """
        from footprintfinder import main
        fits_file = self.fits_file
        if(verbose):
            print("Creating region file....")
        #Call footprint finder
        main('-d ' + os.path.abspath(self.fits_file))

        #Move region file to appropriate directory
        region_filebase = self.fits_file.split(os.path.dirname(self.fits_file)+'/')[1]
        region_filebase = region_filebase.split('.fits')[0]
        region_filename = region_filebase + '_footprint_ds9_image.reg'

        shutil.copy(os.getcwd() + '/' + region_filename,self.region_file)

        if(verbose):
            print("Saved region file in {}".format(self.region_file))

        #Remove other output files 
        if(verbose):
            print("Deleting files created by footprint in local directory.")
        unwanted_file = region_filebase + '_footprint.txt'
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)
        unwanted_file = region_filebase + '_footprint_ds9_linear.reg'
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)
        unwanted_file = region_filename
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)


     ####################################################################
    # Method to compute TPCF for the galaxy
    ####################################################################

    def Compute_TPCF(self,cluster_class=-1,save=False,
        random_method='masked',verbose=False,read_positions=True,read_bins=True,
        tpcf_method='landy-szalay',Nbootstraps=100,age=None,age_cut=1.e7):
        """
        Computes the Two-Point Correlation Function of the star clusters in the 
        galaxy. 
        Parameters
            filename : string
              name of file from which to read galaxy description
            cluster_class : integer
              class of clusters for which TPCF to be computed
              can be 1,2,3 or -1 for 1+2+3 combined
            save : boolean
              flag to save bootstrap realisations of the TPCF
            random_method : string
                method to use to prepare the random catalog
            verbose : boolean
                print out what is being done
            read_positions: boolean
                Reads default position of clusters. Switch this off, 
                if user provided positions for the cluster are provided, 
                as it would be overwrited if true
            read_bins: boolean
                Reads default bins for the TPCF. Switch this off, 
                if user provided bins, as it would be overwrited if true
            tpcf_method : string
                method of computing TPCF. Options are standard or landy-szalay
            Nbootstraps: integer 
                Number of bootstrap samples to compute the TPCF on
            age: string 
                Age group to compute the TPCF on : 'young', or 'old'
            age_cut: float
                The cutoff age below/above which young/old is defined

        Returns
            None

        """
        if(self.catalog_file is None) :
            raise myError("No catalog file defined for galaxy.")

        if(verbose):
            if(cluster_class==-1):
                print("Computing TPCF for {} for cluster class 1+2+3 using {} random method."
                .format(self.name,random_method))
            else:
                print("Computing TPCF for {} for cluster class {} using {} random method."
                    .format(self.name,cluster_class,random_method))
        
        if(verbose and read_positions is True):
            print("Reading cluster positions from cluster catalog in {}"
                .format(self.catalog_file))
        if(read_positions is True):
            self.get_ra_dec(cluster_class=cluster_class,
                age=age,age_cut=age_cut)
        if(read_bins is True):
            self.set_bins()

        # Safety check for masked_radial method
        if random_method in ['masked_radial','masked','exponential']:
            self.region_file = os.path.abspath('../data/') + '/' +\
                self.name + '.reg'

            #If region file not present, create it
            if(not os.path.exists(self.region_file)):
                self.create_region_file(verbose=verbose)
            else :
                if(verbose):
                    print("Region file exists. Using region file {}".format(
                        self.region_file))

        if(random_method == 'exponential'):
            r_c = self.rc
        else:
            r_c = None
        corr,dcorr,bootstraps = bootstrap_two_point_angular(self,
                            method=tpcf_method,Nbootstraps=Nbootstraps,
                            random_method=random_method,r_c=r_c)
        if(verbose):
            print("TPCF computation completed.")
        

        if(age == 'young'):
            self.ycorr = corr
            self.ydcorr = dcorr 
            self.ybootstraps = bootstraps
        elif(age == 'old'):
            self.ocorr = corr
            self.odcorr = dcorr 
            self.obootstraps = bootstraps
        else:
            self.corr = corr 
            self.dcorr = dcorr
            self.bootstraps = bootstraps

    ####################################################################
    # Method to fit power law to TPCF
    #TODO: MCMC fit. 
    ####################################################################

    def fit_TPCF(self,method='mcmc',function='singlepl',age=None,
        plot=False,omega1=True,use_bounds=True,N=1000):
        """
        Fits a functional form to the computed TPCF.
        Parameters
            method: string
                Method used to fit a function to the TPCF. Options are "single"
                for a single instance curve fit, "bootstrap" for performing 
                bootstrap resampling fits, "MCMC" to use an MCMC fitting procedure
            function : string
                function to fit to the TPCF. Can be 'singlepl','piecewise' and 'singletrunc'
            age: string
                Can be 'young', 'old' or None. Chooses the TPCF to fit to.
            plot: Boolean
                flag to plot diagnostic plots for MCMC
            omega1: Boolean
                flag to use 1+Omega as the function to or just Omega
            use_bounds: Boolean
                flag to set bounds on parameters. Not used in MCMC. 
            N : integer
                Number of bootstrap realisations for fitting if using the method. 

        Returns
            None

        """
        from TPCF import linear_function,linear_truncation,\
        smooth_function,onepowerlaw_function,piecewise_truncation
        from MCMCfit import fit_MCMC
        if(method not in ['single','bootstrap','mcmc']) :
            raise ValueError("Method for fitting power law should be one of the"+
            " following: single, bootstrap or mcmc")
        if(function not in ['piecewise','smooth','singletrunc','singlepl']):
            raise ValueError("This funtional form does not exist.")
        
        bins = self.bins_arcsec
        separation_bins = (bins[1:]+bins[:-1])/2
        
        separation_bins = separation_bins.astype(np.float)
        indices = np.where(self.corr>0.0)
        corr_fit = self.corr[indices].astype(np.float)
        dcorr_fit = self.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)


        #Parameter limits
        beta_limits = [50.0, 300.0]

        distance = self.distance*const.Parsec*1.e6
        #beta limits in arcseca
        beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
        beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

        if(function=='piecewise'):
            bounds = ([-10.0,-3.0,-3.0,np.log(beta_limits[0])],[10.0,0.0,0.0,np.log(beta_limits[1])])
        elif(function == 'smooth'):
            bounds = ([-10.0,-3.0,-3.0,beta_limits[0]],[10.0,0.0,0.0,beta_limits[1]])
        elif(function == 'singletrunc'):
            bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])

        if(method == 'single') :
            if(function == 'piecewise'):
               
                #bounds = ([-10,10],[-3.0,0.0],[-3.0,0.0],[np.log(beta_limits[0]),np.log(beta_limits[1])])
                if(omega1):
                    popt,pcov = curve_fit(linear_function,separation_bins,
                        np.log(1+corr_fit),sigma=dcorr_fit/(1+corr_fit),bounds=bounds)
                else:
                    popt,pcov = curve_fit(linear_function,separation_bins,
                        np.log(corr_fit),sigma=dcorr_fit/corr_fit,bounds=bounds)
            elif(function == 'smooth'):
                popt,pcov = curve_fit(smooth_function,separation_bins,
                    corr_fit,sigma=dcorr_fit,bounds=bounds)

            elif(function == 'singletrunc'):
                if(omega1):
                    popt,pcov = curve_fit(linear_truncation,separation_bins,
                                np.log(1+corr_fit),sigma=dcorr_fit/(1+corr_fit),bounds=bounds)
                else:
                    popt,pcov = curve_fit(linear_truncation,separation_bins,
                                np.log(corr_fit),sigma=dcorr_fit/corr_fit,bounds=bounds)
            elif(function == 'singlepl'):
                if(omega1):
                    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                                np.log(1+corr_fit),sigma=dcorr_fit/(1+corr_fit))
                else:
                    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                                np.log(corr_fit),sigma=dcorr_fit/corr_fit)

                  
            self.fit_values = popt
            self.fit_errors = np.sqrt(np.diag(pcov))

        elif(method == 'bootstrap') :
            fit_bootstraps = np.zeros((N,4))
            for i in range(N):
                y_fit = corr_fit + dcorr_fit*np.random.randn(1)

                try:
                #Fit to this
                    if(function == 'piecewise'):
                        if(use_bounds == True):
                            if(omega1):
                                popt,pcov = curve_fit(linear_function,separation_bins,
                                    np.log(1+y_fit),sigma=dcorr_fit/corr_fit,bounds=bounds)
                            else:
                                popt,pcov = curve_fit(linear_function,separation_bins,
                                    np.log(y_fit),sigma=dcorr_fit/corr_fit,bounds=bounds)
                        else:
                            popt,pcov = curve_fit(linear_function,separation_bins,
                                np.log(y_fit),sigma=dcorr_fit/corr_fit)
                    elif(function == 'smooth'):
                        if(use_bounds == True):
                            popt,pcov = curve_fit(smooth_function,separation_bins,
                                corr_fit,sigma=dcorr_fit,bounds=bounds)
                        else:
                            popt,pcov = curve_fit(linear_function,separation_bins,
                                np.log(y_fit),sigma=dcorr_fit/corr_fit)

                    if(function == 'singletrunc'):
                        if(use_bounds == True):
                            popt,pcov = curve_fit(linear_truncation,separation_bins,
                                np.log(y_fit),sigma=dcorr_fit/corr_fit,bounds=bounds)
                        else:
                            popt,pcov = curve_fit(linear_truncation,separation_bins,
                                np.log(y_fit),sigma=dcorr_fit/corr_fit)
                except :
                    continue

                fit_bootstraps[i] = popt

            self.fit_values = np.median(fit_bootstraps,axis=0)
            self.fit_errors = np.std(fit_bootstraps,axis=0)
            self.fit_distribution = fit_bootstraps

        elif(method=='mcmc'):
            fit_MCMC(self,function=function,omega1=True,age=age,
                plot=plot)



    def get_data(self):
        import urllib.request as request
        import urllib.error
        from CombineCatalogs import Combine_NGC0628,Combine_NGC1313
        """
        Get catalogs and images.
        """

        #1. Catalogs

        # First check if required files already present
        # By default we use SB extinction with Average Aperture Correction
        if (not os.path.isfile(self.catalog_file)):
            print("Catalog file not present in ../data, downloading from LEGUS website.....")
            # These galaxies have multiple pointings whose catalogs need 
            # combining, or have special naming crieria            
            if(self.name in ['NGC_0628','NGC_1313','NGC_5194','NGC_7793']):
                
                #Multiple pointings
                if(self.name == 'NGC_1313'):
                    url = "https://archive.stsci.edu/hlsps/legus/ngc1313-e/\
cluster_catalogs/deterministic/hlsp_legus_hst_wfc3_ngc1313-e_multiband_v1_pa\
dagb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-e.tab'.format(self.name))
                    url = "https://archive.stsci.edu/hlsps/legus/ngc1313-w/\
cluster_catalogs/deterministic/hlsp_legus_hst_wfc3_ngc1313-w_multiband_v1_pa\
dagb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-w.tab'.format(self.name))
                    Combine_NGC1313()
                    os.remove('../data/{}-e.tab'.format(self.name))
                    os.remove('../data/{}-w.tab'.format(self.name))
                    # Consider the case where the catalog file name in info file is different
                    os.rename('../data/{}.tab'.format(self.name),self.catalog_file)

                #Multiple pointings
                elif(self.name == 'NGC_0628'):
                    url = "https://archive.stsci.edu/hlsps/legus/ngc628-c/clu\
ster_catalogs/deterministic/hlsp_legus_hst_acs-wfc3_ngc628-c_multiband_v1_pada\
gb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-c.tab'.format(self.name))
                    url = "https://archive.stsci.edu/hlsps/legus/ngc628-e/clu\
ster_catalogs/deterministic/hlsp_legus_hst_acs-wfc3_ngc628-e_multiband_v1_pada\
gb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-e.tab'.format(self.name))
                    Combine_NGC0628()
                    os.remove('../data/{}-c.tab'.format(self.name))
                    os.remove('../data/{}-e.tab'.format(self.name))
                    # Consider the case where the catalog file name in info file is different
                    os.rename('../data/{}.tab'.format(self.name),self.catalog_file)

                #Special name
                elif(self.name == 'NGC_5194'):
                    url = "https://archive.stsci.edu/hlsps/legus/mosaics/ngc5194_ng\
c5195/cluster_catalogs/deterministic/hlsp_legus_hst_acs-wfc3_ngc5194-ngc5195-mosaic_\
multiband_v1_padagb-sbext-avgapcor.tab"
                    request.urlretrieve(url,self.catalog_file)

                #Multiple pointings
                elif(self.name == 'NGC_7793'):
                    url = "https://archive.stsci.edu/hlsps/legus/ngc7793-e/cluster_cata\
logs/deterministic/hlsp_legus_hst_wfc3_ngc7793-e_multiband_v1_padagb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-e.tab'.format(self.name))

                    url = "https://archive.stsci.edu/hlsps/legus/ngc7793-w/cluster_cata\
logs/deterministic/hlsp_legus_hst_wfc3_ngc7793-w_multiband_v1_padagb-sbext-avgapcor.tab"
                    request.urlretrieve(url,'../data/{}-w.tab'.format(self.name))
                    Combine_NGC7793()
                    os.remove('../data/{}-e.tab'.format(self.name))
                    os.remove('../data/{}-w.tab'.format(self.name))
                    # Consider the case where the catalog file name in info file is different
                    os.rename('../data/{}.tab'.format(self.name),self.catalog_file)

            # Dwarfs have been observed with ACS as well, so different naming scheme
            elif(self.name in ['NGC_3738','NGC_4449','NGC_5253']):
                galaxy_id = self.name.split('_')[1]
                url = 'https://archive.stsci.edu/hlsps/legus/ngc{}/cluster_catalogs/det\
erministic/hlsp_legus_hst_acs-wfc3_ngc{}_multiband_v1_padagb-sbext-avgapcor.tab'.format(galaxy_id,galaxy_id)
                try:
                    response = request.urlretrieve(url,self.catalog_file)
                except urllib.error.HTTPError as exception:
                    raise myError("Catalog link does not work. Verify if provided link is correct: {}".format(url))

            #These catalogs not available in public website
            elif(self.name in ['NGC_5457','NGC_3627']):
                raise myError("This galaxy's catalog is not yet available on the public website.\
Please contact Sean Linden at UMass (Amherst) for access to the catalogs.")
            else:
                galaxy_id = self.name.split('_')[1]
                url = "https://archive.stsci.edu/hlsps/legus/ngc{}/\
cluster_catalogs/deterministic/hlsp_legus_hst_wfc3_ngc{}_multiband_v\
1_padagb-sbext-avgapcor.tab".format(galaxy_id,galaxy_id)
                try:
                    request.urlretrieve(url,self.catalog_file)
                except urllib.error.HTTPError as exception:
                    raise myError("The catalog file could not be downloaded. Check if url provided\
is correct. You may need to provide suitable conditions url names in get_data function or \
download data yourself with consistent naming schemes for new galaxies not in Menon et al 2021b. The\
 attempted link is : {}".format(url))

        if(not os.path.isfile(self.fits_file)):
            print("Fits mosaic file not present in ../data, downloading from LEGUS website.....")
            #Now download the fits files
            list_gals = np.array(list_of_galaxies)
            index_listgals = np.where(self.name == list_gals)[0][0]
            mosaicFileLink = mosaic_links[index_listgals]
            fileEnding = mosaicFileLink.split('/')[-1]

            try:
                request.urlretrieve(mosaicFileLink,self.fits_file)                            
            except urllib.error.HTTPError as exception:
                raise myError("The mosaic fits file could not be downloaded. Check if url provided\
is correct. You may need to provide suitable conditions url names in get_data function or \
download data yourself with consistent naming schemes for new galaxies not in Menon et al 2021b. The\
attempted link is : {}".format(mosaicFileLink))

        #3. Region file
        if(not os.path.isfile(self.region_file)):
            print("Region file not present. Creating now.....")
            self.create_region_file()

        #Safety check: Should not generally be ever triggered
        #Check catalog file exists else throw error
        if(not os.path.exists(self.catalog_file)) :
            raise myError("Catalog file " + self.catalog_file +
                " does not exist. This is required to compute the TPCF.")
        if(not os.path.exists(self.fits_file)) :
            raise myError("Fits file " + self.fits_file +
                " does not exist. This is required to compute the TPCF.")

        #Check region file exists else throw error
        if(not os.path.exists(self.region_file)) :
            raise myError("Region file" + self.region_file +
                " does not exist. This is required to compute the TPCF.")

        print("All data required for TPCF computation procured...")


    def compare_AIC(self,version='AICc',age=None):
        """
        Compare AIC of 3 models fitted and choose best fit model
        Parameters:
            version: string
                Version of AIC to use. Can be 'AIC','BIC',or 'AICc'
        """

        if(age == 'young'):
            maxlnprob_singlepl = self.maxlnprob_singlepl_y
            maxlnprob_piecewise = self.maxlnprob_piecewise_y
            maxlnprob_singletrunc = self.maxlnprob_singletrunc_y
            nsamples = np.size(np.where(self.ycorr>0.0))
        elif(age == 'old'):
            maxlnprob_singlepl = self.maxlnprob_singlepl_o
            maxlnprob_piecewise = self.maxlnprob_piecewise_o
            maxlnprob_singletrunc = self.maxlnprob_singletrunc_o
            nsamples = np.size(np.where(self.ocorr>0.0))
        else:
            maxlnprob_singlepl = self.maxlnprob_singlepl
            maxlnprob_piecewise = self.maxlnprob_piecewise
            maxlnprob_singletrunc = self.maxlnprob_singletrunc
            nsamples = np.size(np.where(self.corr>0.0))

        if(version == 'AIC'):
            AIC_single = 2*2.0 - 2*maxlnprob_singlepl
            AIC_piecewise = 2*4.0 - 2*maxlnprob_piecewise
            AIC_singletrunc = 2*3.0 - 2*2*maxlnprob_singletrunc
        elif(version == 'BIC'):
            AIC_single = np.log(nsamples)*2.0 - 2*maxlnprob_singlepl
            AIC_piecewise = np.log(nsamples)*4.0 - 2*maxlnprob_piecewise
            AIC_singletrunc = np.log(nsamples)*3.0 - 2*2*maxlnprob_singletrunc
        elif(version == 'AICc'):
            AIC_single = 2*2.0 - 2*maxlnprob_singlepl +\
             (2*2*(2+1))/(nsamples - 2 - 1)
            AIC_piecewise = 2*4.0 - 2*maxlnprob_piecewise +\
             (2*4*(4+1))/(nsamples - 4 - 1)
            AIC_singletrunc = 2*3.0 - 2*2*maxlnprob_singletrunc +\
             (2*3*(3+1))/(nsamples - 3 - 1)
        else:
            raise myError("AIC type not recognised.")

        galaxy_functions = ['singlepl','piecewise','singletrunc']
        galaxy_AIC = [AIC_single,AIC_piecewise,AIC_singletrunc]
        galaxy_function = galaxy_functions[np.argmin(galaxy_AIC)]

        if(age == 'young'):
            self.aic_singlepl_y = AIC_single
            self.aic_piecewise_y = AIC_piecewise
            self.aic_singletrunc_y = AIC_singletrunc
            self.best_model_y = galaxy_function
            if(self.best_model_y == 'singlepl'):
                self.best_fit_y = self.fit_values_singlepl_y
                self.best_fit_error_y = self.fit_errors_singlepl_y
            elif(self.best_model_y == 'piecewise'):
                self.best_fit_y = self.fit_values_piecewise_y
                self.best_fit_error_y = self.fit_errors_piecewise_y
            else:
                self.best_fit_y = self.fit_values_singletrunc_y
                self.best_fit_error_y = self.fit_errors_singletrunc_y
        elif(age == 'old'):
            self.aic_singlepl_o = AIC_single
            self.aic_piecewise_o = AIC_piecewise
            self.aic_singletrunc_o = AIC_singletrunc
            self.best_model_o = galaxy_function
            if(self.best_model_o == 'singlepl'):
                self.best_fit_o = self.fit_values_singlepl_o
                self.best_fit_error_o = self.fit_errors_singlepl_o
            elif(self.best_model_o == 'piecewise'):
                self.best_fit_o = self.fit_values_piecewise_o
                self.best_fit_error_o = self.fit_errors_piecewise_o
            else:
                self.best_fit_o = self.fit_values_singletrunc_o
                self.best_fit_error_o = self.fit_errors_singletrunc_o
        else:
            self.aic_singlepl = AIC_single
            self.aic_piecewise = AIC_piecewise
            self.aic_singletrunc = AIC_singletrunc
            self.best_model = galaxy_function
            if(self.best_model == 'singlepl'):
                self.best_fit = self.fit_values_singlepl
                self.best_fit_error = self.fit_errors_singlepl
            elif(self.best_model == 'piecewise'):
                self.best_fit = self.fit_values_piecewise
                self.best_fit_error = self.fit_errors_piecewise
            else:
                self.best_fit = self.fit_values_singletrunc
                self.best_fit_error = self.fit_errors_singletrunc



    def getPhysicalProps(self):
        """
        Get physical properties from the TPCF as discussed in Section 4.4 of
        Menon et al 2021b
        """

        # Get alpha and fractal dimension
        if(self.best_model_y == 'singlepl'):
            self.alpha = self.median_fit_singlepl_y[1], \
            self.median_fit_singlepl_y[1]-self.low_error_singlepl_y[1], \
            self.high_error_singlepl_y[1]- self.median_fit_singlepl_y[1]
        elif(self.best_model_y == 'piecewise'):
            self.alpha = self.median_fit_piecewise_y[1], \
            self.median_fit_piecewise_y[1]-self.low_error_piecewise_y[1], \
            self.high_error_piecewise_y[1]- self.median_fit_piecewise_y[1]
        else:
            self.alpha = self.median_fit_singletrunc_y[1], \
            self.median_fit_singletrunc_y[1]-self.low_error_singletrunc_y[1], \
            self.high_error_singletrunc_y[1]- self.median_fit_singletrunc_y[1]

        self.D2 = 2 + self.alpha[0],self.alpha[1],self.alpha[2]


        # Now get lcorr
        if(self.best_model_y == 'singlepl'):
            indices = np.where(self.corr>0.0)
            separation_bins = self.bin_centres*(1./arcsec_to_degree)
            separation_bins = separation_bins.astype(np.float)
            separation_bins = separation_bins[indices].astype(np.float)
            cutoff_scale = np.max(separation_bins)
            cutoff_error = np.diff(separation_bins)[-1]/2.
        
            samples = np.random.uniform(cutoff_scale,cutoff_error,10000)
        elif(self.best_model_y == 'piecewise'):
            samples = np.exp(self.beta_distribution)
        else:
            if(self.name == 'NGC_5253'):
                samples = np.exp(self.beta_distribution)
            else:
                samples = np.zeros(10000)

        distance = self.distance*const.Parsec*1.e6
        distance_err = self.errordist*const.Parsec*1.e6
        distance_samples = np.random.normal(distance,distance_err,
            size=np.size(samples))
        arcsec_to_pc = u.arcsec.to(u.radian)*distance_samples/(const.Parsec)
        samples = samples*arcsec_to_pc

        lcorr_median_fit = np.percentile(samples,50,axis=0)
        lcorr_error_low = np.percentile(samples,16,axis=0)
        lcorr_error_low = lcorr_median_fit-lcorr_error_low
        lcorr_error_high = np.percentile(samples,84,axis=0)
        lcorr_error_high = lcorr_error_high - lcorr_median_fit
        lcorr = lcorr_median_fit, lcorr_error_low, lcorr_error_high
        # NGC 3738 was found to not show any excess fractal structure
        if(self.name == 'NGC_3738'):
            self.lcorr = None 
        else:
            self.lcorr = lcorr


        # Now get rc

        # 5253 has no old TPCF
        if(self.name == 'NGC_5253'):
            corr = self.ycorr 
            dcorr = self.ydcorr

        else:
            corr = self.ocorr 
            dcorr = self.odcorr
        indices = np.where(np.logical_and(corr>0.0,corr>dcorr))
        corr_fit = corr[indices].astype(np.float)
        nsamples = np.size(corr_fit)   

        if(nsamples>5):
            samples = self.theta_distribution
            distance_samples = np.random.normal(distance,distance_err,
                size=np.size(samples))
            arcsec_to_pc = u.arcsec.to(u.radian)*distance_samples/(const.Parsec)
            samples = samples * arcsec_to_pc
            rc_median_fit = np.percentile(samples,50,axis=0)
            rc_error_low = np.percentile(samples,16,axis=0)
            rc_error_low = rc_median_fit-rc_error_low
            rc_error_high = np.percentile(samples,84,axis=0)
            rc_error_high = rc_error_high - rc_median_fit
            rc = rc_median_fit, rc_error_low, rc_error_high
            self.rc_disk = rc
        else:
            self.rc_disk = None


class myError(Exception):
    """
    Class derived from Exception to handle exceptions raised by
    program-specific errors.

    Parameters
       message : string
          the error message
    """

    def __init__(self, message):
        Exception.__init__(self, message)
        self.message = message