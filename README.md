# legus-tpcf
This repository performs the Two-Point Correlation Function (TPCF) Analysis of Star Clusters with the LEGUS Survey, and can reproduce the analysis and plots of
Menon et al 2021b (submitted). I have also added 'summary' class objects for 
each galaxy in `results/Menonetal2021/*.pkl` that contain the final results and values I obtained and report in the paper. 
These files are stored with [Git LFS storage](https://git-lfs.github.com) so do not worry about storing them! 
One could also reproduce the analysis for a specific galaxy, or all galaxies. The catalog files and HST images required for them would automatically be downloaded
from the LEGUS website. Note: The analysis for NGC 3627 and NGC 5457 cannot be reproduced as their star cluster catalogs are not public yet. 

More importantly, the tools developed/written are quite flexible and could be easily extended to more galaxies whenever catalogs are available for them.
Clone this repository and feel free to use it for further analysis. 
Let me know if you have any questions, find any issues, or have further comments by emailing me at shyam.menon@anu.edu.au. 

## Usage

### Reproducing the Plots of Menon et al 2021b
To just reproduce the plots of the paper perform the following steps - 
1. `cd src/`
2. `python Menon21_Plots.py`
This should create the figures and store them in `results/`. This uses the pre-obtained and stored results of the class objects stored in `results/Menonetal2021/*.pkl`. 
You could explore them as well by loading the pickle file and calling `help(file)` to see the list of attributes and methods available for them. 

### Reproducing the Analysis
To run the entire analysis pipeline for all the galaxies perform the following steps - 
1. `cd src/`
2. `python main.py -galaxy all`
3. `python Menon21_Plots.py -indir ../results/` to then plot the Figures in the paper with them. 

To run the pipeline only for a specific galaxy, say NGC 0628, modify step 2 to `python main.py -galaxy NGC_0628`

### Performing analysis for a galaxy outside of the ones in Menon et al 2021b
To perform the analysis for a new galaxy, one would need a star cluster catalog file and a HST image with WCS information of the field-of-view of the observation. 
If these are available and stored in `data/` one can perform the analysis for the galaxy, say NGC 0000 by - 
1. Following the instructions provided in the [infoReadme](data/galaxyinfo/infoReadme) file to provide some metadata for the galaxy. 
2. Open a jupyter notebook or ipython and perform 
```
from Galaxy import *
galaxy_class = Galaxy('NGC_0000')
galaxy_class.Compute_TPCF(age='young') #for young clusters
galaxy_class.Compute_TPCF(age='old') #for old clusters
```
3. For fitting functional forms and further analysis go through the code described in `main.py` to obtain an idea. 


### The gist
I use an object-oriented approach to perform the analysis. Each galaxy is initialised as an object of class `Galaxy` whose attributes 
and methods are defined in [Galaxy.py](src/Galaxy.py). Plotting is carried out with the class objects of type `myPlot` defined in [Plot_Class.py](src/Plot_Class.py).
[TPCF.py](src/TPCF.py) contains the routines for computing the TPCF given a pair of points in space (cartesian or equatorial coordinates). 
[MCMCfit.py](src/MCMCfit.py) contains the routines to fit the various functional models described in the paper. 
All data such as catalog files and HST image files are stored in `data/` with other metadata of galaxies stored in `data/galaxyinfo`.
Resulting Galaxy object files from the analysis would be saved in `results/`. 

## Acknowledgements
We use modified versions of the following open-source software in our implementation - 
1. [FootprintFinder](http://hla.stsci.edu/Footprintfinder/FootprintFinder.html) developed by StSci
2. [astropy](https://github.com/astropy/astropy) 
3. [astroML](http://github.com/astroML/astroML) for pair computatations. Citation: [Vanderplas et al 2012](https://ieeexplore.ieee.org/document/6382200). 
4. [CMasher](https://github.com/1313e/CMasher) for colorbars. 
5. [regions](https://github.com/astropy/regions)
6. [emcee](https://github.com/dfm/emcee) for fitting with an MCMC. 
7. [corner](https://github.com/dfm/corner.py) for plotting corner plots of the posterior distribution. 
