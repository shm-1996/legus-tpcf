from header import *
from Galaxy import Galaxy, myError


def Combine_NGC1313():
    
    print("Combining different mosaic pointing catalogs onto common one for NGC 1313.")

    #Read catalogs
    Catalog_w = '../data/NGC_1313-w.tab'
    Catalog_e = '../data/ngc_1313-e.tab'
    file_w = np.loadtxt(Catalog_w)
    file_e = np.loadtxt(Catalog_e)

    ra_w,dec_w = file_w[:,5],file_w[:,6]
    ra_e,dec_e = file_e[:,5],file_e[:,6]

    #Get intersection
    ra,ra_inde,ra_indw = np.intersect1d(ra_e,ra_w,return_indices=True)
    dec,dec_inde,dec_indw = np.intersect1d(dec_e,dec_w,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_inde,dec_inde)

    #Remove duplicates from east catalog
    file_e_updated = np.delete(file_e,Duplicate_indices,axis=0)

    #Remove 3rd and 4th column to maintain consistency with other galaxy
    # table formats
    file_e_updated = np.delete(file_e_updated,[3,4],axis=1)
    file_w = np.delete(file_w,[3,4],axis=1)

    Combined_file = np.append(file_e_updated,file_w,axis=0)

    np.savetxt('../data/NGC_1313.tab',
           Combined_file)

def Combine_NGC0628():
    print("Combining different mosaic pointing catalogs onto common one for NGC 0628.")

    #Read catalogs
    Catalog_c = '../data/NGC_0628-c.tab'
    Catalog_e = '../data/ngc_0628-e.tab'
    file_c = np.loadtxt(Catalog_c)
    file_e = np.loadtxt(Catalog_e)

    ra_c,dec_c = file_c[:,3],file_c[:,4]
    ra_e,dec_e = file_e[:,3],file_e[:,4]

    #Get intersection
    ra,ra_inde,ra_indc = np.intersect1d(ra_e,ra_c,return_indices=True)
    dec,dec_inde,dec_indc = np.intersect1d(dec_e,dec_c,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_inde,dec_inde)

    #Remove duplicates from east catalog
    file_e_updated = np.delete(file_e,Duplicate_indices,axis=0)

    Combined_file = np.append(file_e_updated,file_c,axis=0)

    np.savetxt('../data/NGC_0628.tab',
           Combined_file)

def Combine_NGC7793():
    print("Combining different mosaic pointing catalogs onto common one for NGC 7793.")

    #Read catalogs
    Catalog_c = '../data/NGC_7793-c.tab'
    Catalog_e = '../data/ngc_7793-e.tab'
    file_c = np.loadtxt(Catalog_c)
    file_e = np.loadtxt(Catalog_e)

    ra_c,dec_c = file_c[:,3],file_c[:,4]
    ra_e,dec_e = file_e[:,3],file_e[:,4]

    #Get intersection
    ra,ra_inde,ra_indc = np.intersect1d(ra_e,ra_c,return_indices=True)
    dec,dec_inde,dec_indc = np.intersect1d(dec_e,dec_c,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_inde,dec_inde)

    #Remove duplicates from east catalog
    file_e_updated = np.delete(file_e,Duplicate_indices,axis=0)

    Combined_file = np.append(file_e_updated,file_c,axis=0)

    np.savetxt('../data/NGC_7793.tab',
           Combined_file)



def Combine_NGC5457():
    
    print("Combining different mosaic pointing catalogs onto common one for NGC 5457.")

    #Read Catalogs
    Catalog_c = '../data/NGC_5457_Catalogs/NGC_5457_c.tab'
    Catalog_nw1 = '../data/NGC_5457_Catalogs/NGC_5457_nw1.tab'
    Catalog_nw2 = '../data/NGC_5457_Catalogs/NGC_5457_nw2.tab'
    Catalog_nw3 = '../data/NGC_5457_Catalogs/NGC_5457_nw3.tab'
    Catalog_se = '../data/NGC_5457_Catalogs/NGC_5457_se.tab'

    file_c = np.loadtxt(Catalog_c)
    file_nw1 = np.loadtxt(Catalog_nw1)
    file_nw2 = np.loadtxt(Catalog_nw2)
    file_nw3 = np.loadtxt(Catalog_nw3)
    file_se = np.loadtxt(Catalog_se)

    #Northwest-pointing has different no of columns (37) vs others (34)
    file_nw1 = np.delete(file_nw1,[34,35,36],axis=1)
    file_nw2 = np.delete(file_nw2,[34,35,36],axis=1)
    file_nw3 = np.delete(file_nw3,[34,35,36],axis=1)


    #Logic: Combine iteratively, removing duplicates with each combine operation
    
    #First set combined file as only centre pointing
    Combined_file = file_c
    
    #Combine this with NW1
    ra_combined,dec_combined = Combined_file[:,3],Combined_file[:,4]
    file_new = file_nw1
    ra_new,dec_new = file_new[:,3],file_new[:,4]
    #Get intersection
    ra,ra_1,ra_2 = np.intersect1d(ra_combined,ra_new,return_indices=True)
    dec,dec_1,dec_2 = np.intersect1d(dec_combined,dec_new,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_2,dec_2)
    file_new_updated = np.delete(file_new,Duplicate_indices,axis=0)
    #Combine
    Combined_file = np.append(Combined_file,file_new_updated,axis=0)

    #Combine this with NW2
    ra_combined,dec_combined = Combined_file[:,3],Combined_file[:,4]
    file_new = file_nw2
    ra_new,dec_new = file_new[:,3],file_new[:,4]
    #Get intersection
    ra,ra_1,ra_2 = np.intersect1d(ra_combined,ra_new,return_indices=True)
    dec,dec_1,dec_2 = np.intersect1d(dec_combined,dec_new,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_2,dec_2)
    file_new_updated = np.delete(file_new,Duplicate_indices,axis=0)
    #Combine
    Combined_file = np.append(Combined_file,file_new_updated,axis=0)

    #Combine this with NW3
    ra_combined,dec_combined = Combined_file[:,3],Combined_file[:,4]
    file_new = file_nw3
    ra_new,dec_new = file_new[:,3],file_new[:,4]
    #Get intersection
    ra,ra_1,ra_2 = np.intersect1d(ra_combined,ra_new,return_indices=True)
    dec,dec_1,dec_2 = np.intersect1d(dec_combined,dec_new,return_indices=True)
    Duplicate_indices = np.intersect1d(ra_2,dec_2)
    file_new_updated = np.delete(file_new,Duplicate_indices,axis=0)
    #Combine
    Combined_file = np.append(Combined_file,file_new_updated,axis=0)


    #Combine this with SE
    # ra_combined,dec_combined = Combined_file[:,3],Combined_file[:,4]
    # file_new = file_se
    # ra_new,dec_new = file_new[:,3],file_new[:,4]
    # #Get intersection
    # ra,ra_1,ra_2 = np.intersect1d(ra_combined,ra_new,return_indices=True)
    # dec,dec_1,dec_2 = np.intersect1d(dec_combined,dec_new,return_indices=True)
    # Duplicate_indices = np.intersect1d(ra_2,dec_2)
    # file_new_updated = np.delete(file_new,Duplicate_indices,axis=0)
    # #Combine
    # Combined_file = np.append(Combined_file,file_new_updated,axis=0)


    np.savetxt('../data/NGC_5457.tab',
           Combined_file)
