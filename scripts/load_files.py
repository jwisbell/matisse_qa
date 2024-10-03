"""
Author: Jacob Isbell 
Date: October 1, 2024

Script which contains functions to load MATISSE fits files and turn them into dictionaries for later usage.
This will unify data formatting.

Specifically
1. Loads TARGET_RAW_INT/CALIB_RAW_INT 
2. Loads OPD files
3. Loads photometric files
"""
from astropy.io import fits 


def load_raw_int(fnames):
    """
    If fname is a single file load it, if fname is an array of files, load each.
    In both cases, return a data_dictionary (or pandas dataframe?)
    """
    if fnames is str:
        fnames = [fnames]

    data = {'vis':{"cflux":[],"cflux_err":[],"u":[], "v":[],"wl_vis":[],"vis2":[], "vis2_err":[],'bcd':[],'vis2_sta':[] },\
            'cphase':{"t3phi":[],"t3phi_err":[], "u1":[],"u2":[], "v1":[],"v2":[],"wl_t3":[], "t3_sta":[], 'bcd':[]},\
            'phot':{"phot_flux":[], "phot_flux_err":[], "phot_tel":[], "phot_wl":[], 'bcd':[]},
            'inst':{'bcd':[] }}

    """
    files_oo = np.sort(glob(fdir+'/*_???_INT_0005.fits'))
        files_ii = np.sort(glob(fdir+'/*_???_INT_0006.fits'))
        files_io = np.sort(glob(fdir+'/*_???_INT_0005.fits'))
        files_oi = np.sort(glob(fdir+'/*_???_INT_0006.fits'))"""

    for f in fnames:
        x = fits.open(f)
        print(f'\t\t {f}')

        bcd = 'oo'
        if '0002.fits' in f:
            bcd = 'ii'
        elif '0003.fits' in f:
            bcd = 'io'
        elif '0004.fits' in f:
            bcd = 'oi'
        elif '0005.fits' in f:
            bcd = 'oo_phot'
        elif '0006.fits' in f:
            bcd = 'ii_phot'
        

        

        oi_key = "OI_VIS"
        #print(x[0].header.keys )
        for i,val in enumerate( x[oi_key].data["visamp"] ):
            data['vis']['cflux'].append(val)
            data['vis']["u"].append( x[oi_key].data["ucoord"][i] ) 
            data['vis']["v"].append( x[oi_key].data["vcoord"][i] ) 
            data['vis']["cflux_err"].append( x[oi_key].data["visamperr"][i] ) 
            data['vis']["wl_vis"].append( x['oi_wavelength'].data["eff_wave"] ) 
            data['vis']['vis2_sta'].append( x[oi_key].data["sta_index"][i] )
            data['vis']['bcd'].append( bcd )


        oi_key = "OI_VIS2"
        for i,val in enumerate( x[oi_key].data["vis2data"] ):
            data['vis']['vis2'].append(val)
            data['vis']["vis2_err"].append( x[oi_key].data["vis2err"][i] )

        oi_key = "OI_T3"
        for i,val in enumerate( x[oi_key].data["t3phi"] ):
            data['cphase']['t3phi'].append(val)
            data['cphase']["t3phi_err"].append( x[oi_key].data["t3phierr"][i] ) 
            data['cphase']["u1"].append( x[oi_key].data["u1coord"][i] ) 
            data['cphase']["v1"].append( x[oi_key].data["v1coord"][i] ) 
            data['cphase']["u2"].append( x[oi_key].data["u2coord"][i] ) 
            data['cphase']["v2"].append( x[oi_key].data["v2coord"][i] ) 
            data['cphase']["t3_sta"].append( x[oi_key].data["sta_index"][i] ) 
            data['cphase']["wl_t3"].append( x['oi_wavelength'].data["eff_wave"])
            data['cphase']['bcd'].append( bcd )

        oi_key = "OI_FLUX"
        try:
            for i,val in enumerate( x[oi_key].data["fluxdata"] ):
                data['phot']['phot_flux'].append(val)
                data['phot']["phot_flux_err"].append( x[oi_key].data["fluxerr"][i] ) 
                data['phot']["phot_tel"].append( x[oi_key].data["sta_index"][i] ) 
                data['phot']["phot_wl"].append( x['oi_wavelength'].data["eff_wave"])
                data['phot']['bcd'].append( bcd )
        except KeyError:
            print('No flux data found for', f)

    return data

def load_opd(fname):
    return None
