from astropy.io import fits
import pandas as pd
import numpy as np

fits_file = "~/Projects/physics/DESI/spectro/redux/iron/zcatalog/v1/zall-pix-iron.fits"

with fits.open(fits_file, memmap=True) as hdu_list:
    hdu = hdu_list[1]
    
    type_mask = (
        (hdu.data['SPECTYPE'] == 'GALAXY') &
        (hdu.data['ZWARN'] == 0) &
        (hdu.data['ZCAT_PRIMARY'] == True)
    )
    selected_cols = [
    'TARGETID',
    'SURVEY',
    'ZCAT_PRIMARY',
    'TARGET_RA',
    'TARGET_DEC',
    'Z',
    'ZERR',
    'ZWARN',
    'DELTACHI2',
    'SPECTYPE',
    'SUBTYPE',
    'FLUX_G',
    'FLUX_R',
    'FLUX_Z',
    'EBV',
    'MORPHTYPE',
    'SERSIC',
    'SHAPE_R'
]
    
    data_dict = {
        col: hdu.data[col][type_mask] 
        for col in selected_cols
    }

df = pd.DataFrame(data_dict)

c = 299792.458
z = df['Z'].to_numpy(dtype=np.float64)
df['LOS_V'] = c * (((1.0 + z)**2 - 1.0)/((1.0 + z)**2 + 1.0))
df['V_ERR'] = c * df['ZERR'].to_numpy(dtype=np.float64)
output = "~/Projects/physics/DESI/galaxy_data_clean.csv"
df.to_csv(output, index=False)
