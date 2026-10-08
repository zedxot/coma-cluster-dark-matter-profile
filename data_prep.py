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
    # Identifiers & Primary Flags
    'TARGETID',         # Index 0   (K)  - Unique target identifier
    'SURVEY',           # Index 1   (7A) - Survey program ('main', 'sv1', etc.)
    'ZCAT_PRIMARY',     # Index 134 (L)  - Best redshift flag (eliminates duplicates)
    
    # Sky Coordinates
    'TARGET_RA',        # Index 16  (D)  - Right Ascension (degrees)
    'TARGET_DEC',       # Index 17  (D)  - Declination (degrees)
    
    # Redshifts & Spectroscopic Quality
    'Z',                # Index 5   (D)  - Spectroscopic redshift
    'ZERR',             # Index 6   (D)  - Redshift measurement error (for error propagation)
    'ZWARN',            # Index 7   (K)  - Warning bitmask (0 = reliable fit)
    'DELTACHI2',        # Index 14  (D)  - Fit significance over second-best template
    'SPECTYPE',         # Index 11  (6A) - Spectral classification ('GALAXY')
    'SUBTYPE',          # Index 12  (20A)- Sub-classification / spectral details
    
    # Optical Photometry & Extinction (Luminosities & Mass-to-Light Ratios)
    'FLUX_G',           # Index 32  (E)  - Total g-band flux (nanomaggies)
    'FLUX_R',           # Index 33  (E)  - Total r-band flux (nanomaggies)
    'FLUX_Z',           # Index 34  (E)  - Total z-band flux (nanomaggies)
    'EBV',              # Index 31  (E)  - Galactic extinction E(B-V) for de-reddening
    
    # Morphological Properties (Light Profile Tracers)
    'MORPHTYPE',        # Index 30  (4A) - Morphological type (PSF, REX, EXP, DEV, COMP)
    'SERSIC',           # Index 49  (E)  - Sersic index n
    'SHAPE_R'           # Index 50  (E)  - Half-light radius in arcsec (for galaxy light profile)
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
