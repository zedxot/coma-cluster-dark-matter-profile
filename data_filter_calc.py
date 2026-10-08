import pandas as pd
import numpy as np

print("Reading csv...")
df = pd.read_csv("~/Projects/physics/DESI/clean_project/galaxy_data_clean.csv")

ra_center_0 = np.radians(194.953)
dec_center_0 = np.radians(27.981)
ra_i = np.radians(df['TARGET_RA'])
dec_i = np.radians(df['TARGET_DEC'])

print("Calculating angular separation")
cos_theta = np.sin(dec_i) * np.sin(dec_center_0) + np.cos(dec_i) * np.cos(dec_center_0) * np.cos(ra_i - ra_center_0)
theta_rad = np.arccos(np.clip(cos_theta,-1.0,1.0))
df['ANGLR_SEP_RAD'] = theta_rad
df['ANGLR_SEP_DEG'] = np.degrees(theta_rad)

print("applying angular and redshift cut")
spatial_n_redshift_cut = (
    (df['ANGLR_SEP_DEG'] <= 2.0) &
    (df['Z'].between(0.013, 0.033))
)

coma_df = df[spatial_n_redshift_cut].copy().reset_index(drop=True)
print(len(coma_df))

output = "~/Projects/physics/DESI/coma_galaxy_data.csv"
coma_df.to_csv(output, index=False)
