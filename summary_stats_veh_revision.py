import os
import numpy as np
import pandas as pd
import sys
import glob

in_dir = '/home/atom/ongoing/work_Veh/specific_thinning_areas/dh_pergla_cut_corrected_withmissing'
list_csv = glob.glob(os.path.join(in_dir, '*_rates.csv'))

list_df = []
for csv in list_csv:
    list_df.append(pd.read_csv(csv))

df_all = pd.concat(list_df)

df_all = df_all[df_all.period=='2000-01-01_2020-01-01']

# Mean statistics, area-weighted, for all glacier dams
mean_valid_obs = np.sum(df_all.valid_obs * df_all.area) / np.sum(df_all.area)
mean_coverage =  np.sum(df_all.perc_area_meas * df_all.area) / np.sum(df_all.area)
mean_dhdt = np.sum(df_all.dhdt * df_all.area) / np.sum(df_all.area)
mean_err_dhdt = 2 * np.sum(df_all.err_dhdt * df_all.area) / np.sum(df_all.area)

# Minimum statistics
min_valid_obs = np.min(df_all.valid_obs)
min_coverage = np.min(df_all.perc_area_meas)

print('SUMMARY STATISTICS')

print('Average number of observations in 20 years: {:.1f} '.format(mean_valid_obs))
print('Average spatial coverage over glacier dams: {:.1f}%'.format(mean_coverage*100))
print('Average thinning rate over glacier dams: {:.12} m yr-1'.format(mean_dhdt))
print('Average error at the 95% confidence level: {:.2f} m yr-1'.format(mean_err_dhdt))

print('Minimum number of observations in 20 years: {:.1f}'.format(min_valid_obs))
print('Minimum spatial coverage over glacier dams: {:.1f}%'.format(min_coverage*100))

print('SUPPLEMENTARY TABLE')

df_all = df_all[['reg', 'rgiid', 'area', 'dhdt', 'err_dhdt', 'perc_area_meas', 'valid_obs']]
regs = [r[6:8] for r in df_all.rgiid]
df_all['reg'] = regs
df_all['area'] /= 1000000
df_all['err_dhdt'] *= 2
df_all['perc_area_meas'] *= 100
covered_by_arcticdem = df_all.reg.isin(['01', '06', '08'])
df_all['arcticdem'] = covered_by_arcticdem
df_all = df_all.round({'area':2, 'dhdt': 2, 'err_dhdt': 2, 'perc_area_meas': 1, 'valid_obs': 1, 'reg': 0})
df_all.columns = ["Region", "RGIId", "Area (km²)", "Mean elevation change (m yr-1)", "Error (2-sigma) in mean elevation change (m yr-1)", "Area measured (%)", "Average valid elevations in 20 years", "Covered by ArcticDEM"]
df_all = df_all.sort_values(by=['Region', 'RGIId'])
df_all.to_csv('/home/atom/ongoing/work_Veh/specific_thinning_areas/veh_2022_supp_table_dem.csv', index=None)
