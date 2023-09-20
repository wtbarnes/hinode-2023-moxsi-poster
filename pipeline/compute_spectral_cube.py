"""
Correct and align data
"""
import pathlib

import aiapy
import astropy.time
import astropy.units as u
import numpy as np
import sunpy.io._fits as sunpy_fits
import sunpy.map

from mocksipipeline.physics.dem.data_prep import DataPrep
from mocksipipeline.physics.spectral import SpectralModel


if __name__ == '__main__':
    obstime = astropy.time.Time(snakemake.config['obstime'])
    temperature_bin_edges = 10**np.arange(
        float(snakemake.config['log_t_left_edge']),
        float(snakemake.config['log_t_right_edge']) + float(snakemake.config['delta_log_t']),
        float(snakemake.config['delta_log_t']),
    ) * u.K

    aia_maps = sunpy.map.Map(list(pathlib.Path(snakemake.input[0]).glob('*.fits')))
    xrt_maps = sunpy.map.Map(list(pathlib.Path(snakemake.input[1]).glob('*.fits')))
    if isinstance(xrt_maps, sunpy.map.GenericMap):
        xrt_maps = [xrt_maps,]

    pointing_table_ar = aiapy.calibrate.util.get_pointing_table(obstime-6*u.h, obstime+6*u.h)
    correction_table = aiapy.calibrate.util.get_correction_table()
    error_table = aiapy.calibrate.util.get_error_table()

    dq = DataPrep(map_list=aia_maps+xrt_maps,
                  aia_error_table=error_table,
                  aia_correction_table=correction_table,
                  aia_pointing_table=pointing_table_ar,
                  temperature_bin_edges=temperature_bin_edges)
    dem = dq.run()

    spectral_model = SpectralModel(spectral_table=snakemake.config['spectral_table'])
    spectral_cube = spectral_model.run(dem, dq.celestial_wcs)
    sunpy_fits.write(snakemake.output[0],
                     spectral_cube.data,
                     spectral_cube.meta,
                     overwrite=True)
