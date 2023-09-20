"""
Download data needed for 
"""
import pathlib

import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.net.attr import or_
from sunpy.time import parse_time

import mocksipipeline.net  # Register XRTSynopticClient
from mocksipipeline.net import FilterWheel1, FilterWheel2


if __name__ == '__main__':
    obstime_window = 1 * u.h
    aia_wavelengths = [94, 131, 171, 193, 211, 335] * u.angstrom
    xrt_filters = [('Be-thin', 'Open')]

    obstime = parse_time(snakemake.config['obstime'])
    time = a.Time(obstime-obstime_window/2, end=obstime+obstime_window/2, near=obstime)
    # Construct AIA query
    aia_query = a.Instrument.aia & or_(*[a.Wavelength(w) for w in aia_wavelengths])
    # Construct XRT query
    fw_combos = [FilterWheel1(fw1) & FilterWheel2(fw2) for fw1, fw2 in xrt_filters]
    xrt_query = a.Instrument.xrt & a.Source.hinode & a.Provider('MSU') & a.Level(2) & or_(*fw_combos)
    query = Fido.search(time, aia_query | xrt_query)

    data_directory = pathlib.Path(snakemake.output[0]).parent
    files = Fido.fetch(query, path=data_directory / '{instrument}', overwrite=True)
