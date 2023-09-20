"""
Sample and map counts from spectral cube to detector plane
"""
import astropy.units as u

from mocksipipeline.util import read_data_cube
from mocksipipeline.detector.component import sample_and_remap_spectral_cube
from mocksipipeline.detector.response import Channel, SpectrogramChannel
from overlappy.io import write_overlappogram


if __name__ == '__main__':
    spectral_cube = read_data_cube(snakemake.input[0])
    if snakemake.params.channel_type == 'filtergram':
        channel = Channel(snakemake.params.channel_name, full_detector=False)
    elif snakemake.params.channel_type == 'dispersed':
        channel = SpectrogramChannel(snakemake.params.spectral_order, full_detector=False)
    else:
        raise ValueError(f'Unrecognized channel type {snakemake.params.channel_type}')
    detector_image = sample_and_remap_spectral_cube(
        spectral_cube,
        channel,
        dt=u.Quantity(snakemake.config['exposure_time'], 's'),
        interval=u.Quantity(snakemake.config['integration_time'], 's'),
        convert_to_dn=True,
        roll_angle=float(snakemake.config['roll_angle'])*u.deg,
        dispersion_angle=float(snakemake.config['dispersion_angle'])*u.deg,
        chunks='auto',
    )
    write_overlappogram(detector_image, snakemake.output[0])
