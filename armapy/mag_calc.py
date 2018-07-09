#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 06:52:28 2017

@author: Patrick Rauer
"""

import os
import warnings
import math
import numpy as np
from astropy.table import Table

from armapy import svo

try:
    from scipy import interpolate
except ImportError:
    interpolate = None

# from . import extinction as extlaws

c = 2.99792458e18  # speed of light in Angstrom/sec

ab_zero = -48.60  # magnitude zero points
st_zero = -21.10

vega_file = os.path.join(os.path.dirname(__file__), r'alpha_lyr_stis_006.fits')


def get_vega_data(file):
    """
    Load Vega spectrum

    :param file: The path to the vega spectrum.
    :type file: str
    :returns: The wavelength and the flux from the vega spectrum.
    :rtype: numpy.array, numpy.array
    """
    vega = Table.read(file)
    return np.array(vega['WAVELENGTH']), np.array(vega['FLUX'])

try:
    vega_wavelength, vega_flux = get_vega_data(vega_file)
except KeyError:
    vega_wavelength = None
    vega_flux = None

filter_dir = os.path.join(os.path.dirname(__file__), 'filters')
list_filters = {}
show_filters = {}


def _list_files(dir_name):
    """
    Create dictionary containing paths to available filters
    """
    for root, dirs, files in os.walk(dir_name):
        f_key = os.path.basename(root)
        f = []
        for name in files:
            key = os.path.splitext(name)[0]
            list_filters[key] = os.path.join(root, name)
            f.append(key)
        if f:
            show_filters[f_key] = f


_list_files(filter_dir)

# Initialize list of available extinction laws - only Cardelli is implemented at the moment
# laws_dic = {'cardelli': extlaws.cardelli}
# list_laws = laws_dic.keys()


def _file_exists(name):
    """
    Check if a file exists and is accessible
    """
    return os.path.exists(name)


class Spectrum(object):
    """
    General spectrum class
    """

    def __init__(self):
        self.wavelength = None
        self.flux = None

    def load_fits(self, filename=None, wave_name='Wavelength', flux_name='Flux'):
        """
        Load a FITS file containing a spectrum
        """
        tab = Table.read(filename)
        self.wavelength = tab[wave_name]
        self.flux = tab[flux_name]

    def load_ascii(self, filename=None):
        """
        Load an ASCII file containing one or two columns

        :param filename: The path to the file
        :type filename: str
        """
        data = np.loadtxt(filename)
        if len(data.shape) == 1:
            self.flux = data
        elif len(data.shape) == 2:
            self.wavelength = data[:, 0]
            self.flux = data[:, 1]

    def set_flux(self, flux):
        """
        Sets a new flux range

        :param flux: The new flux
        :type flux: numpy.array
        """
        self.flux = flux

    def set_wavelength(self, w):
        """
        Sets a new wavelength range

        :param w: The new wavelengths
        :type w: numpy.array
        """
        self.wavelength = w


def _wave_range(wb, ws):
    """
    Checks if the wavelength of the band object is complete covered by the spectra wavelength range

    :param wb: The wavelengths of the BandPass object
    :type wb: numpy.array
    :param ws: THe wavelength of the spectra
    :type ws: numpy.array
    :return: The indices of the covered wavelength range
    """
    wb_max = np.max(wb)
    wb_min = np.min(wb)
    ws_max = np.max(ws)
    ws_min = np.min(ws)
    if wb_min < ws_min or wb_max > ws_max:
        warnings.warn('Spectrum doesn\'t overlap the complete bandpass')
        return np.linspace(0, len(ws)-1, num=len(ws), dtype=np.int32)
    else:
        wave_range = np.where(np.logical_and(ws <= wb_max, ws >= wb_min))[0]
        return wave_range


class StarSpectrum(Spectrum):
    """
    StarSpectrum contains the wavelength and the flux of a source.
    It can load the data from a file directly or you can give the wavelengths and the fluxes as initial arguments.
    Together with :class:`Band` you can calculate the magnitudes from the spectra for different filters

    :param file:
        indicates the name of a file with the spectrum to load.
        It accepts files with the FITS format or plain ASCII
        with two column. If it is a fits file the wavelength-column must have the name 'Wavelength' and
        the flux-column must have the name 'Flux'.
        If it is a ascii file, the first column must be the wavelength values and the second column must be
        the flux values.
        Default is None
    :type file: str
    :param wavelength:
        The wavelength as an array if no file is given, or if the file name is given,
        the name of the wavelength column
    :type wavelength: numpy.array
    :param flux: The flux as an array if no file is given, or if the file name is given,
        the name of the flux column
    :type flux: numpy.array
    """

    def __init__(self, file=None, wavelength=None, flux=None):
        """

        """
        Spectrum.__init__(self)
        self.wavelength = wavelength
        self.flux = flux
        self._file = file
        self.__load__()

    def __load__(self):
        if self._file:
            if os.path.exists(self._file):
                f_name, f_extension = os.path.splitext(self._file)
                if '.fit' in f_extension:
                    if type(self.wavelength) is str and type(self.flux) is str:
                        self.load_fits(self._file, wave_name=self.wavelength, flux_name=self.flux)
                    else:
                        self.load_fits(self._file)
                else:
                    self.load_ascii(self._file)
            else:
                warnings.warn("Warning: Could not find file {} - no spectrum loaded" .format(self._file))

    def ap_mag(self, band, mag='Vega', mag_zero=0.):
        """
        Compute an apparent magnitude in a given system

        :param band:
            function or callable
            a band pass given as a function of wavelength
        :type band: BandPass
        :param mag:
            {'Vega','AB','ST'}
            name of the magnitude system to use
            Default is 'Vega'
        :type mag: str
        :param mag_zero:
            magnitude of the standard star corresponding to
            the zero point (usually 0.) in Vega system
            Default is 0.
        :type mag_zero: float
        :returns: The magnitude in the filter and in the magnitude system
        :rtype float:
        """
        # estimates the needed range
        r = _wave_range(band.wavelength, self.wavelength)
        wr = self.wavelength[r]
        fr = self.flux[r]
        band_interp = band.get_transmission(wr)
        band_interp_wr = band_interp * wr
        integrated_flux = np.trapz(fr * band_interp_wr, x=wr)
        if mag == 'AB':
            f = integrated_flux / np.trapz(band_interp / wr, x=wr) * c
            return -2.5 * math.log10(f) + ab_zero
        elif mag == 'ST':
            f = integrated_flux / np.trapz(band_interp_wr, x=wr)
            return -2.5 * math.log10(f) + st_zero
        elif mag == 'Vega':
            # TODO an error that an ndarray is not callable
            vega_r = _wave_range(band.wavelength, vega_wavelength)
            vega_wr = vega_wavelength[vega_r]
            vega_fr = vega_flux[vega_r]
            f = integrated_flux
            vega_f = np.trapz(vega_fr * band.get_transmission(vega_wr) * vega_wr, x=vega_wr)
            return -2.5 * math.log10(f) + 2.5 * np.log10(vega_f) + mag_zero
        else:
            raise ValueError('Magnitude system is not a valid choice, check input string')


class Band(Spectrum):
    """
    Band passes photometric response curves

    :param band:
        indicates the name of a filter in the internal list (see show_filters function)
        or a filename of a file containing two columns, one with the wavelength and
        one with the normalized response
    :type band: str
    :param smt:
        indicates the type of smoothing to perform on the response curve. Uses same
        names as Scipy's interp1d ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')
        Default is slinear
    :type smt: str

    """

    def __init__(self, band=None, wavelength=None, response=None, smt='slinear'):
        """
        Initialize a band

        """
        Spectrum.__init__(self)
        self.wavelength = wavelength
        self.flux = None
        self.name = None
        self.interpolated_response = None
        self._smt = smt
        self.response = response
        self.interpolated_response = None
        
        if band is not None:
            self._load_(band)
        elif self.wavelength is not None and self.response is not None:
            self.smooth(self._smt)

    def _load_(self, band):
        # Check if input is a band part of our list or an existing file
        if band:
            if self._smt in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
                if band in list_filters:
                    self.load_ascii(list_filters[band])

                elif _file_exists(band):
                    self.load_ascii(band)

                else:
                    warnings.warn("Warning: Could not find filter or file {} ".format(band))
                self.smooth(kind=self._smt)

            else:
                raise ValueError('Unknown type of smoothing')
        elif self.wavelength is not None and self.response is not None:
            self.smooth(self._smt)

    def smooth(self, kind='linear'):
        """
        Interpolate the filter curve

        :param kind: Kind of interpolation
        :type kind: str
        """
        if interpolate is not None:
            self.interpolated_response = interpolate.interp1d(self.wavelength, self.response, kind=kind,
                                                              bounds_error=False, fill_value=0.)
        else:
            warnings.warn("Warning: Smoothing requires Scipy, returning a linear interpolation instead.")
            self.interpolated_response = lambda w: np.interp(w, self.wavelength, self.response, left=0., right=0.)

    def __call__(self, spectra, mag_system='AB'):
        """
        Returns the magnitude of the spectra in this filter

        :param spectra: The spectra
        :type spectra: StarSpectrum
        :param mag_system: The magnitude system ('Vega', 'AB', 'ST')
        :type mag_system: str
        :return: The magnitude in the filter
        :rtype: float
        """
        return spectra.ap_mag(self, mag=mag_system)

    def get_transmission(self, w):
        """
        Returns the interpolated filter curve function for the wavelength range

        :param w: The wavelengths of the requested area
        :type w: StarSpectrum
        :return: The values of the interpolated filter curve for these wavelength
        :rtype: numpy.array
        """
        return self.interpolated_response(w)


class BandSVO(Band):
    telescope = ''
    instrument = ''
    filter_band = ''

    def __init__(self, telescope, instrument, filt, smt='linear'):
        """
        A child class of :class:`Band`. It does exactly the same but it will take
        the filter transmission curves from the `SVO <http://svo2.cab.inta-csic.es/svo/theory/fps3/>`_
        web-page directly.
    
        :param telescope: The name of the telescope/satellite
        :type telescope: str
        :param instrument: 
            The name of the instrument (if the telescope has only one instrument,
            the name of the instrument is mostly equal to the telescope name)
        :type instrument: str
        :param filt: The name of the band
        :type filt: str
        """
        Band.__init__(self, smt=smt)
        self.telescope = telescope
        self.instrument = instrument
        self.filter_band = filt
        filter_curve = svo.get_filter_curve(telescope, instrument, filt)
        self.information = svo.get_filter_information(telescope, instrument, filt)
        self.wavelength = filter_curve['Wavelength']
        self.response = filter_curve['Transmission']
        self.smooth(kind=self._smt)


class ColorSVO:
    """
    Class to calculate the photometric color of a spectra. It uses the filter curves from SVO.

    :param telescope: The name of the telescope/satellite
    :type telescope: str
    :param instrument: The name of the instrument
    :type instrument: str
    :param band1: The name of the first filter
    :type band1: str
    :param band2: The name of the second filter
    :type band2: str
    """
    def __init__(self, telescope, instrument, band1, band2):
        self.filter1 = BandSVO(telescope, instrument, band1)
        self.filter2 = BandSVO(telescope, instrument, band2)

    def __call__(self, spec, mag_system='AB'):
        return self.filter1(spec, mag_system)-self.filter2(spec, mag_system)
