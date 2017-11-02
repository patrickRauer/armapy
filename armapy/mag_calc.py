#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 06:52:28 2017

@author: Patrick Rauer
"""

import os
import warnings

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
    """
    vega = Table.read(file)
    print(vega)
    return vega['WAVELENGTH'], vega['FLUX']

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


class GeneralSpectrum(object):
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
        return np.linspace(0, len(ws)-1, num=len(ws))
    else:
        wave_range = np.where(np.logical_and(ws <= wb_max, ws >= wb_min))[0]
        return wave_range


class StarSpectrum(GeneralSpectrum):
    """
    Load and operate on spectra

    :param file:
        indicates the name of a file with the spectrum to load.
        It accepts files with the FITS format or plain ASCII
        with one (flux) or two columns.
        Default is None
    :type file: str
    :param wavelength: The wavelength as an array if no file is given
    :type wavelength: numpy.array
    :param flux: The flux as an array if no file is given
    :type flux: numpy.array
    """

    def __init__(self, file=None, wavelength=None, flux=None):
        """

        """
        GeneralSpectrum.__init__(self)
        self.wavelength = wavelength
        self.flux = flux
        self._file = file
        self.__load__()

    def __load__(self):
        if self._file:
            if _file_exists(self._file):
                f_name, f_extension = os.path.splitext(self._file)
                if f_extension == ".fits" or f_extension == ".FITS":
                    self.load_fits(self._file)
                else:
                    self.load_ascii(self._file)
            else:
                warnings.warn("Warning: Could not find file {} - no spectrum loaded" .format(self._file))

    def reflux(self, theta=None):
        """
        Scale the flux by theta**2
        """
        self.flux = self.flux * theta ** 2.

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
        if mag == 'AB':
            f = np.trapz(fr * band(wr) * wr, x=wr) / np.trapz(band(wr) / wr, x=wr) * c
            return -2.5 * np.log10(f) + ab_zero
        elif mag == 'ST':
            f = np.trapz(fr * band(wr) * wr, x=wr) / np.trapz(band(wr) * wr, x=wr)
            return -2.5 * np.log10(f) + st_zero
        elif mag == 'Vega':
            vega_r = _wave_range(band.wavelength, vega_wavelength)
            vega_wr = vega_wavelength[vega_r]
            vega_fr = vega_flux[vega_r]
            f = np.trapz(fr * band(wr) * wr, x=wr)
            vega_f = np.trapz(vega_fr * band(vega_wr) * vega_wr, x=vega_wr)
            return -2.5 * np.log10(f) + 2.5 * np.log10(vega_f) + mag_zero
        else:
            raise ValueError('Magnitude system is not a valid choice, check input string')


class BandPass(GeneralSpectrum):
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
        Initialize a bandpass

        """
        GeneralSpectrum.__init__(self)
        self.wavelength = wavelength
        self.flux = None
        self.name = None
        self.interpolated_response = None
        self._smt = smt
        self.response = response
        self.interpolated_response = None
        
        if band is not None:
            self._load_(band)
        else:
            self.smooth(self._smt)

    def _load_(self, band):
        # Check if input is a band part of our list or an existing file
        if band:
            if band in list_filters:
                self.load_ascii(list_filters[band])
                if self._smt in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
                    self.smooth(kind=self._smt)
                else:
                    raise ValueError('Unknown type of smoothing')

            elif _file_exists(band):
                self.load_ascii(band)
                if self._smt in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
                    self.smooth(kind=self._smt)
                else:
                    raise ValueError('Unknown type of smoothing')

            else:
                warnings.warn("Warning: Could not find filter or file {} ".format(band))
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

    def __call__(self, w):
        """
        Returns the interpolated filter curve function for the wavelength range

        :param w: The wavelengths of the requested area
        :type w: numpy.array
        :return: The values of the interpolated filter curve for these wavelength
        :rtype: numpy.array
        """
        return self.interpolated_response(w)


class BandPassSVO(BandPass):
    def __init__(self, telescope, instrument, filt, smt='linear'):
        """
        A child class of BandPass. It does exactly the same but it will take
        the filter transmission curves from the SVO web-page directly.
    
        :param telescope: The name of the telescope/satellite
        :type telescope: str
        :param instrument: 
            The name of the instrument (if the telescope has only one instrument,
            the name of the instrument is equal to the telescope name)
        :type instrument: str
        :param filt: The name of the band
        :type filt: str
        """
        BandPass.__init__(self, smt=smt)
        filter_curve = svo.get_filter_curve(telescope, instrument, filt)
        self.wavelength = filter_curve['Wavelength']
        self.response = filter_curve['Transmission']
        self.smooth(kind=self._smt)
