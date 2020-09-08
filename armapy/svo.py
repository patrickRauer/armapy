#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 11:00:33 2017

Script to download the list of available filter-curves from the SVO 
(http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?&mode=browse)
website and to download the filter-curves itself.

@author: Patrick Rauer
"""
import urllib
from astropy.table import Table, vstack
from bs4 import BeautifulSoup
import pandas as pd
import os
import configparser
import pkg_resources

SVO_URL = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/{}'

resource_package = 'armapy'
resource_path = '/'.join(('local', 'shortcuts.ini'))
shortcut_path = pkg_resources.resource_filename(resource_package, resource_path)

cache_path = '/'.join(('local', 'cache', 'cache.ini'))
CACHE = configparser.ConfigParser()
try:
    CACHE.read(pkg_resources.resource_filename(resource_package, cache_path))
except FileNotFoundError:
    pass


SURVEY_SHORTCUT = configparser.ConfigParser()
SURVEY_SHORTCUT.read(shortcut_path)


def _split_into_properties(data):
    telescope, instrument, band = [], [], []
    for s in data['Filter ID'].values:
        s = s.split('/')
        telescope.append(s[0])
        s = s[-1].split('.')
        instrument.append(s[0])
        band.append(s[-1])

    data['telescope'] = telescope
    data['instrument'] = instrument
    data['band'] = band


def _read_filters(page_content):
    so = BeautifulSoup(page_content)
    tab = so.find_all('table')[-3]
    filter_infos = pd.read_html(tab.decode(), header=0)[0]
    return filter_infos


def get_svo_filter_list(path='', update=False):
    """
    Creates a complete list of all filter-curves, which are available on
    the SVO web page.
    
    :returns: 
        A table with the telescope names, instrument names and the filter names
    :rtype: astropy.table.Table
    """

    # if a path is given and the list should be updated
    if path != '' and not update:
        if os.path.exists(path):
            return pd.read_csv(path)

    data = []
    for link in links:
        url = link['href']
        if 'mode=browse&' in url:
            with urllib.request.urlopen(SVO_URL.format(url)) as response:
                filter_infos = _read_filters(response.read())
                data.append(filter_infos)
    data = pd.concat(data)

    _split_into_properties(data)

    data.to_csv(path)
    return data


def save_list(surveys, path='./filter_list_svo.fits',
              update=False):
    """
    Saves the list of the available SVO filters

    :param surveys: The output of :meth:`get_svo_filter_list`
    :type surveys: list
    :param path: The path to the save place
    :type path: str
    :param update: True if the current file should be overwrite, else false.
    :type update: bool
    :return: A table with all the filters.
    :rtype: astropy.table.Table
    """
    out = []
    for f in surveys:
        for i in f:
            try:
                tab = Table()
                tab['band'] = i['bands']
                tab['telescope'] = i['survey']
                tab['instrument'] = i['name']
                tab = tab[['telescope', 'instrument', 'band']]
                out.append(tab)
            except TypeError:
                pass

    # stack all the data in one table
    out = vstack(out)
    if path != '' and update:
        out.write(path, overwrite=True)
    return out


def download_filter_curve(telescope, instrument, band, 
                          path, file_type='VO'):
    """
    Downloads a filter transmission curve from the SVO web page.
    
    :param telescope: The name of the telescope/ satellite
    :type telescope: str
    :param instrument: 
        The name of the instrument (if the telescope has only one instrument,
        the name of the instrument is equal to the telescope name)
    :type instrument: str
    :param band: The name of the band
    :type band: str
    :param path: The path to the save place
    :type path: str
    :param file_type: 
        The required file type, available formats are 'ascii' or 'VO'
    :type file_type: str
    """
    # basic url of the SVO web-page
    url = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/'
    # if the file type is 'ascii'
    if file_type == 'ascii':
        url += 'getdata.php?format=ascii&id='
    # if the file type is not 'ascii', which means that will be a VO-table
    else:
        url += 'fps.php?ID='
    # add the telescope, instrument and filter names to the url
    url += '{}/{}.{}'.format(telescope, instrument, band)
    urllib.request.urlretrieve(url, path)
    

def get_filter_curve(telescope, instrument, band, directory=''):
    """
    Returns the specific filter curve back.
    To do that it will download the filter curve as a VO-table and read
    the results back in.
    
    :param telescope: The name of the telescope/satellite
    :type telescope: str
    :param instrument: 
        The name of the instrument (if the telescope has only one instrument,
        the name of the instrument is equal to the telescope name)
    :type instrument: str
    :param band: The name of the band
    :type band: str
    :param directory: The ground directory to store the downloaded files
    :type directory: str
    :returns:
        A table with two columns. 'Wavelength' contains the wavelength in angstrom and 'Transmission' contains
        the transmission at the wavelength. The header of the table has additional information about the filter,
        like the effective wavelength or the zero points in the different filter systems.
    :rtype: astropy.table.Table
    """
    save_path = '{}/{}/{}.fits'
    if directory != '':
        # if the last character is not '/', which means that is the name of the directory and not the path to the
        # directory
        if directory[-1] != '/':
            directory += '/'
        # check if the directory exists
        if not os.path.exists(directory+'{}/{}/'.format(telescope,
                                                        instrument)):
            os.makedirs(directory+'{}/{}/'.format(telescope,
                                                  instrument))
        directory += save_path.format(telescope, 
                                      instrument, 
                                      band)
        # if a file with the name exists, return the file as a astropy.table.Table
        if os.path.exists(directory):
            return Table.read(directory)

    # create the temporary name of the file
    path = './{}-{}-{}.xml'.format(telescope, instrument, band)

    # download the file
    download_filter_curve(telescope, instrument, band,
                          path)

    # read the file in and add the meta information
    tab = Table.read(path)
    tab.meta.update(get_filter_information(telescope, instrument, band))

    # remove the temporary file
    os.remove(path)

    # if a directory is given, save the file
    if directory != '':
        tab.write(directory, overwrite=True)
    return tab


def _get_wavelength_properties(soup_filter):
    wavelengths = pd.read_html(soup_filter.find_all('table')[-6].decode(), header=0)[0]
    wavelengths['Property'] = [s.split(' ')[0].replace('λ', 'lambda_') for s in wavelengths['Property'].values]
    wavelengths = wavelengths.set_index('Property')['Calculated'].to_dict()
    return wavelengths


def _get_zero_points(soup_filter):
    zero_points = []
    for i, s in enumerate(['Vega', 'AB', 'ST']):
        ab = pd.read_html(soup_filter.find_all('table')[-5 + i].decode(), header=0)[0][:2].T[2:].T
        ab['system'] = [f'{s}_{u}' for u in ['ergs', 'jy']]
        zero_points.append(ab)
    zero_points = pd.concat(zero_points)
    return zero_points.set_index('system')['Calculated'].to_dict()


def get_filter_information(telescope, instrument, band):
    """
    Collects the additional filter information from a SVO filter web-page

    :param telescope: The name of the telescope
    :type telescope: str
    :param instrument: The name of the instrument
    :type instrument: str
    :param band: The name of the filter
    :type band: str
    :return: A dict with all the additional information from the web-page
    :rtype: dict
    """
    with urllib.request.urlopen(SVO_URL.format(f'index.php?id={telescope}/{instrument}.{band}')) as response:
        soup_filter = BeautifulSoup(response.read())
    try:
        properties = _get_wavelength_properties(soup_filter)
        properties.update(_get_zero_points(soup_filter))

        for old, new in zip(['Weff', 'FWHM', 'Af/AV'], ['w_eff', 'fwhm', 'af_av']):
            properties[new] = properties.pop(old)

        del properties['lambda_ref']
        del properties['Fsun']
        return properties
    except KeyError:
        raise ValueError(f'The combination of {telescope}, {instrument} and {band} is not valid. '
                         f'Check for misspelling '
                         f'(small and capital letter are important) or check the filter list for the required filter.')


def add_shortcut(shortcut_name, telescope, instrument, overwrite=False):
    """
    Adds a new shortcut to the shortcut configuration file, if the shortcut name didn't exists before.
    If the shortcut exists and overwrite is True, it will overwrite the old entry. Otherwise it
    will do nothing.

    :param shortcut_name: The new shortcut value
    :type shortcut_name: str
    :param telescope: The name of the telescope
    :type telescope: str
    :param instrument: The name of the instrument
    :type instrument: str
    :param overwrite: True if the old version should be overwritten, else False. Default is False.
    :type overwrite: bool
    :return:
    """
    if not SURVEY_SHORTCUT.has_section(shortcut_name):
        SURVEY_SHORTCUT.add_section(shortcut_name)
        SURVEY_SHORTCUT[shortcut_name] = {'telescope': telescope,
                                          'instrument': instrument}
    elif overwrite:
        SURVEY_SHORTCUT[shortcut_name] = {'telescope': telescope,
                                          'instrument': instrument}

    wr = open(shortcut_path, 'w')
    SURVEY_SHORTCUT.write(wr)


def get_survey_filter_information(survey, band):
    """
    Short cut for common large surveys.
    If the survey isn't include in the shortcuts, a KeyError will raise.

    To add a shortcut, use :meth:`add_shortcut`

    :param survey: The name of the survey
    :type survey: str
    :param band: The name of the filter band
    :type band: str
    :return: The information of the filter of the survey
    :rtype: dict
    """
    survey = survey.lower()
    if survey in SURVEY_SHORTCUT.keys():
        s = SURVEY_SHORTCUT[survey]

        # if the filter name is in the config-description
        # replace the band name
        if band in s.keys():
            band = s[band]
        cache_name = '{}_{}'.format(survey, band)
        if CACHE.has_section(cache_name):
            return CACHE[cache_name]
        else:
            properties = get_filter_information(s['telescope'], s['instrument'], band)
            if 'vega' in s.keys():
                properties['vega'] = bool(s['vega'])
            else:
                properties['vega'] = False
            CACHE.add_section(cache_name)
            CACHE[cache_name] = properties
            if not os.path.exists(cache_path):
                p = os.path.split(os.path.abspath(cache_path))[0]
                if not os.path.exists(p):
                    os.makedirs(p)
            wr = open(cache_path, 'w')
            CACHE.write(wr)
            return properties
    else:
        raise KeyError('For {} is no short cut available.'.format(survey))


def get_filter_curve_survey(survey, band):

    """
    Short cut for common large surveys.
    If the survey isn't include in the shortcuts, a KeyError will raise.

    To add a shortcut, use :meth:`add_shortcut`

    :param survey: The name of the survey
    :type survey: str
    :param band: The name of the filter band
    :type band: str
    :return: The information of the filter of the survey
    :rtype: dict
    """
    survey = survey.lower()
    if survey in SURVEY_SHORTCUT.keys():
        s = SURVEY_SHORTCUT[survey]

        # if the filter name is in the config-description
        # replace the band name
        if band in s.keys():
            band = s[band]
        return get_filter_curve(s['telescope'], s['instrument'], band)
    else:
        raise KeyError('For {} is no short cut available.'.format(survey))
