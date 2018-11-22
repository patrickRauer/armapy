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
import os
import configparser
import pkg_resources


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
# SURVEY_SHORTCUT = {
#     # '2mass': {'telescope': '2MASS', 'instrument': '2MASS'},
#     # 'sdss': {'telescope': 'SLOAN', 'instrument': 'SDSS'},
#     # 'panstarrs': {'telescope': 'PAN-STARRS', 'instrument': 'PS1'},
#     # 'pan-starrs': {'telescope': 'PAN-STARRS', 'instrument': 'PS1'},
#     # 'ps': {'telescope': 'PAN-STARRS', 'instrument': 'PS1'},
#     # 'kids': {'telescope': 'PARANAL', 'instrument': 'OmegaCam'},
#     # 'viking': {'telescope': 'PARANAL', 'instrument': 'VISTA'},
#     # 'wise': {'telescope': 'WISE', 'instrument': 'WISE'},
#     # 'galex': {'telescope': 'GALEX', 'instrument': 'GALEX'},
#     # 'gaia': {'telescope': 'GAIA', 'instrument': 'GAIA2r'},
#     'ukidss': {'telescope': 'UKRIT', 'instrument': 'UKIDSS'}}


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
            return Table.read(path)
    
    user_agent = ''.join(['Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10_6_4; en-US) ',
                          'AppleWebKit/534.3 (KHTML, like Gecko) Chrome/6.0.472.63 Safari/534.3'])
    headers = {'User-Agent': user_agent}
    
    # download the main web page with all the telescopes
    url = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?&mode=browse'
    req = urllib.request.Request(url, None, headers)
    response = urllib.request.urlopen(req)
    page = response.read()
    response.close()
    page = str(page)
    url_survey = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname={}'

    # split the html-code in parts for the different surveys
    page = page.split('<a href')[2:-3]
    names = []
    surveys = []
    for p in page:
        p = p.split('</a>')[0]
        p = p.split('>')[-1]
        
        # add the survey name to the list
        names.append(p)
        
        # download the telescope web page
        req = urllib.request.Request(url_survey.format(p), None, headers)
        response = urllib.request.urlopen(req)
        survey_page = response.read()
        response.close()
        survey_page = str(survey_page)
        
        # split the telescope web page into the names of the instruments
        filters = survey_page.split('{} filters'.format(p))[-1]
        filters = filters.split('</table')[0]
        filters = filters.split('<a href=\"')[1:]
        filters_links = []
        for f in filters:
            f = f.split('\">')
            filters_links.append({'name': f[1].split('</a')[0],
                                  'survey': p,
                                  'bands': []})
        if len(filters_links) == 0:
            filters_links.append({'name': p,
                                  'survey': p,
                                  'bands': []})
        surveys.append(filters_links)
    instrument_url = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname={}&gname2={}'
        
    # download the html-code of all instruments to get the filter names
    for f in surveys:
        for i in f:
            # download the web-page of an instrument
            req = urllib.request.Request(instrument_url.format(i['survey'], i['name']), None, headers)
            response = urllib.request.urlopen(req)
            instrument_page = response.read()
            response.close()
            instrument_page = str(instrument_page)
            # split the part with the instrument names from the rest of the HTML-code
            instrument_page = instrument_page.split('Description</a>')[-1]
            instrument_page = instrument_page.split('<p><b>Filter Plots')[0]
            instrument_page = instrument_page.split('#filter\">')
            instrument_page = instrument_page[1:]

            # extract the filter names from the code
            for ins in instrument_page:
                ins = ins.split('</a>')
                band = ins[0].split('.')[-1]
                i['bands'].append(band)
                
    surveys = save_list(surveys, path=path, update=update)
    return surveys


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


def __lambda_values__(string):
    """
    Extracts the wavelength information from the string.

    :param string: A part of the HTML-code from a SVO filter web-page
    :type string: str
    :return: The float value of the string
    :rtype: float
    """

    lambda_value = string.split('<td')[2]
    lambda_value = lambda_value.split('</td>')[0]
    lambda_value = lambda_value.split('>')[-1]
    lambda_value = float(lambda_value)
    return lambda_value


def __math_properties__(page):
    """
    Extracts the mathematical properties from a SVO filter web-page

    :param page: The HTML-code of a SVO filter web-page
    :type page: str
    :return: A dict with all the mathematical properties of the filter
    :rtype: dict
    """

    # split the mathematical properties part from the rest of the HTML-code
    math_props = page.split('Mathematical properties')[-1]
    math_props = math_props.split('</table')[0]
    math_props = math_props.split('<tr>')[1:]

    # extract the mathematical properties from the HTML-code
    lambda_mean = __lambda_values__(math_props[0])
    lambda_cen = __lambda_values__(math_props[1])
    lambda_eff = __lambda_values__(math_props[2])
    lambda_peak = __lambda_values__(math_props[3])
    lambda_pivot = __lambda_values__(math_props[4])
    lambda_phot = __lambda_values__(math_props[5])
    lambda_min = __lambda_values__(math_props[6])
    lambda_max = __lambda_values__(math_props[7])
    w_eff = __lambda_values__(math_props[8])
    fwhm = __lambda_values__(math_props[9])
    af_av = __lambda_values__(math_props[10])

    # create the output-dict
    out = {'lambda_mean': lambda_mean,
           'lambda_cen': lambda_cen,
           'lambda_eff': lambda_eff,
           'lambda_peak': lambda_peak,
           'lambda_pivot': lambda_pivot,
           'lambda_phot': lambda_phot,
           'lambda_min': lambda_min,
           'lambda_max': lambda_max,
           'w_eff': w_eff,
           'fwhm': fwhm,
           'af_av': af_av}
    return out


def __zero_points__(calib):
    """
    Extracts the zero points from the HTML-code of a SVO filter web-page

    :param calib: The HTML-code of the zero-point part of a SVO filter web-page
    :type calib: str
    :return: The zero points in units of erg/cm2/s/A and in units of Jy
    :rtype: float, float
    """
    # split the unused part
    calib = calib.split('</table')[0]
    calib = calib.split('Zero Point')[-1]
    calib = calib.split('ZP Type')[0]
    calib = calib.split('<td nowrap>')[-8:-1]

    try:
        # take the erg/cm2/s/A value
        ergs = calib[1].split('</td')[0]
        ergs = float(ergs)
    except ValueError:
        # take the erg/cm2/s/A value
        ergs = calib[0].split('</td')[0]
        ergs = float(ergs)
    try:
        # take the Jy value
        jy = calib[-2].split('</td')[0]
        jy = float(jy)
    except ValueError:
        # take the Jy value
        jy = calib[-3].split('</td')[0]
        jy = float(jy)

    return ergs, jy


def __all_zero_points__(page):
    """
    Extracts all zero points of the different systems from the HTML-code of SVO filter web-page

    :param page: The HTML-code of a SVO filter web-page
    :type page: str
    :return: A dict with the zero points of the different magnitude systems
    :rtype: dict
    """
    vega_ergs, vega_jy = __zero_points__(page.split('Vega System')[-1])
    ab_ergs, ab_jy = __zero_points__(page.split('AB System')[-1])
    st_ergs, st_jy = __zero_points__(page.split('ST System')[-1])

    # prepare the output dict
    out = {'Vega_ergs': vega_ergs,
           'Vega_jy': vega_jy,
           'AB_ergs': ab_ergs,
           'AB_jy': ab_jy,
           'ST_ergs': st_ergs,
           'ST_jy': st_jy}
    return out
    

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
    # create the URL to the SVO filter web-page
    url = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?'
    url += 'id={}/{}.{}&&mode=browse&gname={}&gname2={}#filter'.format(telescope,
                                                                       instrument,
                                                                       band,
                                                                       telescope,
                                                                       instrument)

    # download the HTML-code of the web-page
    user_agent = ''.join(['Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10_6_4; en-US) ',
                          'AppleWebKit/534.3 (KHTML, like Gecko) Chrome/6.0.472.63 Safari/534.3'])
    headers = {'User-Agent': user_agent}
    req = urllib.request.Request(url, None, headers)
    response = urllib.request.urlopen(req)
    page = response.read()
    response.close()
    page = str(page)

    # extract the mathematical properties
    math_props = __math_properties__(page)

    # extract the zero points
    zero_points = __all_zero_points__(page)
    
    math_props.update(zero_points)
    return math_props


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
