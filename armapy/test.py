from armapy.mag_calc import ColorSVO, StarSpectrum, BandSVO
from armapy.svo import SURVEY_SHORTCUT
from astropy.table import Table
from armapy.svo import get_filter_information, get_survey_filter_information
from armapy.plot.plotting import BandPlot


def test_color_calc():
    mass = ColorSVO('2MASS', '2MASS', 'H', 'Ks')

    tab = Table.read('/Users/patrickr/Dropbox/Compact_system/iso_10.11_0.07/spec.fits')
    hcsc = StarSpectrum(wavelength=tab['wavelength'], flux=tab['flux'])
    print('H-Ks=', mass(hcsc, mag_system='Vega'))
    mass = ColorSVO('2MASS', '2MASS', 'J', 'H')

    # tab = Table.read('/Users/patrickr/Dropbox/Compact_system/iso_9.85_0.002/spec.fits')
    hcsc = StarSpectrum(wavelength=tab['wavelength'], flux=tab['flux'])
    print('J-H=', mass(hcsc, mag_system='Vega'))


def test_filter_information():
    try:
        p = get_filter_information('2MASS', '2MASS', 'J')
        print(p)
    except ValueError:
        print('fail')


def test_plotting():
    try:
        tab = Table.read('/Users/patrickr/Dropbox/Compact_system/iso_10.11_0.07/spec.fits')
        hcsc = StarSpectrum(wavelength=tab['wavelength'], flux=tab['flux'])
        band = BandSVO('2MASS', '2MASS', 'J')
        band2 = BandSVO('Paranal', 'VISTA', 'J')
        band3 = BandSVO('2MASS', '2MASS', 'H')
        band4 = BandSVO('Paranal', 'VISTA', 'H')
        band5 = BandSVO('2MASS', '2MASS', 'Ks')
        band6 = BandSVO('Paranal', 'VISTA', 'Ks')
        bp = BandPlot([band, band2,
                       band3, band4,
                       band5, band6], spectra=hcsc)
        bp.show(red_spec=True)
    except ValueError:
        print('fail')


def test_shortcuts():
    print(get_survey_filter_information('gaia', 'BP'))


if __name__ == '__main__':
    test_shortcuts()
