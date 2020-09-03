from distutils.core import setup

setup(
    name='armapy',
    version='0.8.3',
    packages=['armapy', 'armapy.plot'],
    url='',
    license='GPL',
    author='Patrick Rauer',
    author_email='j.p.rauer@sron.nl',
    description='Package to calculate magnitudes from a spectra',
    requires=['astropy', 'numpy', 'scipy'],
    package_data={
          '': ['./alpha_lyr_stis_006.fits']
      },
)
