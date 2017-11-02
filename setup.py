from distutils.core import setup

setup(
    name='armapy',
    version='0.7',
    packages=['armapy'],
    url='',
    license='GPL',
    author='Patrick Rauer',
    author_email='j.p.rauer@sron.nl',
    description='Package to calculate magnitudes from a spectra',
    requires=['astropy', 'numpy', 'scipy']
)
