import pylab as pl
import numpy as np


class BandPlot:
    band = None
    spectra = None

    def __init__(self, band, spectra=None):
        self.band = band
        self.spectra = spectra

    def show(self, red_spec=False, nbins=250):
        """
        Plot the bandpass and the spectra, if a spectra was set.

        :param red_spec:
            True if the set spectra should be rebinned (to make the chart clearer) or
            False if the original spectra should be taken.
            It has no effect if no spectra was set.
        :type red_spec: bool
        :param nbins:
            The number of new bins if the spectra should be rebinned.
            It has no effect if red_spec is False and/or no spectra was set.
            Default is 250 bins.
        :type nbins: int
        :return:
        """
        fig = pl.figure()
        fig.clf()
        fig.subplots_adjust(bottom=0.4,
                            top=0.95,
                            left=0.1,
                            right=0.91)
        sp = fig.add_subplot(111)

        # if there is a list of bandpass
        if type(self.band) == list:
            wave_min = 1e9
            wave_max = 0
            for b in self.band:
                # if telescope name and instrument name are the
                # same, take it one time only
                if b.telescope == b.instrument:
                    label = '{} {}'.format(b.instrument,
                                           b.filter_band)
                # if telescope name and instrument name aren't the same
                # take both
                else:
                    label = '{} {} {}'.format(b.telescope,
                                              b.instrument,
                                              b.filter_band)
                sp.plot(b.wavelength,
                        b.response,
                        label=label,
                        zorder=20, lw=2)
                wave_min = min(wave_min,
                               np.min(b.wavelength))
                wave_max = max(wave_max, np.max(b.wavelength))
        else:
            if self.band.telescope == self.band.instrument:
                label = '{} {}'.format(self.band.telescope,
                                       self.band.filter_band)
            else:
                label = '{} {} {}'.format(self.band.telescope,
                                          self.band.instrument,
                                          self.band.filter_band)
            sp.plot(self.band.wavelength,
                    self.band.response,
                    label=label,
                    zorder=20, lw=2)
            wave_min = np.min(self.band.wavelength)
            wave_max = np.max(self.band.wavelength)

        pl.legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)

        # if a spectra is set
        if self.spectra is not None:
            sp2 = sp.twinx()

            # cut the spectra at the edges of the bandpass curve
            mask = (self.spectra.wavelength > wave_min) & (self.spectra.wavelength < wave_max)
            wavelength = self.spectra.wavelength[mask]
            flux = self.spectra.flux[mask]

            # if a reduced spectra is required
            if red_spec:
                waves = np.linspace(wave_min, wave_max, nbins+1)
                fluxes_new = []
                for wave_start, wave_end in zip(waves[:-1], waves[1:]):
                    m = (wavelength >= wave_start) & (wavelength < wave_end)
                    fluxes_new.append(np.mean(flux[m]))

                # shift the wavelength by a half and cut the last item
                wavelength = waves[:-1]+(waves[1]-waves[0])/2
                flux = np.array(fluxes_new)
            sp2.step(wavelength,
                     flux,
                     color='k',
                     zorder=1, alpha=0.5)
            sp2.set_ylabel("flux")

        sp.set_xlabel('wavelength [angstrom]')
        sp.set_ylabel('efficiency')
        pl.show()
