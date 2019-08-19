"""Definitions for the `Kasen0` class."""
import numpy as np
from mosfit.modules.seds.sed import SED
import pickle 
from astropy import constants as c
from tensorflow import keras
import os
import sys


# Important: Only define one ``Module`` class per file.


class Kasen1c(SED):
    """
    Defining the Kasen-simulation based SED
    for the single component model
    Kamile Lukosiute August 2019
    """

    C_CONST = c.c.cgs.value

    def __init__(self, **kwargs):
        super(Kasen1c, self).__init__(**kwargs)
        self.dir = os.path.dirname(os.path.realpath(__file__))
        # self.kasen_wavs = pickle.load(open(os.path.join(self.dir, 'wavelength_angstroms.p'), 'rb'))
        self.kasen_wavs = np.load(os.path.join(self.dir, 'wavelengths.npy'))
        self.kasen_times = pickle.load(open(os.path.join(self.dir, 'times_days.p'), 'rb'))
        self.KSNN = keras.models.load_model(os.path.join(self.dir, 'weights-5805.hdf5'))

        self.closest_wavs_ind = np.vectorize(self.find_closest)

    def find_closest(self, x):
        return np.abs(self.kasen_wavs - x).argmin()

    @staticmethod
    def m_realtonn(x):
        x_log = np.log(x)
        maxi = -2.3025850929940455
        mini = -6.907755278982137
        return (x_log - mini) / (maxi - mini)

    @staticmethod
    def v_realtonn(x):
        maxi = 0.4
        mini = 0.03
        return (x - mini) / (maxi - mini)

    @staticmethod
    def x_realtonn(x):
        x_log = np.log(x)
        maxi = -2.302585092994047
        mini = -20.72326583694641
        return (x_log - mini) / (maxi - mini)

    @staticmethod
    def t_realtonn(x):
        maxi = 14.34999942779541
        mini = 0.05000000074505806
        return (x - mini) / (maxi - mini)

    @staticmethod
    def m_nntoreal(x):
        maxi = -2.3025850929940455
        mini = -6.907755278982137
        return np.exp(x * (maxi - mini) + mini)

    @staticmethod
    def v_nntoreal(x):
        maxi = 0.4
        mini = 0.03
        return x * (maxi - mini) + mini

    @staticmethod
    def x_nntoreal(x):
        maxi = -2.302585092994047
        mini = -20.72326583694641
        return np.exp(x * (maxi - mini) + mini)

    @staticmethod
    def t_nntoreal(x):
        maxi = 14.34999942779541
        mini = 0.05000000074505806
        return x * (maxi - mini) + mini

    #    def total_x_nntoreal(x):
    #      return np.array([self.m_nntoreal(x[0]), self.v_NNtoreal(x[1]), self.x_NNtoreal(x[2]), self.t_NNtoreal(x[3])])

    @staticmethod
    def lum_nntoreal(x):
        # given neural net result of prediction, returns the physical prediction
        mini = -11.512925464970229
        maxi = -0.028132106114705117
        return (np.exp(x * (maxi - mini) + mini) - .00001) * 6.297913324264987e+38

    def process(self, **kwargs):
        # 57982.529
        kwargs = self.prepare_input(self.key('luminosities'), **kwargs)
        self._luminosities = kwargs[self.key('luminosities')]
        self._my_times = kwargs[self.key('rest_times')]
        self._band_indices = kwargs['all_band_indices']
        self._frequencies = kwargs['all_frequencies']

        self._vk = kwargs[self.key('v')]
        self._xlan = kwargs[self.key('X')]
        self._mass = kwargs[self.key('M')]

        cc = self.C_CONST
        zp1 = 1.0 + kwargs[self.key('redshift')]
        czp1 = cc / zp1

        lums = []
        rest_wavs_dict = {}

        kasen_params = np.array([self.m_realtonn(self._mass), self.v_realtonn(self._vk), self.x_realtonn(self._xlan)])
        times_trans = self.t_realtonn(np.array([self._my_times])).T
        nn_input = np.hstack((np.vstack([kasen_params]*len(self._my_times)), times_trans))
        predictions = self.lum_nntoreal(self.KSNN.predict(nn_input)).clip(0) # because physics

        seds = []
        for ti, time in enumerate(self._my_times):
            bi = self._band_indices[ti]
            if bi >= 0:
                rest_wavs = rest_wavs_dict.setdefault(
                    bi, self._sample_wavelengths[bi] / zp1)
            else:
                rest_wavs = np.array(  # noqa: F841
                    [czp1 / self._frequencies[ti]])

            closest_wavs = self.closest_wavs_ind(rest_wavs)

            sed = predictions[ti][closest_wavs]
            if time <= 1e-8:
                sed[sed >= 0.] = 1.e30  # should be all values (can't have neg luminosity

            seds.append(sed)

            lum = np.trapz(predictions[ti], x=self.kasen_wavs)
            self._luminosities[ti] = self._luminosities[ti] + lum
            lums.append(lum)

        # This line turns all nans to 0s
        seds[-1][np.isnan(seds[-1])] = 0.0
        seds = self.add_to_existing_seds(seds, **kwargs)
        return {'sample_wavelengths': self._sample_wavelengths, 'seds': seds,  'lum': lums}
