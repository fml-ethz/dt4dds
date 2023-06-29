import numpy as np
import pprint

from .helpers import tools
from . import settings

import logging
logger = logging.getLogger(__name__)



class AmplificationEfficiency():

    _efficiency_cache = {}

    _settings = settings.PropertiesSettings()


    @classmethod
    def _efficiency_distribution(cls, size):
        """ get user-specified efficiency distribution """
        if not hasattr(tools.rng, cls._settings.efficiency_distribution):
            raise NotImplementedError(f"A distribution called {cls._settings.efficiency_distribution} is not implemented.")

        return getattr(tools.rng, cls._settings.efficiency_distribution)(**cls._settings.efficiency_params)


    @classmethod
    def get_efficiencies(cls, sequences):
        """  """
        # get efficiencies for known oligos
        efficiencies = np.array([cls._efficiency_cache.get(seq, np.nan) for seq in sequences])

        # create new efficiencies for unknown oligos
        missing_mask = np.isnan(efficiencies)

        logger.info(f"Unknown efficiencies: {np.count_nonzero(missing_mask)}, {100*np.count_nonzero(missing_mask)/len(sequences):.2f}%")

        # short-circuit if no missing efficiencies
        if not np.count_nonzero(missing_mask):
            return efficiencies

        # generate missing efficiencies
        efficiencies[missing_mask] = cls._efficiency_distribution(np.count_nonzero(missing_mask))

        for ix in missing_mask.nonzero()[0]:
            cls._efficiency_cache[sequences[ix]] = efficiencies[ix]

        return efficiencies


    @classmethod
    def get_efficiency(cls, sequence):
        """  """
        if sequence not in cls._efficiency_cache.keys():
            cls._efficiency_cache[sequence] = cls._efficiency_distribution(1)
        
        return cls._efficiency_cache[sequence]

    
    @classmethod
    def duplicate_efficiencies(cls, from_sequence, to_sequences):
        """  """
        # get efficiency for original sequence
        efficiency = cls.get_efficiency(from_sequence)

        # create identical efficiency for other oligos
        for seq in to_sequences:
            cls._efficiency_cache[seq] = efficiency



def set_property_settings(settings):
    AmplificationEfficiency._settings = settings
    if logger.level <= logging.INFO:
        settings_dict = pprint.pformat(settings, compact=True, width=10000000)
        logger.info(f"Updated settings for properties: {settings_dict}")