from .. import datastructures

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)





class AbstractDilution(Step):


    #
    # Basic functions
    #

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def __repr__(self):
        return f"{type(self).__name__}()"


    def process(self, source_pool: datastructures.SeqPool):
        """ Main entry point. Performs the dilution. """

        # set up
        self._run_pre_process_hooks()

        # perform dilution
        diluted_pool = self._dilute(source_pool)

        # finish up
        self._run_post_process_hooks()

        return diluted_pool


    def _dilute(self, source_pool: datastructures.SeqPool):
        raise NotImplementedError("You must call a subclass of dilution.")





class FixedVolumeDilution(AbstractDilution):


    def __init__(self, sample_volume = 1, final_volume = 50):
        super().__init__()

        self.sample_volume = sample_volume
        self.final_volume = final_volume


    def _dilute(self, source_pool: datastructures.SeqPool):
        """  """
        if source_pool.volume is None:
            raise AttributeError("Provided SeqPool does not have a specified volume, cannot perform dilution.")
        
        sampling_ratio = self.sample_volume/source_pool.volume

        diluted_pool = source_pool.sample_by_factor(sampling_ratio, remove_sampled_oligos=True)
        diluted_pool.volume = self.diluted_volume

        return diluted_pool



class FixedMassConcentrationDilution(AbstractDilution):

    def __init__(self, concentration = 1, final_volume = 50):
        super().__init__()

        self.target_conc = concentration
        self.target_volume = final_volume


    def _dilute(self, source_pool: datastructures.SeqPool):
        """ """
        if source_pool.volume is None:
            raise AttributeError("Provided SeqPool does not have a specified volume, cannot perform dilution.")

        init_conc = source_pool.mass_concentration
        dil_factor = self.target_conc/init_conc

        if dil_factor > 1:
            raise AssertionError(f"Seqpool at concentration {init_conc} cannot be diluted to a higher concentration of {self.target_conc}.")

        sampled_volume = self.target_volume*dil_factor
        sampling_ratio = sampled_volume/source_pool.volume

        diluted_pool = source_pool.sample_by_factor(
            sampling_ratio, 
            remove_sampled_oligos=True
        )
        diluted_pool.volume = self.target_volume

        return diluted_pool



class FixedMolarConcentrationDilution(AbstractDilution):

    def __init__(self, concentration = 1, final_volume = 50):
        super().__init__()

        self.target_conc = concentration
        self.target_volume = final_volume


    def _dilute(self, source_pool: datastructures.SeqPool):
        """ """
        if source_pool.volume is None:
            raise AttributeError("Provided SeqPool does not have a specified volume, cannot perform dilution.")

        init_conc = source_pool.molar_concentration
        dil_factor = self.target_conc/init_conc

        if dil_factor > 1:
            raise AssertionError(f"Seqpool at concentration {init_conc} cannot be diluted to a higher concentration of {self.target_conc}.")

        sampled_volume = self.target_volume*dil_factor
        sampling_ratio = sampled_volume/source_pool.volume

        diluted_pool = source_pool.sample_by_factor(
            sampling_ratio, 
            remove_sampled_oligos=True
        )
        diluted_pool.volume = self.target_volume

        return diluted_pool