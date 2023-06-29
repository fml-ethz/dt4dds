from .. import datastructures

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)





class AbstractPurification(Step):


    #
    # Basic functions
    #

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def __repr__(self):
        return f"{type(self).__name__}()"


    def process(self, source_pool: datastructures.SeqPool):
        """ Main entry point. Performs the purification. """

        # set up
        self._run_pre_process_hooks()

        # perform purification
        purified_pool = self._purify(source_pool)

        # finish up
        self._run_post_process_hooks()

        return purified_pool


    def _purify(self, source_pool: datastructures.SeqPool):
        raise NotImplementedError("You must call a subclass of purification.")






class StandardPurification(AbstractPurification):

    def __init__(self, recovery = 0.9, final_volume = 50):
        super().__init__()

        if recovery > 1 or recovery < 0:
            raise ArithmeticError(f"A recovery of {recovery} is impossible (0.0 <= recovery <= 1.0).")

        self.recovery = recovery
        self.final_volume = final_volume


    def _purify(self, source_pool: datastructures.SeqPool):
        """  """

        purified_pool = source_pool.sample_by_factor(self.recovery, remove_sampled_oligos=False)
        purified_pool.volume = self.final_volume

        return purified_pool




class SizePurification(AbstractPurification):

    def __init__(self, recovery = 0.9, final_volume = 50, length_window=[100, 102]):
        super().__init__()

        if recovery > 1 or recovery < 0:
            raise ArithmeticError(f"A recovery of {recovery} is impossible (0.0 <= recovery <= 1.0).")
        
        self.recovery = recovery
        self.final_volume = final_volume
        self.length_window = length_window


    def _purify(self, source_pool: datastructures.SeqPool):
        """  """

        size_pool = source_pool.filter_by_length(self.length_window[0], self.length_window[1])

        purified_pool = size_pool.sample_by_factor(self.recovery, remove_sampled_oligos=False)
        purified_pool.volume = self.final_volume

        return purified_pool