import numpy as np

from . import tools

import logging
logger = logging.getLogger(__name__)


show_progressbars = False
""" Whether to show progress bars for time-consuming operations in the console. """

error_generation_max_coverage_threshold = 10000
""" Maximum threshold coverage (oligos/sequence) to which simulation is limited and result is instead scaled-up. """

error_generation_error_coverage = 5
""" Average number of mutated oligos of one sequence to simulate for any error type. Scales with error rate to determine coverage threshold. """

seqpool_maximum_pool_diversity = 1e7
""" Number of individual oligos any pool is limited to, before a subset of this size is scaled up to represent it instead. """

multiprocessing_sequences_batch_size = 10000
""" Maximum batch size in number of sequences for one multiprocessing batch. """


def set_random_seed(seed):
    """  """
    tools.rng = np.random.default_rng(np.random.SFC64(seed))