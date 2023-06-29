import dt4dds
import sys

TARGET_FOLDER = sys.argv[1]

dt4dds.config.show_progressbars = False
dt4dds.default_logging()

# very skewed synthesis distribution and high errors by electrochemical synthesis
synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_CustomArray(
    oligo_distribution_type='lognormal',
    oligo_distribution_params={'mean': 0, 'sigma': 1.30},
)

# normal Taq polymerase
polymerase_fidelity = 1


from experiments import run_experiments
run_experiments(synthesis_settings, polymerase_fidelity, TARGET_FOLDER)