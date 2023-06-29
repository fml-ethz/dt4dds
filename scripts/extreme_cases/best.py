import dt4dds
import sys

TARGET_FOLDER = sys.argv[1]

dt4dds.config.show_progressbars = False
dt4dds.default_logging()

# narrow synthesis distribution and low errors by Twist
synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Twist(
    oligo_distribution_type='lognormal',
    oligo_distribution_params={'mean': 0, 'sigma': 0.30},
)

# highest fidelity polymerase
polymerase_fidelity = 280


from experiments import run_experiments
run_experiments(synthesis_settings, polymerase_fidelity, TARGET_FOLDER)