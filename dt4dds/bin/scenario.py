import dt4dds
import argparse
import sys
import importlib.util
import pathlib

dt4dds.config.show_progressbars = False
dt4dds.default_logging()

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

AVAILABLE_SCENARIOS = {}
for scenario in pathlib.Path(__file__).parent.glob("scenario-configs/*.py"):
    scenario_name = scenario.stem
    AVAILABLE_SCENARIOS[scenario_name] = scenario.resolve()


def load_and_run_scenario(scenario_name):
    scenario = AVAILABLE_SCENARIOS[scenario_name]
    spec = importlib.util.spec_from_file_location("scenario", scenario)
    module = importlib.util.module_from_spec(spec)
    sys.modules["scenario"] = module
    spec.loader.exec_module(module)
    module.parse(sys.argv[2:])


def main():
    parser = argparse.ArgumentParser(description='Running pre-defined scenarios')
    parser.add_argument('scenario', type=str, help='Name of the scenario to run', choices=AVAILABLE_SCENARIOS.keys())

    args = parser.parse_args(sys.argv[1:2])

    logger.info(f"Running scenario {args.scenario}.")
    load_and_run_scenario(args.scenario)
