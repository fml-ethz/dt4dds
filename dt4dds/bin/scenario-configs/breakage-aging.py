import dt4dds
import argparse
import sys
import pathlib

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

seq_primers = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"]


def run(design_file, target_folder, coverage, seq_depth, n_halflives):
    logger.info(f"Running decay challenge on design file {design_file.resolve()} with output folder {target_folder.resolve()}")
    
    # assign efficiency properties
    dt4dds.properties.set_property_settings(
        dt4dds.settings.defaults.SequenceProperties(
            efficiency_distribution='normal',
            efficiency_params={'loc': 1.0, 'scale': 0.0051},
        )
    )

    # read sequences from reference
    seq_list = dt4dds.tools.txt_to_seqlist(design_file)
    n_sequences = len(seq_list)
    logger.info(f"Total number of sequences read: {n_sequences}")

    # use synthesis without base errors
    array_synthesis = dt4dds.processes.ArraySynthesis(dt4dds.settings.defaults.ArraySynthesis_Twist())
    array_synthesis.process(seq_list)

    # perform synthesis for desired physical redundancy
    pool = array_synthesis.sample_by_counts(coverage*n_sequences)
    pool.is_doublestranded = True
    pool.volume = 5

    # perform the decay aging
    aging = dt4dds.processes.Aging(dt4dds.settings.defaults.Aging(
        breakage_enabled=True,
        breakage_ref_length=len(seq_list[0]),
    ), n_halflives=n_halflives)
    pool = aging.process(pool)

    # post-aging adapter ligations
    prep = dt4dds.processes.ssAdapterLigation(
        final_volume=5,
        repair_3prime=1.0,
    )
    pool = prep.process(pool)

    # add indexed primers by PCR
    pcr = dt4dds.processes.PCR(dt4dds.settings.defaults.PCR_Taq(
        primers=seq_primers, 
        n_cycles=20,
    ))
    pool = pcr.process(pool)

    # dilution to 50 pM * 20 uL
    pool = pool.sample_by_moles(0.05*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(dt4dds.settings.defaults.iSeq100(
        output_directory=f"{target_folder}/",
        n_reads=seq_depth*n_sequences,
        primers_read=['ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC']
    ))
    sbs_sequencing.process(pool)




def parse(args):
    parser = argparse.ArgumentParser(description='Scenario: Breakage-Aging')

    parser.add_argument('design_file', type=str, help='Design file to simulate')
    parser.add_argument('target_folder', type=str, help='Target folder for output')
    parser.add_argument('-c', '--coverage', type=float, default=20, help='Target coverage for the pool during aging relative to the number of sequences in the design file')
    parser.add_argument('-d', '--seqdepth', type=float, default=50, help='Sequencing coverage relative to the number of sequences in the design file')
    parser.add_argument('-t', '--halflives', type=float, default=5, help='Number of halflives to simulate during aging')
    args = parser.parse_args(args)

    design_file = pathlib.Path(args.design_file)
    if not design_file.exists():
        raise FileNotFoundError(f"Design file at {args.design_file} does not exist.")
    
    target_folder = pathlib.Path(args.target_folder)
    if not target_folder.exists():
        target_folder.mkdir(parents=True)
        logger.warning(f"Target folder at {args.target_folder} did not exist and was created.")

    return run(design_file, target_folder, args.coverage, args.seqdepth, args.halflives)


if __name__ == "__main__":
    dt4dds.helpers.logging.default_logging()
    parse(sys.argv[1:])