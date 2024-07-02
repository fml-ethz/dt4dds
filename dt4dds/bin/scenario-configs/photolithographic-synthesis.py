import dt4dds
import argparse
import sys
import pathlib

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def run(design_file, target_folder, coverage, seq_depth):

    logger.info(f"Running scenario on design file {design_file.resolve()} with output folder {target_folder.resolve()}")
    
    seq_primers = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"]

    # assign efficiency properties
    dt4dds.properties.set_property_settings(
        dt4dds.settings.defaults.SequenceProperties(
            efficiency_distribution='normal',
            efficiency_params={'loc': 1.0, 'scale': 0.0051},
        )
    )

    # read sequences from reference
    seq_list = dt4dds.tools.txt_to_seqlist(design_file)
    n_seqs = len(seq_list)
    logger.info(f"Total number of sequences: {n_seqs}")

    # synthesis
    array_synthesis = dt4dds.processes.ArraySynthesis(dt4dds.settings.defaults.ArraySynthesis_Photolithographic())
    array_synthesis.process(seq_list)
    pool = array_synthesis.sample_by_counts(coverage*n_seqs)
    pool.volume = 1

    # post-synthesis adapter ligations
    ligation = dt4dds.processes.AdapterLigation(
        adapter_3prime='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        adapter_5prime='CTCTTTCCCTACACGACGCTCTTCCGATCT',
        final_volume=5,
        tail_bases=[],
        tail_length=0,
        repair_3prime=1.0,
    )
    pool = ligation.process(pool)

    # sequencing PCR with indexed primers
    pcr = dt4dds.processes.PCR(dt4dds.settings.defaults.PCR_Taq(
        polymerase_basesubstitutionrate=0.0,
        primers=seq_primers, 
        n_cycles=15,
    ))
    pool = pcr.process(pool)

    # dilution to 50 pM * 20 uL
    pool = pool.sample_by_moles(0.05*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(
        dt4dds.settings.defaults.iSeq100(
            output_directory=f"{target_folder}/",
            n_reads=seq_depth*n_seqs,
            primers_read=['ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'],
            substitution_rates=[0.0, 0.0],
        )
    )
    sbs_sequencing.process(pool)




def parse(args):
    parser = argparse.ArgumentParser(description='Scenario: Photolithographic Synthesis')

    parser.add_argument('design_file', type=str, help='Design file to simulate')
    parser.add_argument('target_folder', type=str, help='Target folder for output')
    parser.add_argument('-c', '--coverage', type=float, default=100, help='Physical redundancy after synthesis relative to the number of sequences in the design file')
    parser.add_argument('-d', '--seqdepth', type=float, default=30, help='Sequencing coverage relative to the number of sequences in the design file')
  
    args = parser.parse_args(args)
    
    design_file = pathlib.Path(args.design_file)
    if not design_file.exists():
        raise FileNotFoundError(f"Design file at {args.design_file} does not exist.")
    
    target_folder = pathlib.Path(args.target_folder)
    if not target_folder.exists():
        target_folder.mkdir(parents=True)
        logger.warning(f"Target folder at {args.target_folder} did not exist and was created.")

    return run(design_file, target_folder, args.coverage, args.seqdepth)


if __name__ == "__main__":
    dt4dds.helpers.logging.default_logging()
    parse(sys.argv[1:])