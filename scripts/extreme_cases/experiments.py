import dt4dds
import shutil
import gc


primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
primers_1 = ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]
primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]

# assign efficiency properties
dt4dds.properties.set_property_settings(
    dt4dds.settings.defaults.SequenceProperties(
        efficiency_distribution='normal',
        efficiency_params={'loc': 1.0, 'scale': 0.0051},
    )
)

# common sequencing workflow for bunnies
def sequence(seqpool, name, folder, pcr_settings):

    pool = seqpool.copy()
    
    # samples dilution
    pool = pool.sample_by_counts(500*12000)
    pool.volume = 1

    # sequencing PCR with indexed primers
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=10, primers=primers_2))
    pool = pcr.process(pool)

    # dilution to 20 pM * 20 uL
    pool = pool.sample_by_moles(0.02*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(
        dt4dds.settings.defaults.iSeq100(
            output_directory=f"{folder}/{name}/",
            n_reads=2000000,
        )
    )
    sbs_sequencing.process(pool)
    shutil.copyfile('./design_files.fasta', f"{folder}/{name}/design_files.fasta")





def run_experiments(synthesis_settings, fidelity, folder):

    # common settings for PCR
    pcr_settings = dt4dds.settings.defaults.PCR(
        primers=primers_0,
        template_volume=1,
        volume=20,
        efficiency_mean=0.95,
        polymerase_fidelity=fidelity
    )


    # read sequences from reference
    seq_list = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')

    # use synthesis settings as passed
    array_synthesis = dt4dds.processes.ArraySynthesis(synthesis_settings)
    array_synthesis.process(seq_list)

    # perform synthesis
    pool = array_synthesis.sample_by_mass(0.001)
    pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
    pool = pool.sample_by_counts(200*len(seq_list))
    pool.volume = 1

    # amplification PCR
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=20))
    pool = pcr.process(pool)

    # 
    # master pool
    # 
    sequence(pool, f'MasterPool', folder, pcr_settings)
    print("MasterPool done")
    gc.collect()


    # 
    # different coverage experiments
    # 
    for i in [50, 45, 40, 35, 32.5, 30, 27.5, 25, 22.5, 20, 18, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.5]:
        
        # sample exact coverage
        ipool = pool.sample_by_counts(i*len(seq_list))
        ipool.volume = 1

        aging = dt4dds.processes.Aging(
            dt4dds.settings.defaults.Aging(
                fixed_decay_ratio=0.5,
                substitution_rate=0.0,
            )
        )
        ipool = aging.process(ipool)
        ipool.volume = 1

        # amplification PCR
        pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=30))
        ipool = pcr.process(ipool)

        # sequencing
        sequence(ipool, f'Pool_{str(int(10*i)).zfill(3)}', folder, pcr_settings)
        del ipool
        print(f'Pool_{str(int(10*i)).zfill(3)} done')
        gc.collect()
