import dt4dds
import sys
import shutil
import gc

TARGET_FOLDER = sys.argv[1]

primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
primers_1 = ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]
primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]

dt4dds.config.show_progressbars = False
dt4dds.default_logging()

# assign efficiency properties
dt4dds.properties.set_property_settings(
    dt4dds.settings.defaults.SequenceProperties(
        efficiency_distribution='normal',
        efficiency_params={'loc': 1.0, 'scale': 0.012},
    )
)


# common settings for PCR
pcr_settings = dt4dds.settings.defaults.PCR(
    primers=primers_0,
    template_volume=1,
    volume=20,
    efficiency_mean=0.98,
)


# sequencing workflow for master pool
def sequence_master(seqpool, name):

    pool = seqpool.copy()
    
    # initial dilution
    pool = pool.sample_by_mass(0.01)
    pool.volume = 10

    # first sequencing PCR
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=12, primers=primers_1))
    pool = pcr.process(pool)

    # purification, 10x dilution
    pool.volume *= 10

    # second sequencing PCR
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=10, primers=primers_2))
    pool = pcr.process(pool)

    # dilution to 20 pM * 20 uL
    pool = pool.sample_by_moles(0.02*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(
        dt4dds.settings.defaults.iSeq100(
            output_directory=f"{TARGET_FOLDER}/{name}/",
            n_reads=1000000,
        )
    )
    sbs_sequencing.process(pool)
    shutil.copyfile('./design_files.fasta', f"{TARGET_FOLDER}/{name}/design_files.fasta")


# common sequencing workflow for bunnies
def sequence(seqpool, name):

    pool = seqpool.copy()
    
    # samples were equally diluted
    pool = pool.sample_by_mass(0.005)
    pool.volume = 50

    # first sequencing PCR
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=25, primers=primers_1))
    pool = pcr.process(pool)

    # purification, 10x dilution
    pool.volume *= 10

    # second sequencing PCR
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=10, primers=primers_2))
    pool = pcr.process(pool)

    # dilution to 20 pM * 20 uL
    pool = pool.sample_by_moles(0.02*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(
        dt4dds.settings.defaults.iSeq100(
            output_directory=f"{TARGET_FOLDER}/{name}/",
            n_reads=2000000,
        )
    )
    sbs_sequencing.process(pool)
    shutil.copyfile('./design_files.fasta', f"{TARGET_FOLDER}/{name}/design_files.fasta")



# read sequences from reference
seq_list = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')

# synthesis setup
array_synthesis = dt4dds.processes.ArraySynthesis(
    dt4dds.settings.defaults.ArraySynthesis_CustomArray(
        oligo_distribution_type='lognormal',
        oligo_distribution_params={'mean': 0, 'sigma': 0.94},
    )
)
array_synthesis.process(seq_list)

# perform synthesis, sample for total put into wells
pool = array_synthesis.sample_by_mass(1)
pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)

# dilute 1:46 for all 96 wells
pool.volume = 96*46


# 
# BUNNY_M
# 

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=12))
pool = pcr.process(pool)

# sequencing
sequence_master(pool, f'Bunny_M')
print("Bunny_M done")
gc.collect()


# 
# BUNNY_P + BUNNY_Fx
# 

# dilute master library for amplification
pool = pool.sample_by_mass(0.025)
pool.volume = 50

# generations
for i in ["M", "P", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"]:

    # amplification PCR between generations
    pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=12))
    pool = pcr.process(pool)

    # sample 25 pg DNA from filament in 50 uL
    pool = pool.sample_by_mass(0.025)
    pool.volume = 50

    # sequencing
    sequence(pool, f'Bunny_{i}')
    print(f"Bunny_{i} done")
    gc.collect()
