import dt4dds
import shutil
import gc
import sys

TARGET_FOLDER = sys.argv[1]

dt4dds.config.show_progressbars = False
dt4dds.config.enable_multiprocessing = False
dt4dds.config.n_processes = 4
dt4dds.default_logging()


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


# common settings for PCR
pcr_settings = dt4dds.settings.defaults.PCR(
    primers=primers_0,
    template_volume=1,
    volume=20,
    efficiency_mean=0.94,
    polymerase_fidelity=1,
    polymerase_basesubstitutionrate=1.97e-4,
)


# common sequencing workflow
def sequence(seqpool, name):

    # dilution to 50 pM * 20 uL
    pool = seqpool.sample_by_moles(0.05*20e-6)

    # sequencing
    sbs_sequencing = dt4dds.processes.SBSSequencing(
        dt4dds.settings.defaults.iSeq100(
            output_directory=f"{TARGET_FOLDER}/{name}/",
            n_reads=1000000,
        )
    )
    sbs_sequencing.process(pool)
    shutil.copyfile('./design_files.fasta', f"{TARGET_FOLDER}/{name}/design_files.fasta")





# read sequences from reference
seq_list = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')

# use synthesis settings as passed
array_synthesis = dt4dds.processes.ArraySynthesis(dt4dds.settings.defaults.ArraySynthesis_CustomArray_GCfix(
    oligo_distribution_type='lognormal',
    oligo_distribution_params={'mean': 0, 'sigma': 0.58},
))
array_synthesis.process(seq_list)

# perform synthesis
pool = array_synthesis.sample_by_mass(0.0302*50)
pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
pool = pool.sample_by_mass(0.0302*50)
pool.volume = 50

#
# PCR0
#
ipool = pool.sample_by_volume(5)
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=10, template_volume=5))
ipool = pcr.process(ipool)
sequence(ipool, f'PCR0')
print("PCR0 done")
del ipool
gc.collect()


#
# PCR1
#
ipool = pool.sample_by_volume(2)
ipool.volume = 5
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=15, template_volume=1))
ipool = pcr.process(ipool)
sequence(ipool, f'PCR1')
print("PCR0 done")
del ipool
gc.collect()


#
# PCR2
#
pool = pool.sample_by_volume(10)
pool.volume = 500
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=20, template_volume=1))
ipool = pcr.process(pool)
sequence(ipool, f'PCR2')
print("PCR0 done")
del ipool
gc.collect()


#
# PCR3
#
pool = pool.sample_by_volume(2)
pool.volume = 40
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=25, template_volume=1))
ipool = pcr.process(pool)
sequence(ipool, f'PCR3')
print("PCR0 done")
del ipool
gc.collect()