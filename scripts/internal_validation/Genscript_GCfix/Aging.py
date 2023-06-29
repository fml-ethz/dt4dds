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
    template_volume=5,
    volume=20,
    efficiency_mean=0.94,
    polymerase_fidelity=1,
    polymerase_basesubstitutionrate=1.97e-4,
)


# common sequencing workflow for bunnies
def sequence(seqpool, name):

    # samples dilution
    pool = seqpool.sample_by_mass(5*0.5)
    pool.volume = 5

    # sequencing PCR with indexed primers
    pcr = dt4dds.processes.PCR(pcr_settings(primers=primers_2, n_cycles=7))
    pool = pcr.process(pool)

    # dilution to 50 pM * 20 uL
    pool = pool.sample_by_moles(0.05*20e-6)

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
pool = array_synthesis.sample_by_mass(0.0302*1)
pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
pool = pool.sample_by_mass(0.0302*10/50)
pool.volume = 500

# initial amplification
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=20, template_volume=96*1, volume=96*20))
pool = pcr.process(pool)



# 
# 0a
# 

# initial dilution for aging
ipool = pool.sample_by_mass(30)
ipool.volume = 200
ipool = ipool.sample_by_volume(5)
ipool.volume = 50

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=10))
ipool = pcr.process(ipool)

sequence(ipool, "0a")
print(f'0a done')
del ipool
gc.collect()



# 
# 0b
# 

# initial dilution for aging
ipool = pool.sample_by_mass(30)
ipool.volume = 200
ipool = ipool.sample_by_volume(5)
ipool.volume = 540

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=16))
ipool = pcr.process(ipool)

sequence(ipool, "0b")
print(f'0b done')
del ipool
gc.collect()



# 
# 2d
# 

# initial dilution for aging
ipool = pool.sample_by_mass(30)

# aging
aging = dt4dds.processes.Aging(dt4dds.settings.defaults.Aging(fixed_decay_ratio=1-0.0832))
ipool = aging.process(ipool)

# post-aging dilution
ipool.volume = 200
ipool = ipool.sample_by_volume(5)
ipool.volume = 60.5

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=16))
ipool = pcr.process(ipool)

sequence(ipool, "2d")
print(f'2d done')
del ipool
gc.collect()



# 
# 4d
# 

# initial dilution for aging
ipool = pool.sample_by_mass(30)

# aging
aging = dt4dds.processes.Aging(dt4dds.settings.defaults.Aging(fixed_decay_ratio=1-0.0272))
ipool = aging.process(ipool)

# post-aging dilution
ipool.volume = 200
ipool = ipool.sample_by_volume(15)
ipool.volume = 50.7

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=16))
ipool = pcr.process(ipool)

sequence(ipool, "4d")
print(f'4d done')
del ipool
gc.collect()



# 
# 7d
# 

# initial dilution for aging
ipool = pool.sample_by_mass(30)

# aging
aging = dt4dds.processes.Aging(dt4dds.settings.defaults.Aging(fixed_decay_ratio=1-0.00942))
ipool = aging.process(ipool)

# no post-aging dilution
ipool.volume = 200

# amplification PCR
pcr = dt4dds.processes.PCR(pcr_settings(n_cycles=16))
ipool = pcr.process(ipool)

sequence(ipool, "7d")
print(f'7d done')
del ipool
gc.collect()
