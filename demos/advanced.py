import dt4dds

# set up config
dt4dds.default_logging()
dt4dds.config.enable_multiprocessing = False
dt4dds.config.n_processes = 1
dt4dds.config.show_progressbars = True


# define primer sequences for PCR
primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]

# assign efficiency properties for amplification, here a normal distribution with std of 0.0051
dt4dds.properties.set_property_settings(
    dt4dds.settings.defaults.SequenceProperties(
        efficiency_distribution='normal',
        efficiency_params={'loc': 1.0, 'scale': 0.0051},
    )
)


# 
# Electrochemical synthesis with specified coverage bias
# 

# read the sequences from the provided example file
seq_list = dt4dds.tools.txt_to_seqlist("design_sequences.txt")
n_seqs = len(seq_list)

# set up the synthesis by using defaults for electrochemical synthesis
synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Twist()
# settings can be customized further when passing to the process instance
array_synthesis = dt4dds.processes.ArraySynthesis(synthesis_settings(
    oligo_distribution_type='lognormal',
    oligo_distribution_params={'mean': 0, 'sigma': 0.30},
))
array_synthesis.process(seq_list)

# sample with mean coverage of 200
pool = array_synthesis.sample_by_counts(200*n_seqs)
pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
pool.volume = 1


# 
# Aging for one half-live
# 
aging_settings = dt4dds.settings.defaults.Aging()
aging = dt4dds.processes.Aging(aging_settings(
    fixed_decay_ratio=0.5,
))
pool = aging.process(pool)
pool.volume = 1


# 
# PCR with High Fidelity polymerase for 30 cycles at mean efficiency of 95%
# 
pcr_settings = dt4dds.settings.defaults.PCR_HiFi()
pcr = dt4dds.processes.PCR(pcr_settings(
    primers=primers_2,
    template_volume=1,
    volume=20,
    efficiency_mean=0.95,
    n_cycles=30,
))
pool = pcr.process(pool)


# 
# Sequencing-by-synthesis with paired reads and sequencing coverage of 25
# 
synthesis_settings =  dt4dds.settings.defaults.iSeq100()
sbs_sequencing = dt4dds.processes.SBSSequencing(synthesis_settings(
    output_directory=".",
    n_reads=int(25*n_seqs),
    read_length=150,
    read_mode='paired-end',
))
sbs_sequencing.process(pool)