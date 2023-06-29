import dt4dds

# set up config
dt4dds.default_logging()
dt4dds.config.enable_multiprocessing = False
dt4dds.config.n_processes = 1
dt4dds.config.show_progressbars = True


# define primer sequences for PCR
primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]


# 
# Synthesis
# 

# read the sequences from the provided example file
seq_list = dt4dds.tools.txt_to_seqlist("design_sequences.txt")

# create a synthesis instance and process the list of design sequences
array_synthesis = dt4dds.processes.ArraySynthesis()
array_synthesis.process(seq_list)

# sample 100k oligos
pool = array_synthesis.sample_by_counts(100000)
pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
pool.volume = 5


# 
# PCR
# 

# specify primers and total number of cycles, then process the effect of PCR
pcr = dt4dds.processes.PCR(primers=primers_2, n_cycles=30)
pool = pcr.process(pool)


# 
# Sequencing-by-synthesis
# 

# specify current directory as output directory to save the sequencing data
sbs_sequencing = dt4dds.processes.SBSSequencing(output_directory=".")
sbs_sequencing.process(pool)