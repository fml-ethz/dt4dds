# 🧬 dt4dds - Digital Twin for DNA Data Storage


- [Overview](#overview)
- [Web-based Tool](#web-based-tool)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)


# Overview
`dt4dds` is a Python package  providing a customizable, digital representation of the widely-used DNA data storage workflow involving array synthesis, PCR, Aging, and Sequencing-By-Synthesis. By modelling each part of such user-defined workflows with fully customizable experimental parameters, `dt4dds` enables data-driven experimental design and rational design of redundancy. `dt4dds` also includes a pipeline for comprehensively analyzing errors in sequencing data, both from experiments and simulation. This Python package is used in the following publications:

> Gimpel, A.L., Stark, W.J., Heckel, R., Grass R.N. A digital twin for DNA data storage based on comprehensive quantification of errors and biases. Nat Commun 14, 6026 (2023). https://doi.org/10.1038/s41467-023-41729-1

> Gimpel, A.L., Stark, W.J., Heckel, R., Grass R.N. Challenges for error-correction coding in DNA data storage: photolithographic synthesis and DNA decay. bioRxiv 2024.07.04.602085 (2024). https://doi.org/10.1101/2024.07.04.602085

The Jupyter notebooks and associated code used for generating the figures in the manuscript are found in the [dt4dds_notebooks repository](https://github.com/fml-ethz/dt4dds_notebooks).

## New: Challenges for DNA Data Storage

A C++ implementation of the Digital Twin for DNA Data Storage for two current challenges in error-correction coding for DNA - Photolithographic DNA Synthesis and DNA Decay - is available [in this GitHub repository](https://github.com/fml-ethz/dt4dds-challenges). More information is also provided in the following publication:

> Gimpel, A.L., Stark, W.J., Heckel, R., Grass R.N. Challenges for error-correction coding in DNA data storage: photolithographic synthesis and DNA decay. bioRxiv 2024.07.04.602085 (2024). https://doi.org/10.1101/2024.07.04.602085


# Web-based Tool
A web-based version of `dt4dds` with an easy-to-use graphical user interface and the most common workflows is available at [dt4dds.ethz.ch](https://dt4dds.ethz.ch).


# System requirements
## Hardware requirements
This package only requires a standard computer. Depending on the size and complexity of the simulated workflows, sufficient RAM to support the in-memory operations is required. The required amount of RAM can be reduced by decreasing the number of cores used for parallelization (see config below), at the cost of increased run time.

## Software requirements
This package is compatible with Windows, macOS and Linux. The package has been developed and tested on Ubuntu 20.04 using Python 3.10. The Python packages listed in [requirements.txt](/requirements.txt) are required.


# Installation guide
To install this package from PyPi, use
```bash
pip3 install dt4dds
```
To install this package from Github, use
```bash
git clone https://github.com/fml-ethz/dt4dds
cd dt4dds
python3 setup.py install .
```



# Demo 
Exemplary scripts using this package are provided [in the demo folder](demos/):

- [Basic demo](demos/basic.py) showing a simple data storage workflow
- [Advanced demo](demos/advanced.py) showcasing a workflow with custom experimental parameters 

All scripts run within seconds on a standard computer and require less than 1 GB of RAM.


# Scripts for reproduction
All scripts used in the manuscript for internal and external validation, as well as the case study are provided [in the scripts folder](scripts/):

- [Interal validation](scripts/internal_validation/) reproducing the experimental conditions of the manuscript
- [Generational experiments by Koch et al.](scripts/bunny_generations/) for external validation
- [Case study](scripts/extreme_cases/) as an example application for the rational design of redundancy

To run these scripts, provide a valid filepath where output files will be stored as the first argument to the command line. As these scripts run multiple complex workflows with a large number of sequencing endpoints at high coverages sequentially, run-times of around 1-3 h are to be expected when using four cores of a standard desktop CPU. Around 64 GB of RAM are recommended.



# Instructions for use

Import the package in your script:
```python
import dt4dds
```
At this point, you can set custom configuration options and enable logging, if desired:
```python
dt4dds.default_logging()                    # enable logging output
```
For the full documentation of configuration options, refer to [config.py](dt4dds/helpers/config.py).


## Representing sequences and oligo pools

Two datastructures are used to represent individual sequences and oligo pools: `Seq` and `SeqPool`. The `Seq` datastructure encapsulates a string representation of the oligo sequence, while the `SeqPool` datastructure encapsulates a `dict` associating each sequence with an abundance. Both datastructures are accessible in the `dt4dds.datastructures` namespace, however they are usually generated by other functions such as the synthesis process (see below) or those in [the generators module](dt4dds/helpers/generators.py).

### Convenience methods

Both the `Seq` and `SeqPool` datastructures implement convenience methods for common operations. In the case of the `SeqPool` datastructure for example, parameters of the oligo pool such as its mass concentration or the number of individual sequences may be accessed:
```python
print(pool.mass_concentration)  # in ng/uL
print(pool.n_sequences)
```
The pool may also be sampled by weight or volume:
```python
sampled_pool = pool.sample_by_mass(50)    # in ng
sampled_pool = pool.sample_by_volume(50)  # in uL
```
Or exported to other formats:
```python
pool.save_as_csv("my/path.csv")
pool.save_as_fastq("my/path.fastq")
```
For the full documentation of convenience methods, refer to [seq.py](dt4dds/datastructures/seq.py) and [seqpool.py](dt4dds/datastructures/seqpool.py).

## Building a workflow

All processes are provided in the `dt4dds.processes` namespace. In general, all process modules are implemented as classes that require instantiation, during which non-default experimental parameters can be optionally provided (see below). Process instances do not perform their actions directly, but must be executed with their `process` method. With the exception of the synthesis module, all modules' `process` methods require a `SeqPool` datastructure as input, representing the oligo pool the process is performed with. With the exception of the sequencing module, all modules' `process` methods return a `SeqPool` datastructure as output, representing the state of the oligo pool after the process.

Using the `SeqPool` datastructure returned by one process as the input to another process starts a workflow. Any `SeqPool` datastructure may be (re-)used in multiple processes to investigate the impacts of different workflows on the same oligo pool.

### Synthesis

The `dt4dds.processes.ArraySynthesis()` class models array-based oligo pool synthesis. It requires the design sequences (without primer regions) as input, and a `SeqPool` can be sampled by mass or moles after invocation of its `process()` method:
```python
# read sequences from reference
seq_list = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')

# set up synthesis, add non-default parameters as arguments if desired
array_synthesis = dt4dds.processes.ArraySynthesis()

# perform synthesis and sample to receive a SeqPool datastructure
array_synthesis.process(seq_list)
pool = array_synthesis.sample_by_mass(1)
```
For the full documentation of methods, refer to [array_synthesis.py](dt4dds/processes/array_synthesis.py). To see customization options and more details, see the [demo files](dt4dds/demos/) and [reproduction scripts](dt4dds/scripts/).


### PCR

The `dt4dds.processes.PCR()` class models amplification by polymerase chain reaction for a user-defined number of cycles. It requires a `SeqPool` as input, and yields a `SeqPool` representative of the amplified oligo pool by invocation of its `process()` method:
```python
# set up PCR with 20 cycles, add further non-default parameters as arguments if desired
pcr = dt4dds.processes.PCR(n_cycles=20)

# perform PCR and receive the oligo pool representation after PCR
post_PCR_pool = pcr.process(pool)
```
For the full documentation of methods, refer to [pcr.py](dt4dds/processes/pcr.py). To see customization options and more details, see the [demo files](dt4dds/demos/) and [reproduction scripts](dt4dds/scripts/).


### Aging

The `dt4dds.processes.Aging()` class models decay for a user-defined number of half-lives. It requires a `SeqPool` as input, and yields a `SeqPool` representative of the decayed oligo pool by invocation of its `process()` method:
```python
# set up aging for one half-life, add further non-default parameters as arguments if desired
aging = dt4dds.processes.Aging(n_halflives=1)

# perform aging and receive the oligo pool representation after decay
post_decay_pool = aging.process(pool)
```
For the full documentation of methods, refer to [aging.py](dt4dds/processes/aging.py). To see customization options and more details, see the [demo files](dt4dds/demos/) and [reproduction scripts](dt4dds/scripts/).

### Sequencing

The `dt4dds.processes.SBSSequencing()` class models sequencing-by-synthesis based on Illumina's sequencing platforms. It requires the `SeqPool` as input, and invocation of its `process()` method yields a user-defined number of sequencing reads saved at the provided output directory:
```python
# set up sequencing to a specific folder with 1M reads
sbs_sequencing = dt4dds.processes.SBSSequencing(output_directory="./sequencing/", n_reads=100000)

# perform sequencing, FASTQ files will be deposited in the output directory
sbs_sequencing.process(pool)
```
For the full documentation of methods, refer to [sbs_sequencing.py](dt4dds/processes/sbs_sequencing.py). To see customization options and more details, see the [demo files](dt4dds/demos/) and [reproduction scripts](dt4dds/scripts/).

## Setting custom parameters

All process classes support full customization of the experimental settings and error parameters used for simulation. Custom parameters must be supplied during instantiation, either as keyword arguments or by providing the appropriate, customized settings object from `dt4dds.settings`. For example, to customize some parameters for PCR, both options can be used:
```python
# set up PCR with keyword arguments
pcr = dt4dds.processes.PCR(n_cycles=20, polymerase_fidelity=5, efficiency_mean=0.9)

# set up PCR with a re-usable settings object
pcr_settings = dt4dds.settings.PCRSettings(n_cycles=20, polymerase_fidelity=5, efficiency_mean=0.9)
pcr = dt4dds.processes.PCR(pcr_settings)
```
For documentation of all settings and details about the settings objects in  `dt4dds.settings`, refer to [settings.py](dt4dds/settings/settings.py).


### Provided defaults
The `dt4dds.settings.defaults` module provides multiple default settings objects for some processes. For example, for the two different synthesis providers characterised in the manuscript:
```python
# set up electrochemical synthesis
synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_CustomArray()

# set up synthesis by material-deposition
synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Twist()
```
For documentation of all defaults, refer to [defaults.py](dt4dds/settings/defaults.py) and the configuration files in the [settings folder](dt4dds/settings/).


### Documentation of all parameters
For documentation of all settings, refer to [settings.py](dt4dds/settings/settings.py).


# License
This project is licensed under the GPLv3 license, see [here](LICENSE).