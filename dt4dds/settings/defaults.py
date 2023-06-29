import os
DEFAULT_DIR = os.path.dirname(os.path.abspath(__file__))

from .settings import AgingSettings, SynthesisSettings, PCRSettings, SequencingSettings, PropertiesSettings


def iSeq100(**kwargs):
    settings = SequencingSettings.from_yaml(os.path.join(DEFAULT_DIR, 'iseq100.yaml'))
    settings.update(**kwargs)
    return settings

def PCR(**kwargs):
    settings = PCRSettings.from_yaml(os.path.join(DEFAULT_DIR, 'pcr_taq.yaml'))
    settings.update(**kwargs)
    return settings

def PCR_Taq(**kwargs):
    settings = PCRSettings.from_yaml(os.path.join(DEFAULT_DIR, 'pcr_taq.yaml'))
    settings.update(**kwargs)
    return settings

def PCR_HiFi(**kwargs):
    settings = PCRSettings.from_yaml(os.path.join(DEFAULT_DIR, 'pcr_hifi.yaml'))
    settings.update(**kwargs)
    return settings

def PCR_Q5(**kwargs):
    settings = PCRSettings.from_yaml(os.path.join(DEFAULT_DIR, 'pcr_Q5.yaml'))
    settings.update(**kwargs)
    return settings

def PCR_Exo(**kwargs):
    settings = PCRSettings.from_yaml(os.path.join(DEFAULT_DIR, 'pcr_exo.yaml'))
    settings.update(**kwargs)
    return settings

def ArraySynthesis_CustomArray(**kwargs):
    settings = SynthesisSettings.from_yaml(os.path.join(DEFAULT_DIR, 'array_synthesis_customarray.yaml'))
    settings.update(**kwargs)
    return settings

def ArraySynthesis_CustomArray_GCfix(**kwargs):
    settings = SynthesisSettings.from_yaml(os.path.join(DEFAULT_DIR, 'array_synthesis_customarray_GCfix.yaml'))
    settings.update(**kwargs)
    return settings

def ArraySynthesis_CustomArray_GCall(**kwargs):
    settings = SynthesisSettings.from_yaml(os.path.join(DEFAULT_DIR, 'array_synthesis_customarray_GCall.yaml'))
    settings.update(**kwargs)
    return settings

def ArraySynthesis_Twist(**kwargs):
    settings = SynthesisSettings.from_yaml(os.path.join(DEFAULT_DIR, 'array_synthesis_twist.yaml'))
    settings.update(**kwargs)
    return settings

def Aging(**kwargs):
    settings = AgingSettings.from_yaml(os.path.join(DEFAULT_DIR, 'fixed_aging.yaml'))
    settings.update(**kwargs)
    return settings

def SequenceProperties(**kwargs):
    settings = PropertiesSettings.from_yaml(os.path.join(DEFAULT_DIR, 'properties.yaml'))
    settings.update(**kwargs)
    return settings