import yaml
import dataclasses
from typing import List, Dict, Tuple, Union



@dataclasses.dataclass
class AbstractSettings():
    """ Class for keeping track of all settings for a process step. """

    #
    # general methods
    #
    
    @classmethod
    def from_yaml(cls, filepath):
        with open(filepath, 'r') as f:
            kwargs = yaml.safe_load(f)
        return cls(**kwargs)


    def write_yaml(self, filepath):
        with open(filepath, 'w') as f:
            yaml.dump(dataclasses.asdict(self), f, default_flow_style=False)

    
    def __call__(self, **kwds):
        return dataclasses.replace(self, **kwds)


    def update(self, **new):
        """ Sets the settings supplied in the dictionary new. """
        for key, value in new.items():
            if hasattr(self, key):
                setattr(self, key, value)




@dataclasses.dataclass
class SynthesisSettings(AbstractSettings):
    """ Settings for Sequencing. """

    # general synthesis settings

    per_oligo_scale: float = 1.0
    """ Per oligo production scale in fmol. """


    # distribution settings

    oligo_distribution_type: str = 'lognormal'
    """ Type of the distribution for oligo count ratio, e.g. (oligo counts per sequence)/(mean oligo count). """

    oligo_distribution_params: Dict[str, float] = dataclasses.field(default_factory = lambda: {'mean': 0.0, 'sigma': 0.30})
    """ Parameters of the distribution for the oligo count ratio. """
    

    # error-related settings

    substitution_rate: float = 0.0
    """ Mean substitution rate per nucleotide. """

    substitution_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": {"C": 1/12, "G": 1/12, "T": 1/12},
        "C": {"A": 1/12, "G": 1/12, "T": 1/12},
        "G": {"A": 1/12, "C": 1/12, "T": 1/12},
        "T": {"A": 1/12, "C": 1/12, "G": 1/12},
    })
    """ Substition biases expressed as conditional error probability from one base to another. Sum of all values gives unity. """

    substitution_length_bias: str = None
    """ Length bias of the substitutions. None for no bias, or 'points'. """

    substitution_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the substitutions. """

    insertion_rate: float = 0.0
    """ Mean insertion rate per nucleotide. """

    insertion_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    })
    """ Insertion biases expressed as conditional error probability. Sum of all values gives unity. """

    insertion_length_bias: str = None
    """ Length bias of the insertions. None for no bias, or 'points'. """

    insertion_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the insertions. """

    deletion_rate: float = 0.0005695
    """ Mean deletion rate per nucleotide. """

    deletion_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": 0.2468,
        "C": 0.2362,
        "G": 0.2669,
        "T": 0.2500,
    })
    """ Deletion biases expressed as conditional error probability. Sum of all values gives unity. """
    
    deletion_read_bias: str = None
    """ Repeat bias of the deletions. None for no bias, or 'points'. """

    deletion_read_bias_params: List[float] = dataclasses.field(default_factory = lambda: None)
    """ Distribution of deletion lengths expressed as conditional error probability. Sum of all values gives unity. """
    
    deletion_repeat_bias: str = None
    """ Repeat bias of the deletions. None for no bias, or 'points'. """

    deletion_repeat_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the repeat bias of the deletions. """

    deletion_length_bias: str = None
    """ Length bias of the deletions. None for no bias, or 'points'. """

    deletion_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the deletions. """



@dataclasses.dataclass
class PCRSettings(AbstractSettings):
    """ Settings for PCR. """

    # general PCR related

    n_cycles: int = 20
    """ Number of PCR cycles. """

    volume: float = 20
    """ Total volume of the PCR mix in each well in uL. """

    template_volume: float = 5
    """ Volume of the template added to each well in uL. """

    primers: List[str] = dataclasses.field(default_factory = lambda: ("", ""))
    """ Tuple of primer sequences (forward, reverse). """


    # efficiency related

    efficiency_mean: float = 0.9
    """ Mean efficiency used for PCR. """


    # polymerase related

    polymerase_fidelity: float = 1.0
    """ Fidelity of polymerase relative to Taq. """

    polymerase_basesubstitutionrate: float = 1.09e-4
    """ Substitution rate of polymerase in subs/base/doubling for fidelity=1 (Taq). """

    polymerase_substitutionbias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": {"C": 0.0147, "G": 0.3028, "T": 0.0630},
        "C": {"A": 0.0150, "G": 0.0071, "T": 0.0975},
        "G": {"A": 0.0975, "C": 0.0071, "T": 0.0150},
        "T": {"A": 0.0630, "C": 0.3028, "G": 0.0147},
    })
    """ Substition biases of the polymerase expressed as conditional error probability from one base to another. Sum of all values gives unity. """


    # technical parameters
    
    primer_minoverlap: int = 15
    """ Minimum base pair overlap for primers. """




@dataclasses.dataclass
class SequencingSettings(AbstractSettings):
    """ Settings for Sequencing. """

    # general Sequencing related

    n_reads: int = 100000
    """ Number of reads. Set to desired sequencing coverage times pool size to get specific coverage depth. """

    read_mode: str = 'paired-end'
    """ Read mode of the sequencer. Either single-end or paired-end. """

    read_length: int = 150
    """ Read length of the sequencer. If read mode is paired-end, this is for each paired-end read individually. """


    # output related

    output_directory: str = './'
    """ Directory which all fastq files are written to. """


    # primer related

    primers_read: Tuple[str] = dataclasses.field(default_factory = lambda: (
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    ))
    """ Available read primers from 5' to 3'. """


    # error related

    substitution_rates: List[float] = dataclasses.field(default_factory = lambda: [0.001123, 0.002500])
    """ Substitution rate per base of sequencer in first and second read. """

    substitution_bias: List[Dict[str, Dict[str, float]]] = dataclasses.field(default_factory = lambda: [{
        "A": {"C": 0.0029, "G": 0.2065, "T": 0.1684},
        "C": {"A": 0.0246, "G": 0.0139, "T": 0.1594},
        "G": {"A": 0.1761, "C": 0.0184, "T": 0.0377},
        "T": {"A": 0.0203, "C": 0.1060, "G": 0.0657},
    },
    {
        "A": {"C": 0.0027, "G": 0.1147, "T": 0.3691},
        "C": {"A": 0.0174, "G": 0.0080, "T": 0.1038},
        "G": {"A": 0.0957, "C": 0.0421, "T": 0.0315},
        "T": {"A": 0.0227, "C": 0.1144, "G": 0.0777},
    }])
    """ Substition biases of the sequencer expressed as conditional error probability from one base to another. Sum of all values gives unity. """

    substitution_length_bias: str = 'points'
    """ Length bias of the substitutions. None for no bias, or 'points'. """

    substitution_length_bias_params: List[Union[None, List[float]]] = dataclasses.field(default_factory = lambda: [
        [1.9152769677330956, 1.7945321644707304, 1.6541113990588867, 1.538349039641976, 1.4789742148335365, 1.4041497563851684, 1.2209683441713857, 0.9857004935128433, 0.8738697525194686, 0.8353477227670248, 0.8162184068536342, 0.7754188120194149, 0.7607737831229126, 0.7465445553396249, 0.7124274668147668, 0.7020980050655049, 0.6843336093513959, 0.6599526606165358, 0.655701177635213, 0.6465398540603987, 0.6482405213116028, 0.6502575652308121, 0.6613973355932354, 0.6798020248190578, 0.7040217534179172, 0.726742630719294, 0.7475960914289453, 0.7778283578386342, 0.8001923065927612, 0.82436176334884, 0.8319266076345705, 0.8309729817259107, 0.8410600582716268, 0.8402938203377556, 0.8272828744836881, 0.8190913546773111, 0.8096513123416752, 0.8080012787192398, 0.8024735754346608, 0.7982950474506363, 0.7939931954057067, 0.794015529934872, 0.7953641675339886, 0.8025953863472424, 0.7958286916306596, 0.8110676473210389, 0.8163157444879552, 0.8121587351081074, 0.8094554960795335, 0.8145946633921047, 0.8336217457665669, 0.8522665108818539, 0.8362105520391567, 0.8410526774532135, 0.8478763121121393, 0.832197309904992, 0.8386641208175477, 0.8418490355130668, 0.8504787473859995, 0.8500094691488297, 0.8535086740937461, 0.8497033323334158, 0.864863553513162, 0.8752063391927191, 0.895703186986697, 0.9035762815137811, 0.9091517993099093, 0.9148584222259435, 0.9048626786498333, 0.8958209152122546, 0.8915051726073993, 0.8932200978500848, 0.9083691491717437, 0.90169685067718, 0.8950736414796894, 0.9068486655604732, 0.9440413849358568, 0.9560718448211081, 0.9816399018645293, 1.0139564899514353, 1.0253096777027721, 1.0248625533736462, 1.025595002947089, 1.0363749780596363, 1.0505694360254838, 1.064103575633507, 1.0337173613419572, 1.0418010594643927, 1.0523581152496313, 1.053886856035225, 1.0604668082337911, 1.063537028518198, 1.0648732272262385, 1.0853124361309439, 1.0878182398769292, 1.0969994719623415, 1.102221095292257, 1.1094596960884477, 1.1225266038668549, 1.1469787123192878, 1.1615861712121585, 1.2030601020603404, 1.2347810896081075, 1.2479258268089088, 1.269619076848689, 1.3028538812930204, 1.3508289044279016, 1.3819195660900423, 1.4399146481874885, 1.4786121856974266, 1.5410582056860536, 1.5885585542051863, 1.636780901734753, 1.6671807480259677, 1.7133116588454904, 1.7403117429525743, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429, 1.7728532334000429],
        [1.308580648333418, 1.2841700639612994, 1.2397737671088362, 1.2044573326951291, 1.1807810811392523, 1.1562964691363666, 1.0869645733919326, 1.028267777047154, 0.9705590116658364, 0.9165571357767837, 0.8706943366847444, 0.8211101189832697, 0.7932093378070888, 0.7668881841303739, 0.749518498815188, 0.7298271630710012, 0.7164558031572197, 0.7092528848955757, 0.70388252811112, 0.7005682064234784, 0.7039188745623709, 0.7067377228804763, 0.7038782290359263, 0.7090811372687946, 0.701856823548073, 0.7088741072750607, 0.7135890919936948, 0.7249499258556854, 0.7320301808388315, 0.7459894952131761, 0.7481497749333993, 0.7551905714562154, 0.7691430530625505, 0.7739149499912142, 0.7756082605200193, 0.7785467156578694, 0.7861520177759128, 0.7861384784367325, 0.7917974923759014, 0.7957479551165788, 0.7959062986205528, 0.8037487613024167, 0.8037304000420518, 0.802509857434918, 0.8113936380477146, 0.8183477123570493, 0.8162989804755628, 0.8244895878955846, 0.8312580824312576, 0.8369161456020497, 0.8517774484979805, 0.8601241213478212, 0.8617153672057754, 0.8718584128867127, 0.8710025507396273, 0.8619560808122714, 0.8745429415384502, 0.8712995606026342, 0.8707199801533025, 0.8778521073176042, 0.8897863169257194, 0.8893532047476093, 0.8982812543081331, 0.9083663357133692, 0.9191708741046236, 0.931527712832471, 0.9344188162247661, 0.9459608032916877, 0.9512606698342007, 0.9543309255187231, 0.9530492514138376, 0.9644874573023223, 0.9813408619823963, 0.9799551980287062, 0.9913969389964906, 1.006486403720898, 1.0237798109443446, 1.0421209972966092, 1.059209778652834, 1.0625379885322648, 1.0738513759246886, 1.0835057268828363, 1.0945896477807948, 1.114279528970601, 1.1208766122620835, 1.1299113419964526, 1.1301014441614927, 1.1343393417884706, 1.1397114653112377, 1.1525289447199765, 1.1536935195020885, 1.1568227412736742, 1.1661654105264019, 1.1773640695970609, 1.1885731342435455, 1.2094534559536332, 1.2320464073740824, 1.243231436984742, 1.26515483005746, 1.2839127918991446, 1.3031638561309205, 1.3229802792215013, 1.32801068539589, 1.3231260132465335, 1.3369020520151895, 1.330810893219669, 1.3377525692119363, 1.3584342319630345, 1.3807938488390143, 1.400122623950993, 1.420712111763282, 1.4447184353741362, 1.460589789803608, 1.4862068005582383, 1.49936664701115, 1.5235108827057993, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726, 1.539305636561726]
    ])
    """ Parameters for the length bias of the substitutions. """

    insertion_rates: List[float] = dataclasses.field(default_factory = lambda: [0.0, 0.0])
    """ Insertion rate per base of sequencer in first and second read. """

    insertion_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    })
    """ Insertion biases of the sequencer expressed as conditional error probability. Sum of all values gives unity. """

    insertion_length_bias: str = None
    """ Length bias of the insertions. None for no bias, or 'points'. """

    insertion_length_bias_params: List[Union[None, List[float]]] = dataclasses.field(default_factory = lambda: [None, None])
    """ Parameters for the length bias of the insertions. """

    deletion_rates: List[float] = dataclasses.field(default_factory = lambda: [0.0, 0.0])
    """ Deletion rate per base of sequencer in first and second read. """

    deletion_bias: Dict[str, float] = dataclasses.field(default_factory = lambda: {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    })
    """ Deletion biases of the sequencer expressed as conditional error probability. Sum of all values gives unity. """

    deletion_length_bias: str = None
    """ Length bias of the deletions. None for no bias, or 'points'. """

    deletion_length_bias_params: List[Union[None, List[float]]] = dataclasses.field(default_factory = lambda: [None, None])
    """ Parameters for the length bias of the deletions. """


    # technical parameters
    
    primer_minoverlap: int = 15
    """ Minimum base pair overlap for primers. """



@dataclasses.dataclass
class AgingSettings(AbstractSettings):
    """ Settings for Aging. """

    # decay related

    method: str = "fixed"
    """ Method of calculating decay, either 'fixed' or 'arrhenius'. """

    fixed_decay_ratio: float = 0.99
    """ Ratio of reads which are decayed if method == 'fixed'. """

    arrhenius_T: float = 20
    """ Temperature of storage in Celsius if method == 'fixed'. """

    arrhenius_t: float = 1
    """ Time of storage in years if method == 'fixed'. """

    arrhenius_RH: float = 0.5
    """ Relative humidity of storage in decimal if method == 'fixed'. """


    # error related

    substitution_rate: float = 0.000164
    """ Substitution rate per base per half-time of decay. """

    substitution_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": {"C": 0.0, "G": 0.0, "T": 0.0},
        "C": {"A": 0.0, "G": 0.0, "T": 1.0},
        "G": {"A": 0.0, "C": 0.0, "T": 0.0},
        "T": {"A": 0.0, "C": 0.0, "G": 0.0},
    })
    """ Substition bias expressed as conditional error probability from one base to another. Sum of all values gives unity. """

    substitution_length_bias: str = None
    """ Length bias of the substitutions. None for no bias, or 'points'. """

    substitution_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the substitutions. """

    insertion_rate: float = 0.0
    """ Insertion rate per base per half-time of decay. """

    insertion_bias: Dict[str, Dict[str, float]] = dataclasses.field(default_factory = lambda: {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    })
    """ Insertion bias expressed as conditional error probability. Sum of all values gives unity. """

    insertion_length_bias: str = None
    """ Length bias of the insertions. None for no bias, or 'points'. """

    insertion_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the insertions. """

    deletion_rate: float = 0.0
    """ Deletion rate per base per half-time of decay. """

    deletion_bias: Dict[str, float] = dataclasses.field(default_factory = lambda: {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    })
    """ Deletion bias expressed as conditional error probability. Sum of all values gives unity. """

    deletion_length_bias: str = None
    """ Length bias of the deletions. None for no bias, or 'points'. """

    deletion_length_bias_params: Union[None, List[float]] = dataclasses.field(default_factory = lambda: None)
    """ Parameters for the length bias of the deletions. """


    # technical parameters
    
    arrhenius_Ea: int = 155000
    """ Activation energy of arrhenius law for DNA decay in J/mol. """



@dataclasses.dataclass
class PropertiesSettings(AbstractSettings):
    """ Settings for general sequence properties. """


    # efficiency settings

    efficiency_distribution: str = 'normal'
    """ Type of the distribution for efficiency distribution. """

    efficiency_params: Dict[str, float] = dataclasses.field(default_factory = lambda: {'loc': 1.0, 'scale': 0.0051})
    """ Parameters of the distribution for the efficiency. """
