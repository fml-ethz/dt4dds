import rapidfuzz
import edlib
import collections
import ruamel.yaml
import Bio.SeqRecord
import dataclasses
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from . import config

class ErrorMap():

    __slots__ = [
        'n_reads', 'n_bases', 'n_aligned_bases', 'n_matches', 'n_bases_by_position', 'n_bases_by_refposition', 'length_distribution_aligned', 'length_distribution_windowed',
        'n_deletions', 'n_delevents', 'n_deletions_by_type', 'n_deletions_by_length', 'n_deletions_by_position', 'n_deletions_by_refposition', 'n_deletions_by_read', 'n_delevents_by_read', 'n_deletions_by_type_by_position', 'n_deletions_by_type_by_refposition',
        'n_insertions', 'n_insevents', 'n_insertions_by_type', 'n_insertions_by_length', 'n_insertions_by_position', 'n_insertions_by_refposition', 'n_insertions_by_read', 'n_insevents_by_read', 'n_insertions_by_type_by_position', 'n_insertions_by_type_by_refposition',
        'n_substitutions', 'n_subevents', 'n_substitutions_by_type', 'n_substitutions_by_length', 'n_substitutions_by_position', 'n_substitutions_by_refposition', 'n_substitutions_by_read', 'n_subevents_by_read', 'n_substitutions_by_type_by_position', 'n_substitutions_by_type_by_refposition',
        'n_reads_perfect', 'n_reads_noinsertion', 'n_reads_nodeletion', 'n_reads_nosubstitution'
    ]

    @classmethod
    def from_mapping(cls, mapping):
        
        instance = cls()

        instance.n_bases = 0
        instance.n_aligned_bases = len(mapping.aligned_ref_sequence)
        instance.length_distribution_aligned[len(mapping.aligned_ref_sequence)] += 1
        instance.length_distribution_windowed[len(mapping.windowed_read_sequence)] += 1
        instance.n_reads = 1

        for operation in mapping.operations:

            if operation[0] == 'equal':
                op_len = operation[4]-operation[3]
                instance.n_matches += op_len
                instance.n_bases += op_len
                for p_read in range(operation[3], operation[4]): 
                    instance.n_bases_by_position[p_read] += 1
                for p_ref in range(operation[1], operation[2]):
                    instance.n_bases_by_refposition[p_ref] += 1

            elif operation[0] == 'delete':
                op_len = operation[2]-operation[1]
                instance.n_deletions += op_len
                instance.n_bases += op_len
                instance.n_delevents += 1
                instance.n_deletions_by_length[op_len] += 1
                instance.n_deletions_by_position[operation[3]] += op_len
                for p_ref in range(operation[1], operation[2]):
                    instance.n_bases_by_refposition[p_ref] += 1
                    instance.n_deletions_by_refposition[p_ref] += 1
                    instance.n_deletions_by_type_by_position[mapping.full_ref_sequence[p_ref]][operation[3]] += 1
                    instance.n_deletions_by_type_by_refposition[mapping.full_ref_sequence[p_ref]][p_ref] += 1
                    instance.n_deletions_by_type[mapping.full_ref_sequence[p_ref]] += 1

            elif operation[0] == 'insert':
                op_len = operation[4]-operation[3]
                instance.n_insertions += op_len
                instance.n_bases += op_len
                instance.n_insevents += 1
                instance.n_insertions_by_length[op_len] += 1
                instance.n_insertions_by_refposition[max(operation[1]-1,0)] += op_len 
                for p_read in range(operation[3], operation[4]): 
                    instance.n_bases_by_position[max(p_read-1,0)] += 1
                    instance.n_insertions_by_position[max(p_read-1, 0)] += 1
                    instance.n_insertions_by_type_by_position[mapping.full_read_sequence[p_read]][max(p_read-1,0)] += 1
                    instance.n_insertions_by_type_by_refposition[mapping.full_read_sequence[p_read]][max(operation[1]-1,0)] += 1
                    instance.n_insertions_by_type[mapping.full_read_sequence[p_read]] += 1
                    
            elif operation[0] == 'replace':
                op_len = operation[4]-operation[3]
                instance.n_substitutions += op_len
                instance.n_bases += op_len
                instance.n_subevents += 1
                instance.n_substitutions_by_length[op_len] += 1
                for p_ref, p_read in zip(range(operation[1], operation[2]), range(operation[3], operation[4])): 
                    instance.n_bases_by_position[p_read] += 1
                    instance.n_bases_by_refposition[p_ref] += 1
                    instance.n_substitutions_by_position[p_read] += 1
                    instance.n_substitutions_by_refposition[p_ref] += 1
                    instance.n_substitutions_by_type[f"{mapping.full_ref_sequence[p_ref]}2{mapping.full_read_sequence[p_read]}"] += 1
                    instance.n_substitutions_by_type_by_position[f"{mapping.full_ref_sequence[p_ref]}2{mapping.full_read_sequence[p_read]}"][p_read] += 1
                    instance.n_substitutions_by_type_by_refposition[f"{mapping.full_ref_sequence[p_ref]}2{mapping.full_read_sequence[p_read]}"][p_ref] += 1

        if instance.n_insertions == 0: instance.n_reads_noinsertion = 1
        if instance.n_substitutions == 0: instance.n_reads_nosubstitution = 1
        if instance.n_deletions == 0: instance.n_reads_nodeletion = 1
        if instance.n_insertions == 0 and instance.n_substitutions == 0 and instance.n_deletions == 0: instance.n_reads_perfect = 1

        instance.n_deletions_by_read[instance.n_deletions] = 1
        instance.n_insertions_by_read[instance.n_insertions] = 1
        instance.n_substitutions_by_read[instance.n_substitutions] = 1

        instance.n_delevents_by_read[instance.n_delevents] = 1
        instance.n_insevents_by_read[instance.n_insevents] = 1
        instance.n_subevents_by_read[instance.n_subevents] = 1

        return instance



    def __init__(self):

        self.n_reads = 0
        self.length_distribution_aligned = collections.Counter()
        self.length_distribution_windowed = collections.Counter()
        self.n_bases = 0
        self.n_aligned_bases = 0
        self.n_bases_by_position = collections.Counter()
        self.n_bases_by_refposition = collections.Counter()
        self.n_matches = 0
        self.n_insertions = 0
        self.n_insevents = 0
        self.n_deletions = 0
        self.n_delevents = 0
        self.n_substitutions = 0
        self.n_subevents = 0

        self.n_insertions_by_type = collections.Counter({
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            'N': 0,
        })

        self.n_insertions_by_type_by_position = {
            'A': collections.Counter(),
            'C': collections.Counter(),
            'G': collections.Counter(),
            'T': collections.Counter(),
            'N': collections.Counter(),
        }
        self.n_insertions_by_type_by_refposition = {
            'A': collections.Counter(),
            'C': collections.Counter(),
            'G': collections.Counter(),
            'T': collections.Counter(),
            'N': collections.Counter(),
        }

        self.n_insertions_by_length = collections.Counter()

        self.n_insertions_by_position = collections.Counter()
        self.n_insertions_by_refposition = collections.Counter()

        self.n_insertions_by_read = collections.Counter()

        self.n_insevents_by_read = collections.Counter()


        self.n_deletions_by_type = collections.Counter({
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
        })

        self.n_deletions_by_type_by_position = {
            'A': collections.Counter(),
            'C': collections.Counter(),
            'G': collections.Counter(),
            'T': collections.Counter(),
            'N': collections.Counter(),
        }
        self.n_deletions_by_type_by_refposition = {
            'A': collections.Counter(),
            'C': collections.Counter(),
            'G': collections.Counter(),
            'T': collections.Counter(),
            'N': collections.Counter(),
        }

        self.n_deletions_by_length = collections.Counter()

        self.n_deletions_by_position = collections.Counter()
        self.n_deletions_by_refposition = collections.Counter()

        self.n_deletions_by_read = collections.Counter()

        self.n_delevents_by_read = collections.Counter()


        self.n_substitutions_by_type = collections.Counter({
            'A2C': 0,
            'A2G': 0,
            'A2T': 0,
            'C2A': 0,
            'C2G': 0,
            'C2T': 0,
            'G2A': 0,
            'G2C': 0,
            'G2T': 0,
            'T2A': 0,
            'T2C': 0,
            'T2G': 0,
        })

        self.n_substitutions_by_type_by_position = {
            'A2C': collections.Counter(),
            'A2G': collections.Counter(),
            'A2T': collections.Counter(),
            'C2A': collections.Counter(),
            'C2G': collections.Counter(),
            'C2T': collections.Counter(),
            'G2A': collections.Counter(),
            'G2C': collections.Counter(),
            'G2T': collections.Counter(),
            'T2A': collections.Counter(),
            'T2C': collections.Counter(),
            'T2G': collections.Counter(),
            'A2N': collections.Counter(),
            'C2N': collections.Counter(),
            'G2N': collections.Counter(),
            'T2N': collections.Counter(),
            'N2A': collections.Counter(),
            'N2C': collections.Counter(),
            'N2G': collections.Counter(),
            'N2T': collections.Counter(),
        }
        self.n_substitutions_by_type_by_refposition = {
            'A2C': collections.Counter(),
            'A2G': collections.Counter(),
            'A2T': collections.Counter(),
            'C2A': collections.Counter(),
            'C2G': collections.Counter(),
            'C2T': collections.Counter(),
            'G2A': collections.Counter(),
            'G2C': collections.Counter(),
            'G2T': collections.Counter(),
            'T2A': collections.Counter(),
            'T2C': collections.Counter(),
            'T2G': collections.Counter(),
            'A2N': collections.Counter(),
            'C2N': collections.Counter(),
            'G2N': collections.Counter(),
            'T2N': collections.Counter(),
            'N2A': collections.Counter(),
            'N2C': collections.Counter(),
            'N2G': collections.Counter(),
            'N2T': collections.Counter(),
        }
        
        self.n_substitutions_by_length = collections.Counter()

        self.n_substitutions_by_position = collections.Counter()
        self.n_substitutions_by_refposition = collections.Counter()

        self.n_substitutions_by_read = collections.Counter()

        self.n_subevents_by_read = collections.Counter()

        self.n_reads_perfect = 0
        self.n_reads_noinsertion = 0
        self.n_reads_nodeletion = 0
        self.n_reads_nosubstitution = 0


    def __add__(self, other):

        for element in self.__slots__:
            this = getattr(self, element)
            that = getattr(other, element)

            if type(this) == dict:
                for key in this:
                    this[key] += that[key]
            else:
                setattr(self, element, this + that)
        
        return self

    
    def save(self, filename):

        if self.n_bases == 0:
            logger.warning(f"Zero bases logged for file {filename}")
            return

        data = {elem: getattr(self, elem, None) for elem in self.__slots__}

        data['r_matches'] = self.n_matches/self.n_bases
        data['r_deletions'] = self.n_deletions/self.n_bases
        data['r_delevents'] = self.n_delevents/self.n_bases
        data['r_insertions'] = self.n_insertions/self.n_bases
        data['r_insevents'] = self.n_insevents/self.n_bases
        data['r_substitutions'] = self.n_substitutions/self.n_bases
        data['r_subevents'] = self.n_subevents/self.n_bases

        data['r_read_perfect'] = self.n_reads_perfect/self.n_reads
        data['r_read_nodeletion'] = self.n_reads_nodeletion/self.n_reads
        data['r_read_noinsertion'] = self.n_reads_noinsertion/self.n_reads
        data['r_read_nosubstitution'] = self.n_reads_nosubstitution/self.n_reads

        data['p_deletions_by_length'] = {l: n/max(sum(self.n_deletions_by_length.values()),1) for l, n in self.n_deletions_by_length.items()}
        data['p_insertions_by_length'] = {l: n/max(sum(self.n_insertions_by_length.values()),1) for l, n in self.n_insertions_by_length.items()}
        data['p_substitutions_by_length'] = {l: n/max(sum(self.n_substitutions_by_length.values()),1) for l, n in self.n_substitutions_by_length.items()}

        data['p_deletions_by_type'] = {b: n/max(sum(self.n_deletions_by_type.values()),1) for b, n in self.n_deletions_by_type.items()}
        data['p_insertions_by_type'] = {b: n/max(sum(self.n_insertions_by_type.values()),1) for b, n in self.n_insertions_by_type.items()}
        data['p_substitutions_by_type'] = {b: n/max(sum(self.n_substitutions_by_type.values()),1) for b, n in self.n_substitutions_by_type.items()}

        data['p_deletions_by_position'] = {p: n/max(self.n_bases_by_position[p],1) for p, n in self.n_deletions_by_position.items()}
        data['p_insertions_by_position'] = {p: n/max(self.n_bases_by_position[p],1) for p, n in self.n_insertions_by_position.items()}
        data['p_substitutions_by_position'] = {p: n/max(self.n_bases_by_position[p],1) for p, n in self.n_substitutions_by_position.items()}
        
        data['p_deletions_by_refposition'] = {p: n/max(self.n_bases_by_refposition[p],1) for p, n in self.n_deletions_by_refposition.items()}
        data['p_insertions_by_refposition'] = {p: n/max(self.n_bases_by_refposition[p],1) for p, n in self.n_insertions_by_refposition.items()}
        data['p_substitutions_by_refposition'] = {p: n/max(self.n_bases_by_refposition[p],1) for p, n in self.n_substitutions_by_refposition.items()}
        
        data['p_deletions_by_read'] = {p: n/max(sum(self.n_deletions_by_read.values()),1) for p, n in self.n_deletions_by_read.items()}
        data['p_insertions_by_read'] = {p: n/max(sum(self.n_insertions_by_read.values()),1) for p, n in self.n_insertions_by_read.items()}
        data['p_substitutions_by_read'] = {p: n/max(sum(self.n_substitutions_by_read.values()),1) for p, n in self.n_substitutions_by_read.items()}
        
        data['p_delevents_by_read'] = {p: n/max(sum(self.n_delevents_by_read.values()),1) for p, n in self.n_delevents_by_read.items()}
        data['p_insevents_by_read'] = {p: n/max(sum(self.n_insevents_by_read.values()),1) for p, n in self.n_insevents_by_read.items()}
        data['p_subevents_by_read'] = {p: n/max(sum(self.n_subevents_by_read.values()),1) for p, n in self.n_subevents_by_read.items()}
                
        data['p_length_distribution_windowed'] = {l: n/self.n_reads for l, n in self.length_distribution_windowed.items()}
        data['p_length_distribution_aligned'] = {l: n/self.n_reads for l, n in self.length_distribution_aligned.items()}

        data['p_deletions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_deletions_by_type_by_position.items()}
        data['p_insertions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_insertions_by_type_by_position.items()}
        data['p_substitutions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_substitutions_by_type_by_position.items()}

        data['p_deletions_by_refposition_by_type'] = {t: {p: n/max(self.n_bases_by_refposition[p],1) for p, n in d.items()} for t, d in self.n_deletions_by_type_by_refposition.items()}
        data['p_insertions_by_refposition_by_type'] = {t: {p: n/max(self.n_bases_by_refposition[p],1) for p, n in d.items()} for t, d in self.n_insertions_by_type_by_refposition.items()}
        data['p_substitutions_by_refposition_by_type'] = {t: {p: n/max(self.n_bases_by_refposition[p],1) for p, n in d.items()} for t, d in self.n_substitutions_by_type_by_refposition.items()}


        for key, item in data.items():
            if type(item) is collections.Counter:
                item = dict(item)
            if type(item) is dict:
                data[key] = dict(sorted({k: dict(sorted(dict(v).items())) if type(v) in (collections.Counter, dict) else v for k, v in item.items()}.items()))

        with open(filename, 'w') as f:
            yaml = ruamel.yaml.YAML()
            yaml.default_flow_style = False
            yaml.sort_keys = False
            yaml.dump(data, f)




@dataclasses.dataclass
class Mapping():

    read: Bio.SeqRecord.SeqRecord
    reference: Bio.SeqRecord.SeqRecord

    aln_ref_start: int
    aln_ref_end: int
    aln_read_start: int
    aln_read_end: int
    
    read_direction: int = None
    read_number: int = None


    def __post_init__(self):
    
        # initialize full window for read
        self.wnd_read_start: int = 0
        self.wnd_read_end: int = len(self.full_read_sequence)

        # restrict window based on trimming ends
        if config.trim_ends:
            self._trim_ends(config.trim_ends)

        # restrict window based on trimming adapters
        if config.trim_adapter:
            self._trim_adapter(config.trim_adapter)


    @classmethod
    def from_sam(cls, alignment, ref, read_direction=None, read_number=None):
        # create read object from alignment
        read = Bio.SeqRecord.SeqRecord(alignment.get_forward_sequence(), id=alignment.query_name)

        # if it's a reverse read, we need to correct the alignment start and end
        if alignment.is_reverse:
            ref_start, ref_end = len(ref.seq)-alignment.reference_end, len(ref.seq)-alignment.reference_start
            aln_start, aln_end = len(read.seq)-alignment.query_alignment_end, len(read.seq)-alignment.query_alignment_start
        
        # if it's a forward read, everything is fine
        else:
            ref_start, ref_end = alignment.reference_start, alignment.reference_end
            aln_start, aln_end = alignment.query_alignment_start, alignment.query_alignment_end

        # if a global alignment is performed, we want the full reference and the full read up to the adapter
        if not config.local:
            ref_start = 0
            ref_end = len(ref.seq)
            aln_start = 0
            aln_end = aln_end

        # create mapping
        return cls(
            read, 
            ref, 
            ref_start,
            ref_end,
            aln_start,
            aln_end,
            read_direction=read_direction,
            read_number=read_number,
        )

    @property
    def full_read_sequence(self):
        return str(self.read.seq)
    
    @property
    def windowed_read_sequence(self):
        return str(self.read.seq[self.wnd_read_start:self.wnd_read_end])
    
    @property
    def full_ref_sequence(self):
        return str(self.reference.seq)

    @property
    def aligned_ref_sequence(self):
        return self.full_ref_sequence[self.aln_ref_start:self.aln_ref_end]

    @property
    def aligned_read_sequence(self):
        return self.windowed_read_sequence[self.aln_read_start:self.aln_read_end] 

    @property
    def operations(self):

        # generate the operations necessary to transform the aligned reference sequence into the aligned read sequence
        operations = rapidfuzz.distance.Levenshtein.opcodes(self.aligned_ref_sequence, self.aligned_read_sequence).as_list()
        
        if config.local:
            # if the last operation is a deletion in the read, remove that operation because it is not necessary for local alignment
            if operations and operations[-1][0] == 'delete':
                op = operations.pop()
                self.aln_ref_end -= op[2]-op[1]
                self.aln_read_end -= op[4]-op[3]

            # if the first operation is a insertion in the read, remove that operation because it is not necessary for local alignment
            if operations and operations[0][0] == 'insert':
                op = operations.pop(0)
                self.aln_ref_start += op[2]-op[1]
                self.aln_read_start += op[4]-op[3]
                # recalculate the operations because removing the start operation changes the local positions in consequent operations
                operations = rapidfuzz.distance.Levenshtein.opcodes(self.aligned_ref_sequence, self.aligned_read_sequence).as_list()

            # adjust the operations to account for the offset of the alignment
            for i, op in enumerate(operations):
                op = [x for x in op]
                op[1] += self.aln_ref_start
                op[2] += self.aln_ref_start
                op[3] += self.aln_read_start
                op[4] += self.aln_read_start
                operations[i] = op
            
            # we need to add skips if the alignment starts after the beginning of the reference sequence
            if self.aln_ref_start != 0:
                operations.insert(0, ['front_skip', 0, self.aln_ref_start, 0, 0])

            # we need to add skips if the alignment ends before the end of the reference sequence
            if self.aln_ref_end != len(self.full_ref_sequence):
                operations.append(['back_skip', self.aln_ref_end, len(self.full_ref_sequence), self.aln_read_end, self.aln_read_end])

            # we need to add skips if the alignment starts after the beginning of the read sequence
            if self.aln_read_start != self.wnd_read_start:
                operations.insert(0, ['front_skip', 0, 0, self.wnd_read_start, self.aln_read_start])

            # we need to add skips if the alignment ends before the end of the read sequence
            if self.aln_read_end != self.wnd_read_end:
                operations.append(['back_skip', len(self.full_ref_sequence), len(self.full_ref_sequence), self.aln_read_end, self.wnd_read_end])
        
        else:
            # we need to add a delete operation at the beginning of the alignment if the alignment starts after the beginning of the reference sequence
            if self.aln_ref_start != 0:
                operations.insert(0, ['delete', 0, self.aln_ref_start, self.aln_read_start, self.aln_read_start])
                
            # we need to add a delete operation at the end of the alignment if the alignment ends before the end of the reference sequence
            if self.aln_ref_end != len(self.full_ref_sequence):
                operations.append(['delete', self.aln_ref_end, len(self.full_ref_sequence), self.aln_read_end, self.aln_read_end])

            # we need to add a insert operation at the beginning of the alignment if the alignment starts after the beginning of the read sequence
            if self.aln_read_start != self.wnd_read_start:
                operations.insert(0, ['insert', 0, 0, self.wnd_read_start, self.aln_read_start])

            # we move the window back if the alignment ends before the end of the read sequence, because we are likely reading into the sequencing adapter
            if self.aln_read_end != self.wnd_read_end:
                self.wnd_read_end = self.aln_read_end


        # common for both local and global alignment

        # introduce front skips if the windowed read sequence starts after the beginning of the read sequence
        if self.wnd_read_start != 0:
            operations.insert(0, ['front_skip', 0, 0, 0, self.wnd_read_start])

        # introduce back skips if the windowed read sequence ends before the end of the read sequence
        if self.wnd_read_end != len(self.full_read_sequence):
            operations.append(['back_skip', len(self.full_ref_sequence), len(self.full_ref_sequence), self.wnd_read_end, len(self.full_read_sequence)])

        return operations
        

    @property
    def aligned_distance(self):
        return rapidfuzz.distance.Levenshtein.distance(self.aligned_ref_sequence, self.aligned_read_sequence)

    @property
    def aligned_similarity(self):
        return rapidfuzz.distance.Levenshtein.normalized_similarity(self.aligned_ref_sequence, self.aligned_read_sequence)


    def __repr__(self):
        return f"Mapping(d={self.aligned_distance}, sim={self.aligned_similarity:.2f}, ref_id={self.reference.id}, ref_len={len(self.aligned_ref_sequence)}, read_len={len(self.aligned_read_sequence)}, op={self.operations})"


    def _trim_ends(self, trim_ends):
        
        # remove any bases from the end of the alignment that are in the trim_ends list
        to_remove = self.wnd_read_end-self.aln_read_end
        while (to_remove < self.wnd_read_end) and (self.windowed_read_sequence[-(to_remove+1)] in trim_ends):
            to_remove += 1

        # what needs to be cut from the alignment
        to_remove_aln = self.aln_read_end - (self.wnd_read_end-to_remove)

        # short-cut if the we are outside the alignment, we can cut freely
        if to_remove_aln <= 0:
            self.wnd_read_end -= to_remove
            return
        
        # we are inside the alignment and need to check overlap with the reference, starting from the innermost base
        pos_read = self.aln_read_end - to_remove_aln
        pos_ref = self.aln_ref_end - to_remove_aln

        # reassign any bases that were cut but match the reference, in case the end of the read contains the bases in the trim_ends list
        while (pos_read < self.wnd_read_end) and (pos_ref < len(self.full_ref_sequence)) and (self.full_read_sequence[pos_read] == self.full_ref_sequence[pos_ref]):
            pos_read += 1
            pos_ref += 1
            to_remove -= 1

        self.aln_read_end = max(pos_read, self.aln_read_start)
        self.aln_ref_end = max(pos_ref, self.aln_ref_start)
        self.wnd_read_end = max(self.wnd_read_end-to_remove, self.aln_read_end)


    def _trim_adapter(self, trim_adapter):
        read_ingress = 10
        
        # get the end of the read, including the bases before the end of alignment
        sub_read_start = max(0, self.aln_read_end-read_ingress)
        sub_read_end = min(len(self.full_read_sequence), self.aln_read_end+20)
        sub_read = self.full_read_sequence[sub_read_start:sub_read_end]

        # short-cut if the remaining read is too short for accurate alignment
        if len(sub_read) < 20:
            return

        # create a subsequence of the adapter
        adapter = trim_adapter[:10+1]

        # look for the adapter in the subsequence
        alignment = edlib.align(adapter, sub_read, mode="HW", task="locations")

        # if there's no alignment, do nothing
        alignment_length = alignment['locations'][-1][1]-alignment['locations'][-1][0]
        if alignment_length < 6 or alignment['editDistance'] > alignment_length/4.0:
            return
        
        # if there's an alignment, but it's already at the alignment end, do nothing
        offset = alignment['locations'][-1][0]
        if offset == read_ingress:
            return
        
        # check that the beginning of the adapter is present in the read
        if sub_read[offset:offset+4] != adapter[:4]:
            return
        
        # if there's an alignment, change the alignment and window positions
        self.aln_read_end = sub_read_start + offset
        self.wnd_read_end = self.aln_read_end


    @property
    def readmetrics(self):
        return {
            'read_id': self.read.id,
            'ref_id': self.reference.id,
            'read_direction': self.read_direction,
            'read_number': self.read_number,
            'len_ref': len(self.full_ref_sequence),
            'len_ref_aligned': len(self.aligned_ref_sequence),
            'len_read': len(self.full_read_sequence),
            'len_read_windowed': len(self.windowed_read_sequence),
            'len_read_aligned': len(self.aligned_read_sequence),
            'distance_aligned': self.aligned_distance,
            'similarity_aligned': self.aligned_similarity,
        }
    
    @property
    def breakmetrics(self):
        # look for and identify a possible 5' break
        break_5_5base, break_5_3base, break_5_pos = None, None, None
        # a 5' break is present if the read starts somewhere in the middle of the reference
        if self.aln_ref_start > 0:
            break_5_pos = self.aln_ref_start
            break_5_5base = str(self.full_ref_sequence[break_5_pos-1])
            break_5_3base = str(self.full_ref_sequence[break_5_pos])

        # look for and identify a possible 3' break
        break_3_5base, break_3_3base, break_3_pos = None, None, None
        # a 3' break is present if the read ends somewhere in the middle of the reference
        if self.aln_ref_end < len(self.full_ref_sequence):
            break_3_pos = self.aln_ref_end + 1 # +1 because the broken abasic position is after the last aligned base
            break_3_5base = str(self.full_ref_sequence[break_3_pos-1])
            break_3_3base = str(self.full_ref_sequence[break_3_pos]) if break_3_pos < len(self.full_ref_sequence) else None

        return {
            'read_id': self.read.id,
            'ref_id': self.reference.id,
            'read_direction': self.read_direction,
            'read_number': self.read_number,
            'break_5_pos': break_5_pos,
            'break_5_5base': break_5_5base,
            'break_5_3base': break_5_3base,
            'break_3_pos': break_3_pos,
            'break_3_5base': break_3_5base,
            'break_3_3base': break_3_3base,
            'adapter_5': self.aln_read_start > 1,
            'adapter_3': self.aln_read_end < len(self.full_read_sequence)-1,
        }


    def to_comparison_string(self):
        info = f">read_id:{self.read.id}|ref_id:{self.reference.id}|read_dir:{self.read_direction}|read_num:{self.read_number}|len_ref:{len(self.full_ref_sequence)}|len_ref_aligned:{len(self.aligned_ref_sequence)}|len_read:{len(self.windowed_read_sequence)}|len_read_aligned:{len(self.aligned_read_sequence)}|distance_aligned:{self.aligned_distance}|similarity_aligned:{self.aligned_similarity:.2f}"
        ref_string = ""
        read_string = ""

        ind_ref_start = False
        ind_ref_end = False
        ind_read_start = False
        ind_read_end = False

        for operation in self.operations:

            # add an indicator if the alignment window begins
            if not ind_ref_start and operation[1] == self.aln_ref_start and operation[0] != 'front_skip':
                ref_string += "[" 
                ind_ref_start = True
            if not ind_read_start and operation[3] == self.aln_read_start and operation[0] not in {'delete', 'front_skip'}:
                read_string += "[" 
                ind_read_start = True

            # for skips, just print the skipped sequence and blank the other
            if operation[0] == 'front_skip':
                op_len = max(operation[2]-operation[1], operation[4]-operation[3])
                ref_string += self.full_ref_sequence[operation[1]:operation[2]].ljust(op_len) + "|"
                read_string += self.full_read_sequence[operation[3]:operation[4]].ljust(op_len) + "|"
            elif operation[0] == 'back_skip':
                op_len = max(operation[2]-operation[1], operation[4]-operation[3])
                ref_string += "|" + self.full_ref_sequence[operation[1]:operation[2]].ljust(op_len)
                read_string += "|" + self.full_read_sequence[operation[3]:operation[4]].ljust(op_len)

            # for equal operations, print both reference and read directly
            elif operation[0] == 'equal':
                ref_string += self.full_ref_sequence[operation[1]:operation[2]]
                read_string += self.full_read_sequence[operation[3]:operation[4]]

            # for deletions, print the reference and add a hyphen to the read
            elif operation[0] == 'delete':
                op_len = operation[2]-operation[1]
                ref_string += self.full_ref_sequence[operation[1]:operation[2]]
                read_string += "-"*op_len

            # for insertions, print the read and add a hyphen to the reference
            elif operation[0] == 'insert':
                op_len = operation[4]-operation[3]
                read_string += self.full_read_sequence[operation[3]:operation[4]]
                ref_string += "-"*op_len   
            
            # for replacements, print both reference and read in lower case
            elif operation[0] == 'replace':
                op_len = operation[4]-operation[3]
                ref_string += self.full_ref_sequence[operation[1]:operation[2]].lower()
                read_string += self.full_read_sequence[operation[3]:operation[4]].lower()

            # also add an indicator if the alignment window is finished
            if not ind_ref_end and operation[2] == self.aln_ref_end and operation[0] not in {'insert', 'back_skip'}:
                ref_string += "]"
                ind_ref_end = True
            if not ind_read_end and operation[4] == self.aln_read_end and operation[0] != 'back_skip':
                read_string += "]"
                ind_read_end = True

        # check that the alignment is complete
        replacetable = str.maketrans({'|': '', '[': '', ']': '', '-': '', ' ': ''})
        if ref_string.translate(replacetable).upper() != self.full_ref_sequence:
            logger.critical(f"Reference sequence does not match alignment for {self.read.id} on {self.reference.id}.\n{ref_string}\n{read_string}\nFull: {self.full_ref_sequence}\nRec.: {ref_string.translate(replacetable).upper()}\nRead: {self.full_read_sequence}\nOps: {self.operations}")
        if read_string.translate(replacetable).upper() != self.full_read_sequence:
            logger.critical(f"Read sequence does not match alignment for {self.read.id} on {self.reference.id}.\n{ref_string}\n{read_string}\nFull: {self.full_read_sequence}\nRec.: {read_string.translate(replacetable).upper()}\nRef: {self.full_ref_sequence}\nOps: {self.operations}")

        return "\n".join([info, ref_string, read_string, "\n"])
