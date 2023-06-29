import rapidfuzz
import collections
import ruamel.yaml
import edlib
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ErrorMap():

    __slots__ = [
        'n_reads', 'n_bases', 'n_matches', 'n_bases_by_position',
        'n_deletions', 'n_delevents', 'n_deletions_by_type', 'n_deletions_by_length', 'n_deletions_by_position', 'n_deletions_by_refposition', 'n_deletions_by_read', 'n_delevents_by_read', 'n_deletions_by_type_by_position', 'n_deletions_by_type_by_refposition',
        'n_insertions', 'n_insevents', 'n_insertions_by_type', 'n_insertions_by_length', 'n_insertions_by_position', 'n_insertions_by_refposition', 'n_insertions_by_read', 'n_insevents_by_read', 'n_insertions_by_type_by_position', 'n_insertions_by_type_by_refposition',
        'n_substitutions', 'n_subevents', 'n_substitutions_by_type', 'n_substitutions_by_length', 'n_substitutions_by_position', 'n_substitutions_by_refposition', 'n_substitutions_by_read', 'n_subevents_by_read', 'n_substitutions_by_type_by_position', 'n_substitutions_by_type_by_refposition',
        'n_reads_perfect', 'n_reads_noinsertion', 'n_reads_nodeletion', 'n_reads_nosubstitution'
    ]

    @classmethod
    def from_mapping(cls, mapping):
        
        instance = cls()

        instance.n_bases = len(mapping.ref_sequence)
        instance.n_reads = 1

        for operation in mapping.operations:

            if operation[0] == 'equal':
                op_len = operation[4]-operation[3]
                instance.n_matches += op_len
                for p_read in range(operation[3], operation[4]): 
                    instance.n_bases_by_position[p_read] += 1

            elif operation[0] == 'delete':
                op_len = operation[2]-operation[1]
                instance.n_deletions += op_len
                instance.n_delevents += 1
                instance.n_deletions_by_length[op_len] += 1
                instance.n_deletions_by_position[operation[3]] += op_len
                for p_ref in range(operation[1], operation[2]):
                    instance.n_deletions_by_refposition[p_ref] += 1
                    instance.n_deletions_by_type_by_position[mapping.ref_sequence[p_ref]][operation[3]] += 1
                    instance.n_deletions_by_type_by_refposition[mapping.ref_sequence[p_ref]][p_ref] += 1
                    instance.n_deletions_by_type[mapping.ref_sequence[p_ref]] += 1

            elif operation[0] == 'insert':
                op_len = operation[4]-operation[3]
                instance.n_insertions += op_len
                instance.n_insevents += 1
                instance.n_insertions_by_length[op_len] += 1
                instance.n_insertions_by_refposition[operation[1]] += op_len 
                for p_read in range(operation[3], operation[4]): 
                    instance.n_bases_by_position[p_read] += 1
                    instance.n_insertions_by_position[p_read] += 1
                    instance.n_insertions_by_type_by_position[mapping.aligned_read_sequence[p_read]][p_read] += 1
                    instance.n_insertions_by_type_by_refposition[mapping.aligned_read_sequence[p_read]][operation[1]] += 1
                    instance.n_insertions_by_type[mapping.aligned_read_sequence[p_read]] += 1
                    
            elif operation[0] == 'replace':
                op_len = operation[4]-operation[3]
                instance.n_substitutions += op_len
                instance.n_subevents += 1
                instance.n_substitutions_by_length[op_len] += 1
                for p_ref, p_read in zip(range(operation[1], operation[2]), range(operation[3], operation[4])): 
                    instance.n_bases_by_position[p_read] += 1
                    instance.n_substitutions_by_position[p_read] += 1
                    instance.n_substitutions_by_refposition[p_ref] += 1
                    instance.n_substitutions_by_type[f"{mapping.ref_sequence[p_ref]}2{mapping.aligned_read_sequence[p_read]}"] += 1
                    instance.n_substitutions_by_type_by_position[f"{mapping.ref_sequence[p_ref]}2{mapping.aligned_read_sequence[p_read]}"][p_read] += 1
                    instance.n_substitutions_by_type_by_refposition[f"{mapping.ref_sequence[p_ref]}2{mapping.aligned_read_sequence[p_read]}"][p_ref] += 1

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

        self.n_bases = 0
        self.n_bases_by_position = collections.Counter()
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
                new_d = {}
                keys = set(this.keys()).union(set(that.keys()))
                for key in keys:
                    new_d[key] = this.get(key, collections.Counter()) + that.get(key, collections.Counter())
                setattr(self, element, new_d)
            else:
                new_value = this + that
                setattr(self, element, new_value)
        
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
        
        data['p_deletions_by_refposition'] = {p: n/self.n_reads for p, n in self.n_deletions_by_refposition.items()}
        data['p_insertions_by_refposition'] = {p: n/self.n_reads for p, n in self.n_insertions_by_refposition.items()}
        data['p_substitutions_by_refposition'] = {p: n/self.n_reads for p, n in self.n_substitutions_by_refposition.items()}
        
        data['p_deletions_by_read'] = {p: n/max(sum(self.n_deletions_by_read.values()),1) for p, n in self.n_deletions_by_read.items()}
        data['p_insertions_by_read'] = {p: n/max(sum(self.n_insertions_by_read.values()),1) for p, n in self.n_insertions_by_read.items()}
        data['p_substitutions_by_read'] = {p: n/max(sum(self.n_substitutions_by_read.values()),1) for p, n in self.n_substitutions_by_read.items()}
        
        data['p_delevents_by_read'] = {p: n/max(sum(self.n_delevents_by_read.values()),1) for p, n in self.n_delevents_by_read.items()}
        data['p_insevents_by_read'] = {p: n/max(sum(self.n_insevents_by_read.values()),1) for p, n in self.n_insevents_by_read.items()}
        data['p_subevents_by_read'] = {p: n/max(sum(self.n_subevents_by_read.values()),1) for p, n in self.n_subevents_by_read.items()}
        
        data['p_deletions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_deletions_by_type_by_position.items()}
        data['p_insertions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_insertions_by_type_by_position.items()}
        data['p_substitutions_by_position_by_type'] = {t: {p: n/max(self.n_bases_by_position[p],1) for p, n in d.items()} for t, d in self.n_substitutions_by_type_by_position.items()}

        data['p_deletions_by_refposition_by_type'] = {t: {p: n/self.n_reads for p, n in d.items()} for t, d in self.n_deletions_by_type_by_refposition.items()}
        data['p_insertions_by_refposition_by_type'] = {t: {p: n/self.n_reads for p, n in d.items()} for t, d in self.n_insertions_by_type_by_refposition.items()}
        data['p_substitutions_by_refposition_by_type'] = {t: {p: n/self.n_reads for p, n in d.items()} for t, d in self.n_substitutions_by_type_by_refposition.items()}


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






class Mapping():

    @classmethod
    def for_phix(cls, read, ref, read_number=None, trim_seq=None):

        instance = cls(read, ref, read_number=read_number, trim_seq=trim_seq)

        read_end = len(instance.full_read_sequence)
        if trim_seq:
            ed = edlib.align(trim_seq, instance.full_read_sequence, mode = "HW", task = "locations")
            if ed['editDistance'] < len(trim_seq)/5.0 and ed['locations'][-1][0] > 0:
                read_end = ed['locations'][-1][0]
               

        read_seq_wo_adapter = instance.full_read_sequence[:read_end]
        instance.align_start = 0
        instance.align_end = read_end

        ed = edlib.align(read_seq_wo_adapter, instance.ref_sequence, mode = "HW", task = "locations")

        instance.ref_start = ed['locations'][-1][0]
        instance.ref_end = ed['locations'][-1][1]+1

        instance.operations = rapidfuzz.distance.Levenshtein.opcodes(instance.ref_sequence, instance.aligned_read_sequence).as_list()

        return instance



    def __init__(self, read, ref, read_number=None, trim_seq=None, read_length=None):

        self.read = read
        self.full_read_sequence = str(read.seq)
        self.read_number = read_number
        self.read_length = read_length

        self.ref = ref
        self.full_ref_sequence = str(ref.seq)
        self.ref_start = 0
        self.ref_end = len(self.full_ref_sequence)

        # generate a prototype alignment
        self._generate_initial_alignment()

        # trim the adapter from end
        if trim_seq: 
            self._trim_adapter(trim_seq)

        # refine alignment 
        self._trim_excess_insertions()
        self._convert_initial_deletions()
        self._trim_initial_insertions()
        
        # remove deletions which occur because of exceeding read window
        ref_length = self.ref_end-self.ref_start
        if self.read_length and ref_length >= self.read_length:
            self._trim_excess_deletions()

        # remove gapped deletions
        self.operations = rapidfuzz.distance.Levenshtein.opcodes(self.ref_sequence, self.aligned_read_sequence).as_list()
        self._close_deletion_gaps()


    def _generate_initial_alignment(self):
        ed = edlib.align(self.ref_sequence, self.full_read_sequence, mode = "HW", task = "locations")
        self.align_start = ed['locations'][-1][0]
        self.align_end = ed['locations'][-1][1]+1



    def _trim_adapter(self, sequence):
        
        ed = edlib.align(sequence, self.full_read_sequence, mode = "HW", task = "locations")
        if self.align_end > ed['locations'][-1][0] & ed['locations'][-1][0] > 0.8*self.ref_end:
            if ed['editDistance'] < len(sequence)/5.0 and ed['locations'][-1][0] > 0:
                self.align_end = ed['locations'][-1][0]


    def _convert_initial_deletions(self):
        # convert initial deletions due to misaligned alignment window back into substitutions
        operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, self.aligned_read_sequence).as_list()
        i = 0
        while operations and operations[0][0] == 'delete' and self.align_start-i > 0:
            i += 1
            aligned_read_sequence = self.full_read_sequence[self.align_start-i:self.align_end]
            operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, aligned_read_sequence).as_list()    
        
        self.align_start -= i


    def _trim_initial_insertions(self):
        # trim excess inserts from start of mapping, these are mapping random noise because of misalignment
        operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, self.aligned_read_sequence).as_list()
        pos = operations[0][1] if operations else None
        i = 0
        while operations and operations[0][0] == 'insert':
            i += 1
            aligned_read_sequence = self.full_read_sequence[self.align_start+i:self.align_end]
            operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, aligned_read_sequence).as_list()
            if operations and operations[0][0] == 'delete' and operations[0][1] <= pos:
                i -= 1
                break
            pos = operations[0][1] if operations else None

        self.align_start += i


    def _trim_excess_insertions(self):
        # trim excess inserts from end of mapping, these are mapping random noise because of deletions
        operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, self.aligned_read_sequence).as_list()
        pos = operations[-1][1] if operations else None
        i = 0
        while operations and operations[-1][0] == 'insert' and operations[-1][1] == pos:
            i += 1
            aligned_read_sequence = self.full_read_sequence[self.align_start:self.align_end-i]
            operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, aligned_read_sequence).as_list()
            if operations and operations[-1][0] == 'delete':
                i -= 1
                break

        self.align_end -= i


    def _trim_excess_deletions(self):
        # trim excess deletions from end of mapping, these are mapping overhangs from the reference
        ref_length = self.ref_end-self.ref_start
        operations = rapidfuzz.distance.Levenshtein.editops(self.ref_sequence, self.aligned_read_sequence).as_list()
        n_operations = len(operations)

        i = 0
        while operations and operations[-1][0] == 'delete' and (ref_length-i) >= self.read_length:
            i += 1
            ref_sequence = self.ref_sequence[self.ref_start:self.ref_end-i]
            operations = rapidfuzz.distance.Levenshtein.editops(ref_sequence, self.aligned_read_sequence).as_list()
            if operations and len(operations) > (n_operations-i):
                i -= 1
                break

        self.ref_end -= i


    def _close_deletion_gaps(self):
        # remove gaps between continuous deletions, these may be caused by wrong assignment of repeat bases

        # only relevant if there are at least two delevents
        delevents = []
        for i, op in enumerate(self.operations):
            if op[0] == 'delete':
                delevents.append(i)
        if len(delevents) < 2:
            return
        
        # go through all pairs of consecutive delevents
        for i in range(1, len(delevents)):
            event1 = self.operations[delevents[i-1]]
            event2 = self.operations[delevents[i]]

            # check distance between delevents, skip if gap is large
            dist = event2[1] - event1[2]
            if dist > 3:
                continue

            # bases that make up the gap            
            gap_bases = self.ref_sequence[event1[2]:event2[1]]

            # check end of second delevent
            del_bases = self.ref_sequence[event2[2]-len(gap_bases):event2[2]]

            # if gap and deleted bases are identical, perform shift
            if del_bases == gap_bases:
                shift = event2[2] - event1[2] - len(gap_bases)
                self.operations[delevents[i-1]] = ('delete', event1[1], event1[2]+shift, *event1[3:])
                for j in range(delevents[i-1]+1, delevents[i]):
                    event = self.operations[j]
                    self.operations[j] = (event[0], event[1]+shift, event[2]+shift, event[3], event[4])
                del self.operations[delevents[i]]

                # recurse, working on newly shifted operations
                return self._close_deletion_gaps()


    @property
    def ref_sequence(self):
        return str(self.ref.seq[self.ref_start:self.ref_end])

    @property
    def aligned_read(self):
        return self.read[self.align_start:self.align_end]

    @property
    def aligned_read_sequence(self):
        return self.full_read_sequence[self.align_start:self.align_end] 

    @property
    def distance(self):
        return rapidfuzz.distance.Levenshtein.distance(self.ref_sequence, self.aligned_read_sequence)

    @property
    def similarity(self):
        return rapidfuzz.distance.Levenshtein.normalized_similarity(self.aligned_read_sequence, self.ref_sequence, processor=str)


    @property
    def read2ref_assignment(self):
        read2ref_assignment = {i: None for i in range(len(self.aligned_read_sequence))}
        for operation in self.operations:
            if operation[0] == 'equal' or  operation[0] == 'replace':
                for i_read, i_ref in zip(range(operation[3], operation[4]), range(operation[1], operation[2])):
                    read2ref_assignment[i_read] = i_ref
        return read2ref_assignment


    @property
    def read2type_assignment(self):
        read2type_assignment = {i: None for i in range(len(self.aligned_read_sequence))}
        for operation in self.operations:
            if operation[0] == 'equal': 
                for i_read in range(operation[3], operation[4]):
                    read2type_assignment[i_read] = 'E'
            elif operation[0] == 'replace':
                for i_read in range(operation[3], operation[4]):
                    read2type_assignment[i_read] = 'S'
            elif operation[0] == 'insert':
                for i_read in range(operation[3], operation[4]):
                    read2type_assignment[i_read] = 'I'
            # deletions do not have a read assignment
        return read2type_assignment




    def __repr__(self):
        return f"Mapping(d={self.distance}, sim={self.similarity:.2f}, ref_id={self.ref.id}, ref_len={len(self.ref_sequence)}, read_len={len(self.aligned_read_sequence)}, op={self.operations})"




    def to_comparison_string(self):
        info = f">read_id:{self.read.id}|ref_id:{self.ref.id}|ref_len:{len(self.ref_sequence)}|read_len:{len(self.aligned_read_sequence)}|dist:{self.distance}|sim:{self.similarity:.2f}"
        ref_string = ""
        read_string = ""

        if self.align_start != 0:
            ref_string = " "*self.align_start + "|"
            read_string = self.full_read_sequence[0:self.align_start] + "|"

        for operation in self.operations:

            if operation[0] == 'equal':
                ref_string += self.ref_sequence[operation[1]:operation[2]]
                read_string += self.aligned_read_sequence[operation[3]:operation[4]]

            elif operation[0] == 'delete':
                op_len = operation[2]-operation[1]
                ref_string += self.ref_sequence[operation[1]:operation[2]]
                read_string += "-"*op_len

            elif operation[0] == 'insert':
                op_len = operation[4]-operation[3]
                read_string += self.aligned_read_sequence[operation[3]:operation[4]]
                ref_string += "-"*op_len   
            
            elif operation[0] == 'replace':
                op_len = operation[4]-operation[3]
                ref_string += self.ref_sequence[operation[1]:operation[2]].lower()
                read_string += self.aligned_read_sequence[operation[3]:operation[4]].lower()

        remaining_ref = self.ref_sequence[operation[2]:]
        if remaining_ref:
            ref_string += "|" + remaining_ref
        remaining_read = self.full_read_sequence[(self.align_start + operation[4]):]
        if remaining_read:
            read_string += "|" + remaining_read

        return "\n".join([info, ref_string, read_string, "\n"])
