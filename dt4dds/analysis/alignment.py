import multiprocessing
import collections
import rapidfuzz
import edlib
import functools
import numpy as np
import Bio.Seq
import Bio.SeqRecord
import logging

from . import mapping


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


alignment = collections.namedtuple("alignment", ['read', 'ref', 'similarity', 'data'])







def align_single_reads_init(dict_fw, dict_rv, readwindow, scorecutoff, trimseq, readoffset, refoffset, readlength):

    global read_window, score_cutoff, trim_seq, references_dict_fw, references_dict_rv, read_offset, ref_offset, read_length
    read_window = readwindow
    score_cutoff = scorecutoff
    trim_seq = trimseq
    read_offset = readoffset
    ref_offset = refoffset
    read_length = readlength
    references_dict_fw = dict_fw
    references_dict_rv = dict_rv

    # create dicts of smaller ref sequences according to the read window used for matching
    global windowed_ref_dict_fw
    if read_window:
        windowed_ref_dict_fw = {ref_id: str(ref.seq[0+ref_offset:read_window+ref_offset]) for ref_id, ref in references_dict_fw.items()}
    else:
        windowed_ref_dict_fw = {ref_id: str(ref.seq) for ref_id, ref in references_dict_fw.items()}



def align_single_reads(readobject):

    # unpack the forward and reverse reads
    read = str(readobject.seq)

    # if requested, shorten the reads to the desired window length
    if read_window:
        short_read = read[0+read_offset:read_window+read_offset]
    else:
        short_read = read

    # extract the best-matching reference
    bestMatch = rapidfuzz.process.extractOne(
        short_read,
        windowed_ref_dict_fw,
        scorer=rapidfuzz.distance.Levenshtein.normalized_similarity,
        score_cutoff=score_cutoff,
        processor=None
    )

    # compile forward and reverse alignments
    readmapping = mapping.Mapping(readobject, references_dict_fw[bestMatch[2]], read_number=1, trim_seq=trim_seq, read_length=read_length)

    return readmapping




def align_single_reads_caller(reads, references_dict_fw, references_dict_rv, read_window=None, score_cutoff=0, trim_seq=None, read_offset=0, ref_offset=0, read_length=None):

    logger.info(f"Aligning {len(reads)} single reads with read windows {read_window} and score cutoff {score_cutoff}.")

    
    with multiprocessing.Pool(
        4, 
        initializer=align_single_reads_init, 
        initargs=(references_dict_fw, references_dict_rv, read_window, score_cutoff, trim_seq, read_offset, ref_offset, read_length),
    ) as pool:
        mappings = pool.imap(
            align_single_reads,
            reads,
            chunksize=min(10000, -(-len(reads) // 4)),
        )

        for mapping in mappings:
            yield (mapping,)










def align_paired_reads_init(dict_fw, dict_rv, readwindow, scorecutoff, trimseq, readoffset, refoffset, readlength):

    global read_window, score_cutoff, trim_seq, references_dict_fw, references_dict_rv, read_offset, ref_offset, read_length
    read_window = readwindow
    score_cutoff = scorecutoff
    trim_seq = trimseq
    read_offset = readoffset
    ref_offset = refoffset
    read_length = readlength
    references_dict_fw = dict_fw
    references_dict_rv = dict_rv

    # create dicts of smaller ref sequences according to the read window used for matching
    global windowed_ref_dict_fw
    if read_window:
        windowed_ref_dict_fw = {ref_id: str(ref.seq[0+ref_offset:read_window+ref_offset]) for ref_id, ref in references_dict_fw.items()}
    else:
        windowed_ref_dict_fw = {ref_id: str(ref.seq) for ref_id, ref in references_dict_fw.items()}




def align_paired_reads(reads):

    # unpack the forward and reverse reads
    readobject_fw, readobject_rv = reads
    read_fw = str(readobject_fw.seq)

    # if requested, shorten the reads to the desired window length
    if read_window:
        short_read = read_fw[0+read_offset:read_window+read_offset]
    else:
        short_read = read_fw

    # extract the best-matching reference
    bestMatch = rapidfuzz.process.extractOne(
        short_read,
        windowed_ref_dict_fw,
        scorer=rapidfuzz.distance.Levenshtein.normalized_similarity,
        score_cutoff=score_cutoff,
        processor=None
    )

    # compile forward and reverse alignments
    mapping_fw = mapping.Mapping(readobject_fw, references_dict_fw[bestMatch[2]], read_number=1, trim_seq=trim_seq, read_length=read_length)
    mapping_rv = mapping.Mapping(readobject_rv, references_dict_rv[bestMatch[2]], read_number=2, trim_seq=trim_seq, read_length=read_length)
    
    return mapping_fw, mapping_rv




def align_paired_reads_caller(reads, references_dict_fw, references_dict_rv, read_window=None, score_cutoff=0, trim_seq=None, read_offset=0, ref_offset=0, read_length=None):

    logger.info(f"Aligning {len(reads)} paired reads with read windows {read_window} and score cutoff {score_cutoff}.")

    
    with multiprocessing.Pool(
        4, 
        initializer=align_paired_reads_init, 
        initargs=(references_dict_fw, references_dict_rv, read_window, score_cutoff, trim_seq, read_offset, ref_offset, read_length),
    ) as pool:
        mappings = pool.imap(
            align_paired_reads,
            reads,
            chunksize=min(10000, -(-len(reads) // 4)),
        )

        for mapping_fw, mapping_rv in mappings:
            yield mapping_fw, mapping_rv







def categorize_mappings(mappings, min_match_similarity=0.7, medium_match_similarity=0.85, paired=False):
     
    logger.info(f"Categorizing reads in paired={paired} mode, using min. match similarity of {min_match_similarity}, and min. medium match similarity of {medium_match_similarity}.")
    
    matches = []
    low_scores = []
    medium_scores = []

    for mapping in mappings:

        if paired:
            
            # get the similarity score of read
            scores = np.array([m.similarity for m in mapping])

             # this is a low scoring match which might be a random simple collision
            if min(scores) < min_match_similarity:
                low_scores.append(mapping)
                continue

            # this is a medium scoring match which might be a corrupted read
            if min(scores) < medium_match_similarity:
                medium_scores.append(mapping)
                continue

            # this is a good match
            matches.append(mapping)

            
        else:
            
            for i_mapping in mapping:

                # get the similarity score of read
                score = i_mapping.similarity

                # this is a low scoring match which might be a random simple collision
                if score < min_match_similarity:
                    low_scores.append(i_mapping)
                    continue

                # this is a medium scoring match which might be a corrupted read
                if score < medium_match_similarity:
                    medium_scores.append(i_mapping)
                    continue

                # this is a good match
                matches.append(i_mapping)


    total_alignments = len(matches) + len(low_scores) + len(medium_scores)
    logger.info(f"Categorized {total_alignments} reads into {len(matches)} matches ({len(matches)/total_alignments*100:.2f}%), {len(medium_scores)} borderline matches ({len(medium_scores)/total_alignments*100:.2f}%), and {len(low_scores)} unknowns ({len(low_scores)/total_alignments*100:.2f}%).")
    return matches, medium_scores, low_scores




def compare_paired_mappings(mappings):

    paired_mappings = []

    for mapping_fw, mapping_rv in mappings:

        read_fw = str(mapping_fw.aligned_read_sequence)
        read_rv_rvcmp = str(Bio.Seq.Seq(mapping_rv.aligned_read_sequence).reverse_complement())

        operations = rapidfuzz.distance.Levenshtein.editops(read_fw, read_rv_rvcmp)

        ref_baselist_fw = list(read_fw)
        ref_baselist_rv = list(read_fw)

        fw_read2ref = mapping_fw.read2ref_assignment

        for op in operations:

            pos_fw = op.src_pos
            pos_rv = op.dest_pos

            if op.tag == "replace":

                mapped_pos_fw = fw_read2ref[pos_fw]
                if mapped_pos_fw is not None:
                    base_ref = mapping_fw.ref_sequence[mapped_pos_fw]
                    ref_baselist_fw[pos_fw] = base_ref
                    ref_baselist_rv[pos_fw] = base_ref
                else:
                    # there is no reference base, therefore assign as substitution for both reads
                    ref_baselist_fw[pos_fw] = read_rv_rvcmp[pos_rv]
                    ref_baselist_rv[pos_fw] = read_fw[pos_fw]


            elif op.tag == "delete":

                mapped_pos_fw = fw_read2ref[pos_fw]
                if mapped_pos_fw is not None:
                    # this was a deletion on reverse read
                    base_ref = mapping_fw.ref_sequence[mapped_pos_fw]
                    ref_baselist_fw[pos_fw] = base_ref
                    ref_baselist_rv[pos_fw] = base_ref
                else:
                    # this is an insertion on the forward read
                    ref_baselist_fw[pos_fw] = ''
                    ref_baselist_rv[pos_fw] = ''


            elif op.tag == "insert":

                rv_rvcmp_ref = mapping_rv.ref[:]
                rv_rvcmp_ref.seq = str(Bio.Seq.Seq(mapping_rv.ref_sequence).reverse_complement())
                read_rv_rvcmp_seqobject = Bio.SeqRecord.SeqRecord(id=mapping_rv.read.id, seq=read_rv_rvcmp)
                mapping_rv_revcomp = mapping.Mapping(read_rv_rvcmp_seqobject, rv_rvcmp_ref)

                mapped_pos_rv = mapping_rv_revcomp.read2ref_assignment[pos_rv] if pos_rv in mapping_rv_revcomp.read2ref_assignment else None
                if mapped_pos_rv is not None:
                    # this was a deletion on forward read
                    base_ref = mapping_rv_revcomp.ref_sequence[mapped_pos_rv]
                    ref_baselist_fw[pos_fw-1] += base_ref
                    ref_baselist_rv[pos_fw-1] += base_ref
                else:
                    # this is an insertion on the reverse read, we keep the forward read as reference
                    pass


        ref_sequence_fw = mapping_fw.ref[:]
        ref_sequence_fw.seq = ''.join(ref_baselist_fw)
        ref_sequence_rv = mapping_rv.ref[:]
        ref_sequence_rv.seq = str(Bio.Seq.Seq(''.join(ref_baselist_rv)).reverse_complement())


        map_fw = mapping.Mapping(mapping_fw.read, ref_sequence_fw, read_number=mapping_fw.read_number)
        map_rv = mapping.Mapping(mapping_rv.read, ref_sequence_rv, read_number=mapping_rv.read_number)
        paired_mappings.append([map_fw, map_rv])
    
    return paired_mappings














def align_phix_reads(ref_fw, ref_rv, reads, max_read_length=None, read_window=None):

    mappings = []

    for i, read in enumerate(reads):

        if max_read_length:
            read = read[0:max_read_length]

        read_seq = str(read.seq)

        # if requested, shorten the reads to the desired window length
        if read_window:
            short_read = read_seq[0:read_window]
        else:
            short_read = read_seq

        # try to align to each reference direction
        short_aln_fw = edlib.align(short_read, str(ref_fw.seq), mode = "HW", task = "locations")
        short_aln_rv = edlib.align(short_read, str(ref_rv.seq), mode = "HW", task = "locations")

        # choose best-fitting direction
        if short_aln_fw['editDistance'] < short_aln_rv['editDistance']:
            ref = ref_fw 
            short_aln = short_aln_fw
        else:
            ref = ref_rv
            short_aln = short_aln_rv

        # ensure that the alignment does not cross the split in the sequence
        shift = read_window // 2 + 1 if read_window else 20
        n = short_aln['locations'][-1][0] - shift
        ref = ref[n:] + ref[:n]

        # create full alignment
        errormap = mapping.Mapping.for_phix(read, ref, read_number=i+1, trim_seq="AGATCGGAAGAGCG")
        mappings.append(errormap)
        
    return mappings




def align_phix_reads_caller(reads, ref_fw, ref_rv, max_read_length=None, read_window=None):

    logger.info(f"Aligning {len(reads)} phix reads with max. read length {max_read_length}, and read window {read_window}.")

    with multiprocessing.Pool(4) as pool:
        mappings = pool.imap(
            functools.partial(
                align_phix_reads, 
                *[ref_fw, ref_rv], 
                **{'max_read_length': max_read_length, 'read_window': read_window}
            ),
            reads,
            chunksize=1000,
        )

        for mapping_fw, mapping_rv in mappings:
            yield mapping_fw, mapping_rv







def categorize_phix_mappings(mappings, min_match_similarity=0.70, medium_match_similarity=0.90):
     
    logger.info(f"Categorizing phix mappings using min. match similarity of {min_match_similarity}, and medium alignment similarity of {medium_match_similarity}.")
    
    matches = []
    low_scores = []
    medium_scores = []
    ambiguous_matches = []

    for mapping_fw, mapping_rv in mappings:

        similarity_alignments = np.array([match.similarity for match in (mapping_fw, mapping_rv)])
        match = min(similarity_alignments)

        # this is a low scoring match which might be a random simple collision
        if match < min_match_similarity:
            low_scores.extend((mapping_fw, mapping_rv))
            continue

        # this is a low scoring match which might be a random simple collision
        if match < medium_match_similarity:
            medium_scores.extend((mapping_fw, mapping_rv))
            continue

        # sort out reads where one read is badly matched, but the other is matchable
        if np.any(similarity_alignments < min_match_similarity):
            ambiguous_matches.extend((mapping_fw, mapping_rv))
            continue

        # this is a good match without a close second match candidate 
        matches.extend((mapping_fw, mapping_rv))


    total_alignments = len(matches) + len(low_scores) + len(medium_scores) + len(ambiguous_matches)
    logger.info(f"Categorized {total_alignments} reads into {len(matches)} matches ({len(matches)/total_alignments*100:.2f}%), {len(ambiguous_matches)} ambiguous matches ({len(ambiguous_matches)/total_alignments*100:.2f}%), {len(medium_scores)} borderline matches ({len(medium_scores)/total_alignments*100:.2f}%), and {len(low_scores)} unknowns ({len(low_scores)/total_alignments*100:.2f}%).")
    return matches, medium_scores, low_scores, ambiguous_matches





def seqmap_paired_mappings(mappings):
    mappings_by_refid = [{}, {}]
    for map in mappings:
        i = int(map.read_number)-1
        if map.ref.id not in mappings_by_refid[i].keys(): mappings_by_refid[i][map.ref.id] = []
        mappings_by_refid[i][map.ref.id].append(map)

    return mappings_by_refid


def seqmap_single_mappings(mappings):
    mappings_by_refid = [{}]
    for map in mappings:
        if map.ref.id not in mappings_by_refid[0].keys(): mappings_by_refid[0][map.ref.id] = []
        mappings_by_refid[0][map.ref.id].append(map)

    return mappings_by_refid