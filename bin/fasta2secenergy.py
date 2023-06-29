import sys
import RNA
import Bio.SeqIO 
import tqdm

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]


fasta_sequences = Bio.SeqIO.parse(open(INPUT_FILE), 'fasta')

with open(INPUT_FILE) as fi, open(OUTPUT_FILE, "w") as fo:
    for record in tqdm.tqdm(Bio.SeqIO.parse(fi, 'fasta')):

        # compute minimum free energy (MFE) 
        _, mfe_fw = RNA.fold(str(record.seq))
        _, mfe_rv = RNA.fold(str(record.seq.reverse_complement()))

        fo.write(",".join([record.id, str(mfe_fw), str(mfe_rv)]))
        fo.write("\n")