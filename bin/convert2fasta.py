import sys

from_file = sys.argv[1]
to_file = sys.argv[2]

sequences = []

with open(from_file) as f:
    for seq in f.readlines():
        sequences.append(seq)

with open(to_file, "w") as f:
    for i, seq in enumerate(sequences):
        txt = ">" + str(i).zfill(6) + "\n" + seq.strip("\n") + "\n\n"
        f.write(txt)
