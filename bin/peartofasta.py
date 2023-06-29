import sys
print(f"Reading in {sys.argv[1]}")
with open(f'{sys.argv[1]}/pear.all.good') as fi, open(f'{sys.argv[1]}/pear.all.good.fasta', 'w') as fo:
    for i, seq in enumerate(fi.readlines()):
        fo.writelines(f">Seq{str(i).zfill(9)}\n{seq}\n")
