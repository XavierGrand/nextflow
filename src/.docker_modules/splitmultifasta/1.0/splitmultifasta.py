import sys
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    (derived from https://biopython.org/wiki/Split_large_file)
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                break # EOF = end of file
            batch.append(entry)
        if batch:
            yield batch

if(len(sys.argv) != 3):
        sys.exit("usage: split_fasta.py MULTI_FASTA_FILE OUTPUT_Path")

ffile=sys.argv[1]  # fasta file
output=sys.argv[2] # output path

for record in SeqIO.parse(open(ffile), "fasta"):
    seqName = record.id.split("|" or " ")
    record.id=seqName[2]
    fastaFileName = output + record.id + ".fasta"
    SeqIO.write(record, fastaFileName, "fasta")

sys.exit(0)
