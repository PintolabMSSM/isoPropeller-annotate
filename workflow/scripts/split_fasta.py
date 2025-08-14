from pathlib import Path
from itertools import islice
from Bio import SeqIO

outdir = Path(snakemake.output[0])
outdir.mkdir(parents=True, exist_ok=True)

seqs_per_chunk = int(snakemake.params.seqs_per_chunk)
records = SeqIO.parse(snakemake.input.fasta, "fasta")

chunk_ids = []
chunk_num = 0
while True:
    chunk = list(islice(records, seqs_per_chunk))
    if not chunk:
        break
    cid = f"{chunk_num:05d}"
    with (outdir / f"chunk_{cid}.fasta").open("w") as fh:
        SeqIO.write(chunk, fh, "fasta")
    chunk_ids.append(cid)
    chunk_num += 1

with (outdir / "chunk_ids.txt").open("w") as fh:
    fh.write("\n".join(chunk_ids))
