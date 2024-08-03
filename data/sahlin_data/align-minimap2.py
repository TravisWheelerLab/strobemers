import subprocess

reference_file = "references/SIRV_isoforms_multi-fasta_170612a.fasta"
query_file = "cDNA/ERR3588905_1_full.fastq"
output_file = "aligned_reads.sam"

cmd = [
    "minimap2", 
    "-ax", "splice",
    reference_file,
    query_file
]

# Running the command and capturing the output
with open(output, "w") as out_file:
    subprocess.run(cmd, stdout=out_file)

