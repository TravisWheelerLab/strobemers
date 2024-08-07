import subprocess

files = ["ERR3588905_1", "ERR3588903_1"]
subdir = "cDNA"
for f in files:
    subprocess.run(
        ["pychopper",
        f"{subdir}/{f}.fastq",
        f"{subdir}/{f}_full.fastq",
        "-S", f"{subdir}/{f}_stats.tsv",
        "-r", f"{subdir}/{f}_pychopper.pdf"]
    ) 