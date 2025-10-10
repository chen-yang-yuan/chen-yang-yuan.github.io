---
layout: page
title: HGCC Guidelines
---

# HGCC Beginner Guidelines

> This note is about the basic use of the Emory Human Genetics Compute Cluster (HGCC). It comes from a realistic case in my research, where I couldn’t read a large raw data file locally and needed to use the HPC cluster to execute my data cleaning script. For detailed HGCC guidelines, please check [Dr. Yang's notes 1](https://yanglab-emory.github.io/assets/ComputationSlides/HGCC_StartGuide_1.html) and [2](https://yanglab-emory.github.io/assets/ComputationSlides/HGCC_StartGuide_2.html).

Here is my local project folder structure:

~/Dropbox/SG/Xenium_5K_OC_FF/ \
|---- raw_data/ \
|---- data/ \
|---- output/ \
|---- 1_clean.py

My `1_clean.py` looks like this:

```python
import numpy as np
import pandas as pd

# Read transcripts
transcripts = pd.read_parquet("raw_data/transcripts.parquet", columns = ["cell_id", "is_gene", "overlaps_nucleus", "feature_name", "x_location", "y_location", "z_location"])
transcripts = transcripts[transcripts["is_gene"] == True]

# Process transcripts
transcripts = transcripts[["cell_id", "overlaps_nucleus", "feature_name", "x_location", "y_location", "z_location"]]
transcripts = transcripts.rename(columns = {"overlaps_nucleus": "in_nucleus", "feature_name": "target", "x_location": "global_x", "y_location": "global_y", "z_location": "global_z"})

# Relative position to cell and nucleus
transcripts["in_cell"] = (transcripts["cell_id"] != "UNASSIGNED").astype(int)
transcripts["overlaps_nucleus"] = (transcripts["in_nucleus"] == transcripts["in_cell"]).astype(int)

print(f"In-nucleus ratio: {100 * np.sum(transcripts['in_nucleus'] == 1) / transcripts.shape[0]:.2f}%")
print(f"In-cytoplasm ratio: {100 * np.sum(transcripts['overlaps_nucleus'] == 0) / transcripts.shape[0]:.2f}%")

# Save transcripts
transcripts = transcripts[["cell_id", "overlaps_nucleus", "target", "global_x", "global_y", "global_z"]]
transcripts.to_parquet("data/transcripts.parquet")

# In-cytoplasm ratio
gene_means = transcripts.groupby("target")["overlaps_nucleus"].mean().reset_index()
gene_means.columns = ["gene", "in_nucleus_ratio"]
gene_means = gene_means.sort_values(by = "in_nucleus_ratio", ascending = True)
gene_means["in_cytoplasm_ratio"] = 1 - gene_means["in_nucleus_ratio"]
gene_means.to_csv("output/in_cytoplasm_ratio.csv", index = 0)
```

Essentially, this script reads the raw `transcript.parquet` from the `raw_data/` folder, performs essential cleaning, saves the cleaned data to the `data/` folder, and outputs summary statistics to the `output/` folder.

Below are the complete steps I use to run the script on HGCC.

## 1. First login

First, make sure you connect to [Emory VPN](https://it.emory.edu/security/vpn.html).

Login by running the bash command:

```bash
ssh cyuan36@hgcc.emory.edu
```

Input Emory password.

## 2. Transfer files to HGCC

Create the directory structure that mirrors the local version:

```bash
mkdir -p ~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/{raw_data,data,output}
```

Exit HGCC:

```bash
exit
```

Transfer the data and script from the local computer to HGCC by running:

```bash
rsync -avP ~/Dropbox/SG/Xenium_5K_OC_FF/raw_data/transcripts.parquet \
cyuan36@hgcc.emory.edu:~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/raw_data/

rsync -avP ~/Dropbox/SG/Xenium_5K_OC_FF/1_clean.py \
cyuan36@hgcc.emory.edu:~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/
```

You need to input Emory password in the process.

## 3. Install necessary packages

My script contains two Python packages that have to be installed beforehand: `numpy` and `pandas`.

First, login and enter the interactive node:

```bash
ssh cyuan36@hgcc.emory.edu
srun -N 1 -n 1 --pty --preserve-env bash
```

Figure out the available Python module:

```bash
module available
```

There is a line indicating the Python version:

```
python/3.11.7-gcc-13.2.0-ljcvqdx
```

Next, load the module by running:

```bash
module load python/3.11.7-gcc-13.2.0-ljcvqdx
```

Confirm it worked:

```bash
python3 --version
```

You should see something like:

```
Python 3.11.7
```

Finally, you are ready to install the packages for *that* Python:

```bash
pip install --user numpy pandas
```

Exit the interactive node and HGCC:

```bash
exit
exit
```

## 4. Create SLURM job script

On the local computer, create a new text file named `1_clean.sh` in the same project folder:

```bash
#!/bin/bash
#SBATCH --job-name=clean_xenium
#SBATCH --output=output/clean_xenium_%j.out
#SBATCH --error=output/clean_xenium_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# Load Python 3.11 module (same one you used interactively)
module load python/3.11.7-gcc-13.2.0-ljcvqdx

# Print some debugging info
echo "Running on $(hostname)"
python3 --version
python3 -m site --user-site

# Move to your project directory
cd ~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF

# Run your script
python3 1_clean.py

echo "Job finished at $(date)"
```

To make sure the partition name valid, run:

```bash
ssh cyuan36@hgcc.emory.edu
sinfo
```

It will output something like:

```
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
nodes*    up     30-00:00:0     2  drain node[02,04]
nodes*    up     30-00:00:0     5  mix   node[01,03,05,09,12]
nodes*    up     30-00:00:0     1  alloc node06
```

The only partition on the cluster is called `nodes` (the asterisk `*` means it’s also the default partition), so the partition name should be it, or simply remove that line entirely.

Also, since the raw data is very large, I set ```--mem=1000G``` to avoid the system running out of memory (OOM). The memory limit can be checked by:

```bash
sinfo -o "%P %m %G"
```

The output is:

```
PARTITION MEMORY GRES
nodes* 1024000 (null)
```

The total memory I can request is `1024000` MB, so setting the parameter to `1000G` is fine (or also `0G` for no upper limit).

Exit the interactive node and HGCC:

```bash
exit
exit
```

Transfer the job script to HGCC in a similar way:

```bash
rsync -avP ~/Dropbox/SG/Xenium_5K_OC_FF/1_clean.sh \
cyuan36@hgcc.emory.edu:~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/
```

## 5. Submit the job on HGCC

Login to HGCC:

```bash
ssh cyuan36@hgcc.emory.edu
cd ~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF
```

Before submitting, confirm the script is executable. This step shouldn't print any output:

```bash
chmod +x 1_clean.sh
```

Alternatively, you can confirm it worked by checking:

```bash
ls -l 1_clean.sh
```

You should see something like:

```
-rwxr-xr-x 1 cyuan36 hulab  512 Oct  8 21:42 1_clean.sh
```

Then submit the job:

```bash
sbatch 1_clean.sh
```

You will see something like:

```
Submitted batch job 487120
```

where `487120` is the job ID. Check your job status:

```bash
squeue -u cyuan36
```

Once it's done, you will find:
* **Standard output:** `output/clean_xenium_487120.out`
* **Error logs:** `output/clean_xenium_487120.err`
* **Generated files:** in `data/` and `output/`

You can preview logs:

```bash
cat output/clean_xenium_487120.out
cat output/clean_xenium_487120.err
```

Exit HGCC:

```bash
exit
```

## 6. Retrieve output back to local computer

From the local terminal, run:

```bash
rsync -avP cyuan36@hgcc.emory.edu:~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/output/ \
~/Dropbox/SG/Xenium_5K_OC_FF/output/

rsync -avP cyuan36@hgcc.emory.edu:~/hulab/projects/SG/Xenium_5K_Ovarian_Cancer_FF/data/ \
~/Dropbox/SG/Xenium_5K_OC_FF/data/
```