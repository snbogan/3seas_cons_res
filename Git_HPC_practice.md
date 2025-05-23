Github and HPC basics
================

# Part 1 - Using the Explorer Cluster

### Logging into Explorer (NU’s instructional HPC)

Open the terminal and enter the following command

``` bash

### For Mac
# Log in using ssh
ssh [username]@login.explorer.northeastern.edu
# type yes and hit enter if prompted

# Type and enter your password (no text will show while this happens)
```

For PC users, try the command above. If it does not work, use the
computer lab iMac in the meantime. Afterward, try downloading the
[Windows Subsystem Linux
(WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).

Create a working directory for yourself

``` bash
# Make new directory for our class
mkdir crms
```

Look at the SLURM queue to see what jobs are running

``` bash
# Peak at queue
squeue
```

### Load software modules

There are two ways to use software on Explorer: loading and installing

Look at what software is preinstalled on Explorer and can be loaded

``` bash
# Peak at software modules
module avail

# Load some software as a test
module load vcftools

# Check that vcftools succesfully loaded
vcftools --help
```

### Install software using a Conda environment

Load conda, create a conda environment, activate, and install software.
We’ll install vcftools, which we will use in later modules. Vcftools is
a software package for analyzing genotype data stored in variant calling
format (vcf) files.

``` bash
# From any directory, load conda
module load miniconda3

# Navigate to your home directory
cd /path/to/home/directory

# Create conda environment
conda create --name vcftools # Respond yes to all prompts
```

Close terminal and then log back in.

``` bash
# Activate the new environment
conda activate vcftools

# Install vcf tools
conda install bioconda::vcftools # Respond yes to all prompts

# Deactivate environment after installation
conda deactivate
```

### Write and submit a SLURM job

Finish the job script below by completing the missing SBATCH headers.

``` bash

#!/bin/bash
#SBATCH --job-name=[name] # Give your job a name
#SBATCH --output=[name].out # Name a file to contain code outputs
#SBATCH --error=[name].err # Name a file to contain error messages
#SBATCH --time=0-01:00:00 # Edit this parameter to specify a time
#SBATCH --mail-user=[email] # Add your email
#SBATCH --ntasks=1 # Let's start with one task
#SBATCH --ntasks-per-node=1 # Let's start with just one task per node
#SBATCH --mail-type=ALL
#SBATCH --mem=[number]G # How many gigabytes of memory does this job need?
#SBATCH --cpus-per-task= [number] # How many CPUs do you think this code needs?

# Make a new directory
mkdir /home/[username]/crms/day2_test

# Navigate to that new directory
cd /home/[username]/crms/day2_test

# Print a messge in a new file to confirm that this job ran
echo "message of your choosing" > day2_test.txt
```

Create a new directory for code and store your script there

``` bash

# In your crms directory, enter
mkdir code

# Navigate to new directory
cd code

# Then, create a new .sh file for the job above
nano day2_test.sh
```

Paste your edited script into the open nano editor.

After pasting, press Ctrl+O to write (save) the file, then hit Enter to
confirm the file name.

Press Ctrl+X to exit nano.

### Run your job script

Evaluate and queue your job to run using sbatch

``` bash
# Take a peak at your new script
cat day2_test.sh

# If everything looks good, queue the job to run
sbatch day2_test.sh

# Check if your job succesfully ran!
cat /home/[username]/crms/day2_test/day2_test.txt
```

# Part 2 - Using git

Before beginning part 2 make sure you:

1.  Have a GitHub account
2.  Have R Studio installed
3.  Install Filezilla client (software for moving files between local
    and remote directories)

[Install Filezilla client here](https://filezilla-project.org/)

## Git through R Studio

### Getting started

1.  Make a GitHub account
2.  Set up login key
3.  Log into your GitHub account in your terminal

### Create an R Project that is version controlled with git

1.  Create a new repository on Github called ‘crms_test’ and include a
    README
2.  Copy the url to the repo
3.  Open R Studio \> File \> New Project \> Version Controlled and enter
    repository name, url, and local directory
4.  Edit the README to say ‘Hello world!’

### Practice git commands in R Studio environment

1.  Open the R Project file for CRMS_test
2.  Hit the pull icon (down arrow)
3.  Stage the README.md file
4.  Hit the commit button (check mark), which will queue updates to push
5.  Write a commit message. For example, “Add text to README”
6.  Hit ‘Push’ to update the repository on GitHub
7.  Open the webpage for your Github repository

## Git on the terminal

1.  Edit your README.md file again in R Studio to add text of your
    choosing
2.  Open Terminal
3.  Navigate to the repo then push your change using the code below:

``` bash
# Navigate to your repo
cd /path/to/my/repo

# Always pull before you push!
git pull

# Stage you README
git add README.md

# Commit the change
git commit -m "Push change from terminal"

# Push the change
git push
```
