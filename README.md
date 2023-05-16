# morphOMICs_v2

`morphOMICs_v2` is a Python package containing tools for analyzing microglia morphology using a topological data analysis approach. Note that this algorithm is designed not only for microglia applications but also for any dynamic branching structures across natural sciences.

- [Overview](#overview)
- [Required Dependencies](#required-dependencies)
- [Installation Guide](#installation-guide)
- [Usage](#usage)

# Overview
`morphOMICs_v2` is a topological data analysis approach which combines the Topological Morphology Descriptor (TMD) with bootstrapping approach, dimensionality reduction strategies to visualize microglial morphological signatures and their relationships across different biological conditions.


# Required Dependencies
Python : 3.8

numpy : 1.8.1+, scipy : 0.13.3+, pickle : 4.0+, enum34 : 1.0.4+, scikit-learn : 0.19.1+, matplotlib : 3.2.0+

Additional dependencies:
umap-learn : 0.3.10+, fa2 (https://github.com/bhargavchippada/forceatlas2)

# Installation Guide
Because morphOMICs relies on forceatlas2, which has problems with python3.9+, we will need to create a virtual environment for morphOMICs:

```
conda create -n morphology python=3.8 scipy numpy pandas scikit-learn matplotlib tomli networkx Cython
conda activate morphology
pip install fa2

git clone https://git.ista.ac.at/rcubero/morphomics_v2.git
cd morphomics_v2
python3 setup.py install
```

# Usage
To run a typical morphOMICs pipeline, create a parameter file with filename Morphomics.Parameters.[Parameters ID].toml (see examples).
The parameter file is build such that it modularizes the steps required to generate the phenotypic spectrum.
Once you have completed filling up the necessary information in the parameter file, run 
`python3 run_morphomics.py [path-to-parameter-file]`

To get started, download a demo folder containing sample traces and a parameter file from:
https://seafile.ist.ac.at/f/7dc51d9b63d24fb28eba/?dl=1
Once downloaded, place the extracted folder `morphomics_v2_demo` inside `morphomics_v2`.
On the terminal, run
`python3 run_morphomics.py morphomics_v2_demo/parameters/Morphomics.Parameters.1234.toml`

<!-- 
The easiest way to navigate through `MorphOMICs` is to run the `Morphomics_demo notebook`:
  - download `demo.zip` from https://seafile.ist.ac.at/f/eb13e707041749269ff9/?dl=1
  - unzip `demo.zip`
  - `cd demo`
  - `jupyter notebook`
  - Copy the url it generates, it looks something like this: `http://127.0.0.1:8888/?token=a4d016c37e162499e17b2993e69073fac0018bd9a779b762`
  - Open it in your browser
  - Then open `Morphomics_demo.ipynb` -->
