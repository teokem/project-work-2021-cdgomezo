<a href="https://zenodo.org/badge/latestdoi/351733532"><img src="https://zenodo.org/badge/351733532.svg" alt="DOI"></a>
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/teokem/project-work-2021-cdgomezo/blob/main/Jupyter%20Course%20Project.ipynb)


## Inverse modeling for calculating the land carbon cycle
Inverse modeling is a commonly used method and a formal approach to estimate the variables driving the evolution of a system, e.g. greenhouse gases (GHG) sources and sinks, based on the observable manifestations of that system, e.g. GHG concentrations in the atmosphere. This has been developed and applied for decades and it covers a wide range of techniques and mathematical approaches as well as topics in the field of the biogeochemistry. In this Jupyter Notebook is a lecture for students interested in learning about inverse modeling. It contains all the theoretical background around the concept of inverse modeling from a beginners level. At the end of the notebook you can find an application of inverse modeling to retrieve sourface fluxes of carbon dioxide.

## What you will find in this GitHub repository

The purpose of this GitHub repository is to host the required files to run a Jupyter notebook. You will find the file ``` environment.yml ```, which sets the enviroment needed to run the Jupyter notebook 

## Instructions for running the notebook

1. Install [miniconda3](https://docs.conda.io/en/latest/miniconda.html).

2. Download the files from this repository and unzip

3. In the terminal, navigate to the folder you downloaded from GitHub

4. Install the ``` Inverter ``` environment by running the follwing lines  
	
  ```
  conda env create -f environment.yml
  conda activate Inverter	  
  ```
5. Run the notebook by typing
```
jupyter notebook
```

## Environment packages
<div style="text-align: justify">
<ul>
  <li><code>numpy</code> and its submodule <code>numpy.linalg</code> will be the main library to read, format and manipulate datasets as well as perform mathematic and linear algebra operations. </li>
  <li><code>matplotlib</code> will be the main library to plot the different results and datasets.</li>
  <li><code>copy</code> is used in the <code>inverter</code> package to copy the format of the datasets.</li>
  <li><code>matplotlib.dates</code>, <code>datetime</code> and <code>calendar</code> are used to format time parameter of fluxes and observations.</li>
  <li><code>inverter</code> is a Python class containing the functions to perform the inversions.</li>
  <li><code>modelCO2</code> contains the Python code of the transport model.</li>
</ul> 
</div>
