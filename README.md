<p align="center"><img src="logo.png" align="left" alt="MOSFiT" width="300"/></p>
<a href="https://travis-ci.org/guillochon/MOSFiT"><img src="https://img.shields.io/travis/guillochon/MOSFiT.svg?branch=master" alt="Build Status"></a>
<a href="https://coveralls.io/github/guillochon/MOSFiT?branch=master"><img src="https://coveralls.io/repos/github/guillochon/MOSFiT/badge.svg?branch=master" alt="Coverage Status"></a>
<a href="https://www.python.org"><img src="https://img.shields.io/badge/python-2.7%2C%203.4%2C%203.5%2C%203.6-blue.svg" alt="Python Version"></a>
<a href="https://badge.fury.io/py/mosfit"><img src="https://badge.fury.io/py/mosfit.svg" alt="PyPI version"></a>
<a href="http://mosfit.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/mosfit/badge/?version=latest" alt="Documentation Status"></a>
<a href="http://ascl.net/1710.006"><img src="https://img.shields.io/badge/ascl-1710.006-blue.svg?colorB=262255" alt="ascl:1710.006" /></a>
<a href="https://slack.astrocats.space"><img src="https://slack.astrocats.space/badge.svg" alt="Currently logged-in users in MOSFiT Slack channel" /></a>

`MOSFiT` (**M**odular **O**pen-**S**ource **Fi**tter for **T**ransients) is a Python 2.7/3.x package for fitting, sharing, and estimating the parameters of transients via user-contributed transient models. Data for a transient can either be provided by the user in a wide range of formats (JSON, ASCII tables, CDS, LaTeX), or can be pulled automatically from one of the Open Catalogs (e.g. the [Open Supernova Catalog](https://sne.space), [Open TDE Catalog](https://tde.space), and [Open Nova Catalog](https://opennova.space)) by its name. With the use of an optional upload flag, fits performed by users can then be uploaded back to the Open Catalogs for the benefit of the transient community.<br clear="all">

This version of MOSFiT was modified by Kamile Lukosiute to include a kilonova model based on [D. Kasen's Kilonova SEDs](https://github.com/dnkasen/Kasen_Kilonova_Models_2017). 

## Installation
First, you will need to set up a new conda environment (please just trust me and make a new one) using the .yml file included (edit to name your conda environment what you would like it to be). 

Then, you will need to download the reformatted SED pickles found [here](https://github.com/SSantosLab/kasen_seds) and put the files (they are big) somewhere that makes sense. Next,in kasen0.py and kasen1.py, you will need to change `self._dir_path` to point to the directory containing the kasen_seds directory. If you are using the DES machines, put them on the big disks. 

Finally, you will be able to install MOSFiT. Use the setup.py script:
```bash
git clone https://github.com/guillochon/MOSFiT.git
cd MOSFiT
python setup.py develop
```

## Using MOSFiT

For detailed instructions on using MOSFiT, please see our documentation on RTD: <http://mosfit.readthedocs.io/>

For more details about the changes Kamile has made, ask Kamile. 
