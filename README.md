# d3v-gsd-occ
Extension of Design visualizer for general ship design that uses python-occ and PyGEM
## Recomended instalation procedure
### Prepare environment
Run miniconda or anaconda command prompt/power shell  
conda create -n d3vgsdocc_env python=3.9 pythonocc-core -c conda-forge  
conda activate d3vgsdocc_env   
pip install pyside6  
pip install openmesh  
pip install pyopengl  
pip install scipy
### Install PyGEM
Folow the instalation procedure given in:
https://github.com/mathLab/PyGeM/blob/master/README.md#dependencies-and-installation
### Prepare d3v-gsd
Download or clone d3v: 
https://github.com/linaetal-fsb/d3v.git  
(or use some specific fork)  
Download or clone d3v-gsd: 
https://github.com/pprebeg/d3v-gsd
(currently work using branch )
Download or clone d3v-gsd-occ: 
https://github.com/pprebeg/d3v-gsd-occ

### How to run
Run as: Path\to\d3v\src\main.py -a Path\to\d3v-gsd -a Path\to\d3v-gsd-occ 
