
# <center>MIRISim TSO</center>
<center>Post-Treatment of MIRISim Images for Time-Series Observations</center>


MRISim_TSO is a Python 3 package. As MIRISim has not been designed for Time-Series Observations (TSO), in order to create those, one has to start a new simulations at each time step. This implies getting rig of all the systematic noises and drifts. The present package aims at reintroducing the detector systematics in the files resulting of a MIRISim computation.

The package have been designed exclusively to treat MIRISim output files, hence we try to be up to date on the data model. Please make sure you comply to the format and folder arborescence.

Here, we create new files in the MIRISim result folder containing a new `det_image` accounting for instrumental systematics.

## I . Getting Started
> _Important Note_ : Current version is still under development.

These instructions will get you to have a version of the package on your computer up and running. Note that a complete and detailed user manual is available in the documentation. It contains some important usage rules, and we strongly advise a careful reading before using this treatment.

Every step described hereafter have to be done in the directory containing the `setup.py` file (usually, the downloaded folder).

###  1 . Requirements
The package is developed in `Python 3.6`. The required python libraries are specified in the file 'requirements.txt'. To ensure all dependencies are satisfied, run the following command line:
```bas
pip install -r requirements.txt
```


### 2 . Installation
You can install mirisim_tso as a standard Python module. To do that, in the root directory of the package, do:
```bash
python setup.py install
```
In your python script or in a python environment, import the module with:
```python
import exonoodle
```

For developers, create a symbolic link so that each time the code is changed, the imported module is changed accordingly (still in the root directory):
```bash
pip install -e .
```

### 3 . Run in native script



## II . Author
   - ***Marine Martin-Lagarde*** *(Corresponding author)* - CEA-Saclay - marine.martin-lagarde@cea.fr
   - ***Christophe Cossou*** - IAS-Orsay - *Python support and package architecture*

## III . Licence
**<center>Work in progress</center>**

## IV . Acknowledgments
If you want to use this code in a scientific publication, it would be appreciated if you cite us. *No referenced article yet*   
The author is partly funded by a CNES grant. The research leading to these development has received funding from the European Union’s Horizon 2020 Research and Innovation Programme, under Grant Agreement 776403.
**<center>Work in progress</center>**

## V . References
The reference articles for the calculations used in the code are the following :   
- D.Dicken et al. _in prep_
