# my_functions
Python personal functions for solving FEM and localisation

#### Download miniconda
https://docs.conda.io/en/latest/miniconda.html#linux-installers

#### Install miniconda
```console
chmod +x Miniconda...
./Miniconda...
```

(Use the TAB key to autocomplete...)

#### Create a virtual environement named 'fb' with spdyer (recomended for debugg)
```console
conda create -n fb spyder
```

#### Or create a virtual environement named 'fb' 
```console
conda create -n fb 
```

#### Install mamba on the base environement
```console
conda install mamba -n base -c conda-forge
```

#### Activate the environement
```console
conda activate fb
```

#### Install fenics with mamba
```console
mamba install -c conda-forge fenics 
```

#### Install some packages with mamba
```console
mamba install -c conda-forge scikit-image matplotlib pandas pytables=3.6.1 numba pycairo
```

#### Install some packages with pip
```console
pip install gmsh==4.10 meshio==4.3.8 mpld3==0.5.2 scipy numpy notebook 
```

#### If h5py is missing, try (inside the environement)
```console
git clone https://github.com/h5py/h5py.git
cd h5py
pip install -v .
```

#### Clone fenics-shells
```console
git clone https://github.com/charlie-garcia/my_functions.git
cd my_functions
python setup.py install
```

#### Clone bempp-cl
```console
git clone https://github.com/charlie-garcia/my_functions.git
cd my_functions
python setup.py install
```

#### Clone this repository and install it
```console
git clone https://github.com/charlie-garcia/my_functions.git
cd my_functions
python setup.py install
```

* Tested on jul 4 2020, conda -V, python -V 
