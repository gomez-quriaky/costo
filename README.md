# README COSTO

This is an algorithm that computes teleseismic joint inversion.

## Requriments

### Intel oneAPI  basic toolkit. 

[intel toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html#gs.jhm495)

The requarement to installed are 
 
- Linux OS
- Intel core procesor

Follow the intruction to download the toolkit.

Notice that the `ifx` compiler is installed within.


## Source Code

The project repository is held in *GitHub*. 
To obtain the project, do

        git clone https://github.com/gomez-quriaky/costo.git


To go to the source code

        cd costo

To change of branch (arch_test)

        git checkout arch_test

## Compilation


The project is build using CMAKE


    mkdir build
    cd build

To built the project, execute the following command:

    cmake ..
    make

It is posible to chose the compiler: `gfortran` or `ifx`.
To do so:

        cmake ..  -DCMAKE_Fortran_COMPILER=/path/to/the/compiler

## Running Costo

The inputs Files are in `fichier_input`. Therefore

        cd fichier_input/case
        ../../build/costo


