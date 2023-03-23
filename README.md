# kMC-DD
kinetic Monte Carlo / Dislocation Dynamics code for edge dislocation motion with super-jogs

The project consists of three parts: 

## Table of Contents

- [Install](#install)
- [Usage](#usage)
- [Examples](#example)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)


## Install

This project is developed in [MATLAB](https://www.mathworks.com/products/matlab.html). 


## Usage

### input parameters
[input_HEA.m](input_HEA.m) is the file where simulation parameters and material parameters are specified, start with this file to make sure the set up is correct.

### kMC/DD simulation
[dd3d_hea.m](dd3d_hea.m) is the main file, it is overall a kMC algorithm, and between each kMC events, DD simulations serve to fill the time interval. As the simulation runs, a movie of dislocation line motion will be created in real time. Also, the strain-vs-simulation time will be written into txt files for later analysis.



## Example
https://user-images.githubusercontent.com/56003395/227333746-c0249328-7a25-4444-899a-64329008e719.mp4


## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by research group at UCLA, Johns Hopkins University, Hongkong City University and Pennsylvania University.


## License

[MIT](LICENSE) © Sicong He
