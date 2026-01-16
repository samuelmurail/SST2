[![Documentation Status](https://readthedocs.org/projects/sst2/badge/?version=latest)](https://sst2.readthedocs.io/en/latest/?badge=latest)
[![DOI:10.1021/acs.jctc.5c00950](http://img.shields.io/badge/DOI-10.1021/acs.jctc.5c00950-B31B1B.svg)](https://doi.org/10.1021/acs.jctc.5c00950)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13772542.svg)](https://doi.org/10.5281/zenodo.13772542)

<img src="https://raw.githubusercontent.com/samuelmurail/SST2/master/docs/source/logo.jpeg" alt="AF2 Analysis Logo" width="400" style="display: block; margin: auto;"/>

# Simulated Solute Tempering 2 (SST2)

## Description

This repository contains the source code for the Simulated Solute Tempering 2 (SST2) algorithm, as described in our paper ["Simulated Solute Tempering 2: An Efficient and Practical Approach to Protein Conformational Sampling and Binding Events."](https://doi.org/10.1021/acs.jctc.5c00950) SST2 is a novel enhanced sampling method for molecular dynamics (MD) simulations that combines the strengths of Simulated Tempering (ST) and Replica Exchange with Solute Tempering 2 (REST2).

* Source code repository on [gihub](https://github.com/samuelmurail/SST2)
* Documentation on readthedocs [![Documentation Status](https://readthedocs.org/projects/sst2/badge/?version=latest)](https://sst2.readthedocs.io/en/latest/?badge=latest)
* Article published in *J. Chem. Theory Comput.* [![DOI:10.1021/acs.jctc.5c00950](http://img.shields.io/badge/DOI-10.1021/acs.jctc.5c00950-B31B1B.svg)](https://doi.org/10.1021/acs.jctc.5c00950)
* Trajectories on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13772542.svg)](https://doi.org/10.5281/zenodo.13772542)

## Algorithm Overview

SST2 aims to overcome limitations in conformational space exploration often encountered in traditional MD simulations, especially for systems with high energy barriers. SST2 builds on previous methods like REST2 and ST:

- **REST2:** Like REST2, SST2 scales solute-solute and solute-solvent interactions to enhance sampling.
- **ST:**  SST2 adopts the concept of single simulation traveling across different temperatures (or scaling factors) from ST, but applies it specifically to solute tempering as in REST2.

These features enable SST2 to effectively sample conformational space and provide valuable insights into the thermodynamics and kinetics of biomolecular systems.


## library Main features:

- ST simulation using Park and Pande weights calculation, and Phong et al. on the flight weights calculation.
- REST2 potential energy implementation
- **Efficient SST2 Sampling:** By selectively scaling solute interactions and exchanging replicas, SST2 significantly enhances sampling efficiency compared to traditional MD simulations.
- Binary scripts to run ST and SST2 simulations starting from an amino acid sequence or a `pdb` file.
- **Flexibility:** The algorithm allows users to adjust various parameters, including the number of rungs, the range of scaling factors, and the exchange frequency, to optimize performance for different systems.
- **Open Source:** The code is made available as open source, facilitating further development, adaptation, and application by the research community.

> [!WARNING]  
> **Charmm36 derived forcefields are not yet supported**: As charmm36 forcefield use specific dihedral terms like CMAP correction, charmm forcefield is not yet implemented in SST2, this feature will be added in future releases. If you need this feature, implementation should be straightforward, please contact the author.

## Implementation

This implementation of SST2 utilizes the [OpenMM](https://openmm.org/) molecular dynamics library. The code is written in Python and utilizes OpenMM's custom forces and integrators to achieve the desired functionality.

## Installation

Please refer to the [installation guide](https://sst2.readthedocs.io/en/latest/installation.html) in the documentation for detailed instructions on how to install SST2 and its dependencies.

To install the last SST2 version, you can use pip:

```bash
pip install git+https://github.com/samuelmurail/SST2.git
```

## Contributing

`SST2` is an open-source project and contributions are welcome. If
you find a bug or have a feature request, please open an issue on the GitHub
repository at [https://github.com/samuelmurail/SST2](https://github.com/samuelmurail/SST2). If you would like
to contribute code, please fork the repository and submit a pull request.

## Author

* [Samuel Murail](https://samuelmurail.github.io/PersonalPage/), Associate Professor - [Université Paris Cité](https://u-paris.fr), [BFA](https://bfa.u-pariscite.fr/).

See also the list of [contributors](https://github.com/samuelmurail/SST2/contributors) who participated in this project.

## License

This project is licensed under the GNU General Public License v2.0 - see the ``LICENSE`` file for details.

## Acknowledgments

The authors would like to express their gratitude to Fabio
Sterpone and Samuela Pasquali for their valuable contributions
to the discussion surrounding the theoretical aspects and
implementation of SST2. 
We'd like to thanks Benjamin Ye ([mdcraft](https://mdcraft.readthedocs.io)) which mdcraft's topology module was used in this project.
Furthermore, we extend our appreciation to [Peter Eastman](https://github.com/peastman) and, more broadly, to the
entire [OpenMM](https://openmm.org/) development team, for their efforts in making OpenMM and the
simulated tempering implementation accessible as an open-source resource.

## Citation

If you use this code in your research, please cite our paper:

Stratmann D, Moroy G, Tuffery P and Murail S. Simulated Solute Tempering 2: An Efficient and Practical Approach to Protein Conformational Sampling and Binding Events. [J. Chem. Theory Comput. (2025). 21(21), 10705–10718.](https://doi.org/10.1021/acs.jctc.5c00950)

```
@article{stratmann_simulated_2025,
        title = {Simulated Solute Tempering 2: An Efficient and Practical Approach to Protein Conformational Sampling and Binding Events},
        volume = {21},
        issn = {1549-9618},
        url = {https://doi.org/10.1021/acs.jctc.5c00950},
        doi = {10.1021/acs.jctc.5c00950},
        shorttitle = {Simulated Solute Tempering 2},
        pages = {10705--10718},
        number = {21},
        journaltitle = {J. Chem. Theory Comput.},
        author = {Stratmann, Dirk and Moroy, Gautier and Tuffery, Pierre and Murail, Samuel},
        urldate = {2025-11-12},
        date = {2025-11-11},
        note = {Publisher: American Chemical Society},
}
```