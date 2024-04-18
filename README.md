# ![AutoPas](https://raw.githubusercontent.com/AutoPas/AutoPas/master/docs/graphics/AutoPasLogo_Large.svg "Title")

AutoPas is a node-level auto-tuned particle simulation library developed
in the context of the [**TaLPas**](http://www.talpas.de) project.
[![CI Status](https://github.com/AutoPas/AutoPas/actions/workflows/TestSuites.yaml/badge.svg)](https://github.com/AutoPas/AutoPas/actions/workflows/TestSuites.yaml)

## Documentation
The documentation can be found at our website:
 <https://autopas.github.io/doxygen_documentation/git-master/>

Alternatively you can build the documentation on your own:
* requirements: [Doxygen](http://www.doxygen.nl/)
* `make doc_doxygen`

## Examples
As AutoPas is only a library, it is not able to run simulations by itself.
We have, however, included a few example proxy applications in the **examples** directory.
The examples include:
* [md-flexible](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible): Molecular dynamics simulations with single centered Lennard-Jones particles.
* [Smoothed particle hydrodynamics simulations](https://github.com/AutoPas/AutoPas/blob/master/examples/sph)

## Using AutoPas
Please look at our [user documentation pages](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc).

## Developing AutoPas
Please look at our [contribution guidelines](https://github.com/AutoPas/AutoPas/blob/master/.github/CONTRIBUTING.md).

## Acknowledgements
This work was financially supported by:
* the Federal Ministry of Education and Research, Germany, project “Task-based load balancing and auto-tuning in particle simulations” (TaLPas) 8 , grant numbers 01IH16008A and 01IH16008B.
* TODO: 3xa
* TODO: MaST
* TODO: WindHPC

## Papers to cite
* F. A. Gratl, S. Seckler, H.-J. Bungartz and P. Neumann: [N Ways to Simulate Short-Range Particle Systems: Automated Algorithm Selection with the Node-Level Library AutoPas](https://www.sciencedirect.com/science/article/abs/pii/S001046552100374X), In Computer Physics Communications, Volume 273, 2022. ([BibTeX](https://mediatum.ub.tum.de/export/1638766/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1638766))
* F. A. Gratl, S. Seckler, N. Tchipev, H.-J. Bungartz and P. Neumann: [AutoPas: Auto-Tuning for Particle Simulations](https://ieeexplore.ieee.org/document/8778280), In 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), Rio de Janeiro, May 2019. ([BibTeX](https://mediatum.ub.tum.de/export/1535848/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1535848))
* S. Seckler, F. Gratl, M. Heinen, J. Vrabec, H.-J. Bungartz, P. Neumann: [AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning](https://www.sciencedirect.com/science/article/abs/pii/S1877750320305901), In Journal of Computational Science, Volume 50, 2021. ([BibTeX](https://mediatum.ub.tum.de/export/1595680/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1595680))
