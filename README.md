<p align="center">
  <a href="https://imr.org">
    <img src="docs/imr.jfif" alt="IMR Banner" width="500"/>
  </a>
</p>

<p align="center">
  <a href="https://doi.org/10.1016/j.jmps.2017.12.006" target="_blank">
    <img src="https://zenodo.org/badge/doi/10.1016/j.jmps.2017.12.006.svg" />
  </a>
  <a href="https://www.gnu.org/licenses/gpl-3.0.html">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" />
  </a>
</p>

# Welcome to IMRv2

IMRv2 is an advanced computational tool for **Inertial Microcavitation Rheometry (IMR)**, enabling the characterization of soft materials under high strain-rate conditions. IMR correlates the evolution of bubble pressure and stress fields in a material with kinematic observations obtained from high-speed videography.

IMRv2 originated from the **Computational Flow Group at Brown University** and is now actively developed by researchers at multiple institutions:
- [Christian Franck](https://directory.engr.wisc.edu/me/Faculty/Franck_Christian/) (University of Wisconsin-Madison)
- [David Henann](https://vivo.brown.edu/display/dhenann) (Brown University)
- [Eric Johnsen](https://me.engin.umich.edu/people/faculty/eric-johnsen/) (University of Michigan)
- [Jon Estrada](https://me.engin.umich.edu/people/faculty/jon-estrada/) (University of Michigan)
- [Mauro Rodriguez](https://vivo.brown.edu/display/mrodri97) (Brown University) - **Lead Developer**
- [Spencer Bryngelson](https://comp-physics.group/) (Georgia Tech)
- [Jin Yang](https://sites.utexas.edu/yang) (University of Texas at Austin)

For questions, contact [Mauro Rodriguez](mailto:mrodri97@brown.edu) or request to join the IMR Slack workspace.

## Features

IMRv2 offers several key capabilities for modeling and analyzing inertial microcavitation:

- **Bubble Dynamics Simulation**: Implements both finite difference and spectral methods for discretizing partial differential equations, as detailed in:
  - Estrada et al., "High Strain-rate Soft Material Characterization via Inertial Cavitation," [Journal of the Mechanics and Physics of Solids, 2017](https://doi.org/10.1016/j.jmps.2017.12.006).
  - Warnez and Johnsen, "Numerical modeling of bubble dynamics in viscoelastic media with relaxation," [Physics of Fluids, 2015](https://doi.org/10.1063/1.4928860).
- **High Strain-rate Characterization**: Analyzes high-speed video data to extract material properties under dynamic loading.
- **MATLAB Integration**: Provides scripts and functions for simulation control, data analysis, and visualization.
- **Experimental Validation Support**: Facilitates comparison between simulation results and experimental observations.

## Getting Started

### Prerequisites

IMRv2 is implemented in **MATLAB** and requires the following dependencies:

- MATLAB (version R2021a or newer recommended)
- Optimization Toolbox (optional but recommended for parameter fitting)
- Image Processing Toolbox (for analyzing high-speed video data)

### Installation

Clone the repository and navigate to the IMRv2 directory:

```bash
git clone https://github.com/InertialMicrocavitationRheometry/IMRv2.git
cd IMRv2
```

Add IMRv2 to your MATLAB path:

```matlab
addpath(genpath('IMRv2'))
savepath
```

### Running an Example

To run a basic simulation of bubble dynamics, execute the following command in MATLAB:

```matlab
run_example.m
```

For additional examples, see the `examples/` directory.

## Documentation

Comprehensive documentation, including a user guide, can be found in the `docs/` folder and on the IMRv2 website.

## Citation

If you use IMRv2 in your research, please cite:

```bibtex
@article{Estrada_2017,
  title   = {High Strain-rate Soft Material Characterization via Inertial Cavitation},
  author  = {J. B. Estrada and C. Barajas and D. L. Henann and E. Johnsen and C. Franck},
  journal = {Journal of the Mechanics and Physics of Solids},
  year    = {2017},
  volume  = {112},
  pages   = {291–317},
  doi     = {10.1016/j.jmps.2017.12.006}
}
```

For additional references, see the **IMR Bibliography** section above.

## License

IMRv2 is released under the [GNU General Public License v3.0](LICENSE).

## Acknowledgments

IMRv2 development has been supported by funding from the National Science Foundation (NSF), the Department of Defense (DOD), and other research institutions.

## Contact

For issues, bug reports, or feature requests, please open an issue on GitHub or contact the maintainers.

