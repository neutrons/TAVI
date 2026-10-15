# TAVI, Triple-Axis data VIsualization toolkit
TAVI can be installed locally as a standard conda package and has been deployed on ORNL's analysis cluster.

<img width="2562" height="1193" alt="image" src="https://github.com/user-attachments/assets/143b0ca6-9cb9-45dc-8112-9b4d6e542720" />

## Installation

Create and activate a virtual environment with [Pixi](https://pixi.sh/).
Prerequisites: Pixi installation e.g. for Linux:

```bash
curl -fsSL https://pixi.sh/install.sh | sh

```

Download the repository. Setup/Update the environment

```bash
pixi install
```

Enter the environment

```bash
pixi shell

```

The Tavi environment is activated and the application is ready to use.

*Alternatively, stable versions of tavi are provided as conda package: [Tavi Package Installation Instructions](https://anaconda.org/neutrons/tavi)

Start the application

```bash
tavi
```

Documentation [tavi.readthedocs.io](https://tavi.readthedocs.io/)

[![CI](https://github.com/neutrons/TAVI/actions/workflows/test_and_deploy.yml/badge.svg?branch=next)](https://github.com/neutrons/TAVI/actions/workflows/test_and_deploy.yml)
[![codecov](https://codecov.io/gh/neutrons/TAVI/graph/badge.svg?token=AYB1X932FV)](https://codecov.io/gh/neutrons/TAVI)
