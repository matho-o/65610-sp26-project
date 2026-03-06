# MIT 6.5610 Applied Cryptography - Spring 2026 Project
Team members: Adhitya Mangudy Venkata Ganesh, Nicola Lawford, Ishan Satish Pednekar, Li Xuan Tan

**Problem:** Performing singular value decomposition on a matrix encrypted using fully homomorphic encryption.

# Setup
<TODO: figure out OS quirks?>

This script builds and installs OpenFHE and OpenFHE-Python for FHERMA challenges.
```sh
git clone --recursive https://github.com/Fherma-challenges/toolchain
cd toolchain
sh prepare_env_debian.sh
```

<TODO: figure out if pip install works?>

# Solution validation
- Install Docker Desktop 4.37.1 or later.
- Run ```sh docker pull yashalabinc/fherma-validator```.
- `cd` into `fherma-svd` directory, which contains a `app` directory where the solution is located, and a `tests` directory with `.json` files representing test cases.
- Run ```sh docker run -ti -v /:/fherma yashalabinc/fherma-validator --project-folder=/app --testcase=/tests/test_case.json```.
- A `result.json` file will be generated in the project folder.

# Resources

- `openfhe-python` implementation template from [OpenFHE](https://github.com/fairmath/fherma-challenges/tree/main/templates/openfhe-python)
- [FHERMA Challenge - Singular Value Decomposition](https://fherma.io/challenges/68fb3d896f81f4f6f684aac2/overview)
- [Docker validator](https://hub.docker.com/r/yashalabinc/fherma-validator) provided by FHERMA
- [OpenFHE-Python library](https://github.com/openfheorg/openfhe-python)
- [FHERMA challenge toolchain](https://github.com/Fherma-challenges/toolchain)
- [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/index.html)
