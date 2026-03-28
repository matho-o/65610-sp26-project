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

If you are using Windows, run on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) and use the following before running `prepare_env_debian.sh` to ensure line endings are correct:

```sh
sed -i 's/\r$//' prepare_env_debian.sh
```

<TODO: figure out if pip install works?>

# Solution validation
- Install Docker Desktop 4.37.1 or later. If you are using windows, ensure [WSL backend is configured](https://docs.docker.com/desktop/features/wsl/#prerequisites)
- Run ```docker pull yashalabinc/fherma-validator```.
- From the root of this repo, run ```sudo docker run -ti -v $(pwd):/fherma yashalabinc/fherma-validator --project-folder=/fherma/openfhe-python --testcase=/fherma/fherma-svd/tests/test_case.json```.
- A `result.json` file will be generated in the project folder.
- **Troubleshooting**: The `fherma-svd` directory should contian a `tests` directory with `.json` files representing test cases. The `openfhe-python` directory should contain your python solution. See [Docker validator](https://hub.docker.com/r/yashalabinc/fherma-validator) provided by FHERMA.
- **Note**: We may need to edit the project structure to match the challenge format.

# Resources

- `openfhe-python` implementation template from [OpenFHE](https://github.com/fairmath/fherma-challenges/tree/main/templates/openfhe-python)
- [FHERMA Challenge - Singular Value Decomposition](https://fherma.io/challenges/68fb3d896f81f4f6f684aac2/overview)
- [Docker validator](https://hub.docker.com/r/yashalabinc/fherma-validator) provided by FHERMA
- [OpenFHE-Python library](https://github.com/openfheorg/openfhe-python)
- [FHERMA challenge toolchain](https://github.com/Fherma-challenges/toolchain)
- [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/index.html)
