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

# ORCD Engaging Cluster
- Get an account following the [startup guide](https://orcd-docs.mit.edu/getting-started/) and SSH into a Centos 7 node.
- `salloc -p mit_normal -c 12 --mem=54G` to request a machine with 12 CPU/54GB RAM just like the FHERMA configuration, on the `mit_normal` ORCD cluster.
- Once logged in, run the following:
    - `module purge`
    - `module load gcc/12.2.0`
    - `module load cmake/3.27.9`
- Install OpenFHE by compiling from source as per official instructions. For CMake, use the command
```
CC=$(which gcc) CXX=$(which g++) cmake \
    -DCMAKE_CXX_FLAGS="-static-libstdc++ -static-libgcc" \
    -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
    -DCMAKE_SHARED_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
    -DCMAKE_INSTALL_PREFIX=~/openfhe-install ..
    ..
```
where the last line points to a userspace install location (since we don't have sudo access on Engaging). Run an example or two to verify successful installation.
- Default ORCD modules don't have recent-enough Python versions. Load ORCD's conda distribution with `module load miniforge/25.11.0-0` and verify that `python --version` returns 3.10+.
- Create a virtual environment, activate it and `pip install "pybind11[global]"`.
- Install OpenFHE-Python by compiling from source as per official instructions. For CMake, use 
```
CC=$(which gcc) CXX=$(which g++) cmake \
   -DCMAKE_CXX_FLAGS="-static-libstdc++ -static-libgcc" \
   -DCMAKE_SHARED_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
   -DCMAKE_PREFIX_PATH=~/openfhe-install ..
```
- Add OpenFHE-Python install to PYTHONPATH: `export PYTHONPATH=/home/<  your login here  >/openfhe-python/build:$PYTHONPATH`. Python examples should run fine now.
- To use the validator, use the cluster's [Singularity/Apptainer](https://orcd-docs.mit.edu/software/apptainer/) install (compatible with Docker images). Run `singularity pull fherma-validator.sif docker://yashalabinc/fherma-validator` to pull the Dockerfile and save it locally, then `singularity run --bind $(pwd):/fherma ./fherma-validator.sif --project-folder=/fherma/openfhe-python --testcase=/fherma/fherma-svd/tests/test_case.json` to run the validation script.


# Resources

- `openfhe-python` implementation template from [OpenFHE](https://github.com/fairmath/fherma-challenges/tree/main/templates/openfhe-python)
- [FHERMA Challenge - Singular Value Decomposition](https://fherma.io/challenges/68fb3d896f81f4f6f684aac2/overview)
- [Docker validator](https://hub.docker.com/r/yashalabinc/fherma-validator) provided by FHERMA
- [OpenFHE-Python library](https://github.com/openfheorg/openfhe-python)
- [FHERMA challenge toolchain](https://github.com/Fherma-challenges/toolchain)
- [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/index.html)
