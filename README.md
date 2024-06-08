# CoopeRIS

CoopeRIS is a Plexe/Veins framework enabling the simulation of reconfigurable intelligent surfaces (RIS) for connected vehicles.

> [!NOTE]
> The software is currently in its infancy. It has been recently released to the public but proper documentation and example are still being developed.

## Building

CoopeRIS depends on the Plexe software ecosystem, so SUMO, OMNeT++, Veins, and of course Plexe.
With respect to SUMO, please consider installing version 1.18.0, whereas for OMNeT++ use 6.0.1.

### GNU scientific library
CoopeRIS requires you need to install the GNU scientific library.
On a Ubuntu/Debian system, please install it with
```bash
sudo apt install libgsl-dev
```
while on macOS you can either use MacPorts
```bash
sudo port install gsl
```
or homebrew
```bash
brew install gsl
```

### Veins
CooperRIS works on a customized version of Veins.
Clone it using [this link](https://github.com/michele-segata/veins/tree/cooperis) and be sure to choose the `cooperis` branch:
```bash
git clone https://github.com/michele-segata/veins/tree/cooperis
cd veins
git checkout cooperis
```
Then simply compile Veins as usual:
```bash
./configure
make
```

### Plexe
CooperRIS example is included within Plexe (starting from version 3.1.3).
Clone it using [this link](https://github.com/michele-segata/plexe) and be sure to choose the 3.1.3 version:
```bash
git clone https://github.com/michele-segata/plexe
cd plexe
git checkout -b plexe-3.1.3-work plexe-3.1.3
```
Then simply compile Plexe as usual:
```bash
./configure
make
```

### CoopeRIS (Multithread support, default)

CoopeRIS can be built with multithread or GPU support to accelerate the
computation of the RIS. By default, CoopeRIS is built with multithread support.
First, clone the repository:

```bash
git clone https://github.com/michele-segata/cooperis
cd cooperis
```
To build it, please configure it indicating the path to the GSL include and lib folders.
You can do so in the following way:
```bash
./configure --with-gsl-include=/opt/local/include --with-gsl-lib=/opt/local/lib
```

Please make sure to change the GSL paths to match your owns.

Finally, simply type

```bash
make
```

to build CoopeRIS.

> [!NOTE]
> You can specify
> the number of compute threads to use for the computation of the RIS gain by
> using the `*.**.nicRis.phyRis.maxWorkerThreads` param in the `omnetpp.ini` file
> of your simulation. If left unspecified, the number of threads will be
> set to the number of available cores on your machine.

> [!WARNING]
> If your machine has a very large number of cores, you might want to limit the
> number of threads to the number of real cores, excluding hyperthreading ones to
> achieve best performance.

### CoopeRIS (GPU support)

CoopeRIS can be built with GPU support to accelerate the computation of the RIS.
Currently, CoopeRIS supports Cuda and OpenCL frameworks.

#### Cuda support

To enable GPU acceleration with Cuda, you need to have the Cuda toolkit installed on your
machine. Please refer to the
[official NVIDIA website](https://developer.nvidia.com/cuda-downloads) to
download the toolkit. Once installed, you can build CoopeRIS with Cuda support
by specifying the `--with-cuda` and the path to the Cuda include and lib
folders, as follows:

```bash
./configure --with-gsl-include=/opt/local/include --with-gsl-lib=/opt/local/lib --with-cuda --with-cuda-include=/opt/local/include --with-cuda-lib=/opt/local/lib
```

Please make sure to change the GSL and Cuda paths to match your owns.

Finally, simply type.

```bash
make
```

to build CoopeRIS with Cuda acceleration.

> [!NOTE]
> If you have multiple Cuda-enabled devices on your machine, you can specify the
> number of the Cuda device to use by using the `*.**.nicRis.phyRis.cudaDeviceId`
> param in the `omnetpp.ini` file of your simulation. If left unspecified, device
> 0 will be used. You can discover the available devices on your machine by using
> the `nvidia-smi` utility command.

#### OpenCL support

To enable GPU acceleration with OpenCL, you need to have the OpenCL framework
installed on your machine. Please refer to the
[official Khronos website](https://www.khronos.org/opencl/) to download the
framework. Once installed, you can build CoopeRIS with OpenCL support by
specifying the `--with-opencl` and the path to the OpenCL include and lib
folders, as follows:

```bash
./configure --with-gsl-include=/opt/local/include --with-gsl-lib=/opt/local/lib --with-opencl --with-opencl-include=/opt/local/include --with-opencl-lib=/opt/local/lib
```

Please make sure to change the GSL and OpenCL paths to match your owns.

Finally, simply type.

```bash
make
```

to build CoopeRIS with OpenCL acceleration.

> [!NOTE]
> If you have multiple OpenCL-enabled devices or platforms on your machine, you
> can specify the number of the OpenCL device and platform by using the
> `*.**.nicRis.phyRis.openclDeviceId` and `*.**.nicRis.phyRis.openclPlatformId`
> params in the `omnetpp.ini` file of your simulation. If left unspecified, device
> 0 and platform 0 will be used. You can discover the
> available devices and platforms on your machine by using the `clinfo` utility command.

> [!WARNING]
> CoopeRIS includes a large set of unit tests to ensure that all mathematical computations are properly working.
> Such unit tests are located in the `subprojects/cooperis_catch` folder.
> We will soon add instructions on how to run them.

### Plexe example

The `plexe_cooperis` subproject example includes a simple intersection scenario where a car is static on a road, while a second one is travelling being tracked by the RIS.
To run the project simply build the subproject:
```bash
cd plexe/subprojects/plexe_cooperis
source setenv
./configure
make
```
Then run the example with:
```bash
cd examples/plexe_cooperis
plexe_cooperis_run -u Cmdenv -c TrackingTIntersection -r 0
```

### Matlab

> [!WARNING]
> This section is still incomplete, we are working on the full release.

### Scientific documentation

> [!WARNING]
> We are in the process of publishing an accepted paper on CoopeRIS, which will explain to user the model in detail.
> We will update the repository as soon as the article will be available online.
