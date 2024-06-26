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

### CoopeRIS
To build CoopeRIS, clone first the repository:
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
