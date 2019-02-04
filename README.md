# SNDrift readme


Yorck Ramachers (Warwick)
Last updated February 4, 2019

The SNDrift code platform attempts to simulate electron drift in the
SuperNEMO tracker gas, defined as Helium (95%), Ethanol (4%) and
Argon (1%) gas mixture. A simplified tracker geometry is available
as a ROOT file using ROOT geometry objects. That geometry maps 
exactly to the geometry used in the finite element code COMSOL
which calculates the fine-grained electrostatic fields in that
tracker geometry. These fields are also provided in form of a 
ROOT file, to be read by the code. These electric fields
can be scaled to realistic bias values on anode wires with
a command line option.
Transport follows the same algorithm used in Garfield++ and other
microscopic charge transport codes, see below. Cross sections 
have been obtained from MagBoltz and are available for reading
in form of a ROOT file.

## Description

To build it, do

``` console
$ ls
CMakeLists.txt  data  examples  FindROOT.cmake  include  src  testing  utils

$ mkdir build
$ cd build
$ cmake ..
...
$ make
...
... If you are developing the module, you can test it by calling
$ make test
...
... or obtain more detail on the tests and launch in the build directory
$ ctest -V
```

The build will create the `transportlib.so` shared library and (currently)
the scan.exe executable in the build directory (from which you can run
the executable, no problem).

## Utilities and Data

scan.exe: comand line options are:

``` console
$ ./build/scan.exe --help
collection scan command line option(s) help
	 -x , --xstart <x-coordinate start [cm]>
	 -y , --ystart <y-coordinate start [cm]>
	 -b , --bias <Anode bias in Volt>
	 -s , --seed <random number seed offset>
$
```

The data directory contains all the required input data in ROOT file
format. 

``` console
$ ls data/
sntracker_driftField.root  trackergasCS.root  trackergeom.gdml
```

Scripts to make (and view) the ROOT geometry and convert the COMSOL 
field output to a ROOT file are in the utils/ folder. The original 
COMSOL file is too large to ship with SNDrift which somewhat restricts
experimentation with this code.

The geometry chosen for SNDrift should reflect all important cases 
for drifting charge in the tracker. Three complete tracker rows of 
9 cells each are implemented following the coordinates from the 
Monte-Carlo geometry. Field-wire coordinates and anode wires as well
as dimensions as programmed there in geometry version 4 have been 
reproduced. The distances to the calorimeter wall and the source 
foil are also following the MC geometry. The distances to the 
X-Wall lack 4 mm each side with this SNDrift distance being smaller
than the MC. That was a mistake but should not cause critical issues
when starting charges in that (almost) field-free space.

This simplified 2D geometry should allow to start ionization charges at 
any point in the geometry, i.e. around the perimeter of the tracker
and at any point in between. Considering only one out of 111 cell 
rows between the perimeter rows should not matter due to symmetry. Also
the de-facto removal of the z-ccordinate should not matter given the 
symmetry. However, the code is inherently working in 3D, i.e. the ROOT
geometry is in fact 2 cm thick in z such that diffusion in z is 
permitted to roughly the tracker z-resolution. Still, starting points
are all at z=0 by default.

The script to run MagBoltz is also included but clearly requires a 
Garfield++ installation. It is adapted from the file Examples/GasFile/ 
example shipped with Garfield++. If copied in that directory and 
compiled with the existing makefile in that directory, one would 
reproduce the text file containing all the requested microscopic cross 
sections. Converting those cross sections to a suitable ROOT file 
can then be achieved with the script cs2root.py in the utils folder.

## Implementation

The 3D charge transport uses the null collision method from Skullerud
like most transport codes. Disadvantages considered so far 
when using tabulated cross sections are the lack of calculated 
angular distributions for elastic scattering (assumed isotropic hence)
and ignoring excitation and Penning transfer cross sections so far.
These could be obtained from Magboltz but are left out for simplicity
in the current version. Only elastic scattering and ionization is 
considered here. [WIP: more here if needed].

Elastic scattering should completely dominate scattering and hence 
contribute most to drift speed simulations. The ionization was included
in order to see whether Geiger avalanches could be obtained, allowing
for maybe better determination of drift speeds. Neglecting 
additional scatter processes and better angular scattering behaviour
risks systematic deviations from absolute drift speeds. Relative 
drift speeds as a function of charge origin might still be fine
regardless. The latter is the main target for this simulation code.

Transporting charge avalanches can quickly get out of hands 
computationally hence this code uses multi-threading with a 
hard-coded restriction to four threads maximum. That can be 
changed in the code if faster execution is required. Likely, this is 
not needed though.

The scan.exe application code is in the examples/ directory and represents 
a typical example of using the transport library. Other applications can be 
considered and likely will be created later on. Output to disk would 
be mandatory at some point and creating a set of starter charges 
initially could also lead to more efficient computations. One might wish 
to run many starter charges from roughly the same location to maybe obtain 
firmer results on drift times and repeat to cover several locations of
interest. The combined set of results would then result in a tracker 
drift model to be compared to previous work on the drift model. Special
cases could then be considered separately like exotic starting points
clustered around the circumference of a field-wire. How much slower, if 
at all, arrive charges from outside a cell, near the calo walls or from 
between foil and first tracker cells?

