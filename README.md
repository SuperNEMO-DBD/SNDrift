# SNDrift readme


Yorck Ramachers (Warwick)

Last updated February 18, 2019

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
have been obtained from Magboltz and are available for reading
in form of a ROOT file.

## Description

Depends on: ROOT, tested with version 6.14.xx

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
two executables in the build directory (from which you can run
the executables, no problem).

## Utilities and Data

scan.exe: comand line options are:

``` console
$ ./build/scan.exe --help
collection scan command line option(s) help
	 -x , --xstart <x-coordinate start [cm]>
	 -y , --ystart <y-coordinate start [cm]>
	 -b , --bias <Anode bias in Volt>
	 -p , --pressure <tracker gas pressure [mbar]>
	 -s , --seed <random number seed offset>
	 -o , --outputFile <FULL PATH ROOT FILENAME>
$
```

mcdrift.exe: command line options are:

``` console
$ ./build/mcdrift.exe --help
Monte-Carlo scan command line option(s) help
	 -x , --xstart <x-coordinate start [cm]>
	 -y , --ystart <y-coordinate start [cm]>
	 -c , --ncharges <number of starter charges at x,y>
	 -b , --bias <Anode bias in Volt>
	 -s , --seed <random number seed offset>
	 -n , --nsim <number of Monte Carlo simulations>
	 -p , --pressure <tracker gas pressure [mbar]>
	 -d , --dataDir <FULL PATH Directory to data file>
	 -o , --outputFile <FULL PATH ROOT FILENAME>
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

The script to run Magboltz is also included but clearly requires a 
Garfield++ installation. It is adapted from the file Examples/GasFile/ 
example shipped with Garfield++. If copied in that directory and 
compiled with the existing makefile in that directory, one would 
reproduce the text file containing all the requested microscopic cross 
sections. Converting those cross sections to a suitable ROOT file 
can then be achieved with the script cs2root.py in the utils folder.

## Implementation

The 3D charge transport uses the null collision method from Skullerud
like most transport codes. Disadvantages considered so far 
when using tabulated cross sections are that excitation and Penning 
transfer cross sections are ignored so far.
These could be obtained from Magboltz but are left out for simplicity
in the current version. Only elastic scattering and ionization is 
considered here. [WIP: more here if needed].

Elastic scattering should completely dominate scattering and hence 
contribute most to drift speed simulations. The ionization was included
in order to see whether Geiger avalanches could be obtained, allowing
for maybe better determination of drift speeds. Better angular 
scattering behaviour for Helium has been achieved by implementing
the formulae from Phys. of Plasmas, 19 (2012) 093511. Scattering 
from rare gas components, Ethanol and Argon, remains isotropic.

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

Allowing for multiple charges at the start of a Monte-Carlo run improves 
efficiency for cluster jobs. The option '-c' for mcdrift.exe allows for 
this feature. The requested number of charges is randomly placed in a cube 
of 2 micro metre side length (+-1 mum) around the starter position, i.e. 
small enough not to impact drift times for those N charges. The N longest 
drift times are then stored in the results. This way, one can get a statistic 
of N charges started times a simulation repetition of nsim times the number 
of jobs on a cluster to achieve a better statistics. Tests appear to 
confirm that this mode of running saves time compared to single 
starter charges.

## Simulation geometry

As described above, three complete tracker cell units of 9 cells each 
are in the geometry with the correct distances around the perimeter to 
the main calorimeter and the source foil. Just the distance to the X-wall 
is 4 mm short on either side. In coordinates, the COMSOL geometry dictates 
what the ROOT geometry has to follow. The project student who made the 
COMSOL model had for some reason a preference for the lower right corner 
in the Cartesian plane...

Therefore, the absolute locations of boundaries and tracker wires sit at 
positive x-coordinate values but negative y-coordinate values. The top left 
corner of the COMSOL volume sits at the absolute origin (0,0). The COMSOL 
chamber extends in y-direction to 43.4 cm. The centres of the nearest 
field-wires to the main calorimeter wall sit at y = -7 mm. The centres of 
the field-wires nearest to the source foil sit at y = -40.3 cm, after which 
there is a 2 mm gap and a 29 mm gap to the source foil, like in the Falaise 
MC geometry version 4.x. The 9 cells of a tracker cell unit 
extend in -y direction. 

In x-direction, one X-wall is at x = 0, the nearest field-wires sit at x = 14 mm. 
A good guide-point for orientation can also be the most top-left anode wire 
at x = 3.6 cm and y = -2.9 cm. Staying close to it is the default location for 
the executable example. The three layers of 9-cell units extend in x-direction.

Two image files are available in the utils folder, made in COMSOL, and show the 
geometry and electrostatic field structure in more detail than can be achieved
with the descriptions above.

Ref: H. Schindler, Garfield++ User manual (2019), http://garfieldpp.web.cern.ch/garfieldpp/documentation/

S. F. Biagi, Magboltz, http://magboltz.web.cern.ch/magboltz

