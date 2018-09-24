# HHQTemplate #
This is a repository with the relevant code for creating Glauber initial conditions, hydrodynamically evolving those conditions forward with SONIC, and simulating heavy quark diffusion through the hydrodynamic medium. The code is structured in such a way that creating the initial conditions and simulating diffusion is performed with ROOT, whereas the hydro evolution is in pure C++. To this end, additional routines are required to "translate" ROOT histogram files for temperature, fluid velocity, etc. information to a .dat format used and stored by hydro. Each routine is described below in the order it should be run in; condor_run.csh is the primary job driver and can illustrate this process as well.

I'm using this version of the sPHENIX software:
source /opt/sphenix/core/bin/sphenix_setup.csh -n pro.538

## TGlauberMC ##
This directory contains runglauber_v3.1.C, which a modified version of the code provided at https://tglaubermc.hepforge.org/. The main difference is the addition of a function "runAndSmearNucleons" which generates a TNTuple with event properties and a 2D histogram with the wounded nucleon density, smeared out by a Gaussian. The TNTuple can be used for getting parameters like Psi2 or Epsilon2, whereas the TH2's are used by hydrodynamics as the initial energy density. When generating initial conditions, in addition to specifying the number of events and the collision system, the nucleon-nucleon cross section and nucleon smearing parameter must be set. If a minimum or maximum impact parameter is also desired, those can be specified in the last 2 arguments.
After saving the Glauber initial conditions (e.g., to "outFile_pPb.root"), Hist2Txt.C must be run as a ROOT macro. This is a simple script that saves each histogram to its own ROOT file in the appropriate initedFiles directory, as well as an identical copy in .dat format. There are a few parameters at the top of the script to edit, but otherwise this *should* be pretty straightforward.

TODO: You need to specify your grid parameters in Glauber before running. This can be done in the code itself. In runAndSmearNucleons, the number of bins and maximum x, y values need to be specified. (Note - xmax = 5 means your grid will go from x=-5 to x=5, for a total width in x of 10. The units are fm.) The energy scaling parameter will be appropriately changed when you modify these parameters. SONIC takes in energy dimensionlessly, i.e. it takes the energy density times (bin width)^4. It is important to calculate the lattice spacing again in you SONIC input - see setup_events.sh below.

## setup_events.sh ##
setup_events.sh will setup $2-$1 events in the directory it is executed in. For event$i, it will create the directory event$i, in which all of the processing will occur. The main script and code will be copied over and should be modified to match the params.txt file, which should be modified in turn to match the initial conditions. This includes changing the system "SPC", number of bins, but also the grid spacing/grid size, nucleon-nucleon cross-section. The smoothing option is also important for larger systems with lumpy initial conditions, e.g. He3+Au or He4+Au. The script requires 5 arguments total however; $3 tells the script if SONIC is to be run, $4 is if Txt2Hist will be run, and $5 is if diffusion will be run.

TODO: Before setting up your events, this needs to happen:<br>
        TGlauberMC: this should already have been run, but if you haven't generated initial conditions and converted to .dat, this is a good time to do so<br>
        Edit setup_events.sh: user-defined directories should be changed (see TODO items in file!)<br>
      After setting up your events, these files need to be appropriately modified:<br>
        condor_run.csh<br>
        params.txt: values for SIGMANN, NUMT, AT, SNAPUPDATE should be modified. ***This is probably the most important file to modify correctly!*** <br>
          SIGMANN = nucleon nucleon cross-section, I'm not sure this really is important but it can't hurt to make sure it is correct.<br>
          NUMT = just the number of spatial bins in x,y; e.g. 200.<br>
          AT = this is the spatial lattice spacing, in units of GeV^-1. Following the notation from TGlauberMC, AT should be given by AT = 5.0677 * 2*xmax / bins.
          SNAPUPDATE = this tells SONIC how many times to save the current temperature, velocity fields, etc. If you only want dNch_dy for example, set this number high to save disk space (e.g. 20000), but if you do care about what the temperature looks like at each timestep, you can set it lower to, e.g., 100.
          SMOOTHING = whether SONIC should do smoothing. Generally this should be 1 for Glauber "lumpy" IC's. I've found this isn't needed for 8TeV p+Pb, but it was critically important for 200GeV He3+Au and larger systems.

        Txt2Hist.C: values for N
        diffusion.C: values for NGridX, NGridY, timestep_length_fm should be modified (see TODO)

## condor_run.csh ##
This is the primary job driver script. It executes each module in the correct order, organizing input files appropriately and copying over needed file. For instance, if a custom file (such as an initial condition) is provided, a line in condor_run.csh can be added to copy over the corresponding .dat file. (Alternatively, of course, you can modify the line below to copy over the .dat file when instantiating the event.) Other files, like diffusion.C and Txt2Hist.C, are copied at run time from the parent directory so that necessary changes can be made once and propagate to all events.
This script takes 6 arguments; the first is simply the run directory, which must be in the same directory that condor_run.csh is in. $2 is the event number, so if the folder is event510, then $2 should be 510. $3 is the job number, so if you are running events 40 to 80, then event40 would have $3 = 0 (this is only for plotting purposes in Txt2Hist, so if you don't care about generating a .gif animation, this isn't important.) Finally, $4, $5, and $6 tell the script whether to run SONIC, Txt2Hist, and diffusion, respectively, and can take the values 0 (don't run) or 1 (do run).
One last comment: there is a parameter at the top of the script called MOVE_TO_SCRATCHDISK which tells the code whether to move the event to a scratchdisk on the node the job is being executed on. This can help reduce stress on the sPHENIX disk and is recommended for bulk jobs.

## init, generate, vh2, b3d, and b3d-analyze ##
These are the main sonic routines. b3d requires access to hdf5, which is installed on my own home directory ~jouellette/ and should be publically accessible. These are setup to output in the subdirectory event$i/logdir. vh2 is the most intensive since it runs the main hydro evolution. It can be monitored in the log file by controlling the UPDATE parameter in params.txt.

## Txt2Hist.C ##
This file is meant to be run after the entire sonic routine has been completed. It takes the output from vh2 and generates a root file from the FOdata, Tcontour, and (optionally) EDprofile with temperature and fluid velocity distributions at each timestep which is intended for use in the diffusion routine. It can also save .gif images of the fluid at each timestep which can be later combined into animations.

## diffusion.C ##
The main diffusion routine is meant to be run after Txt2Hist.C. As input, it requires several different pieces of information. Primarily, it requires the hydro output which comes from Txt2Hist.C. It also needs a pT spectrum for heavy quarks (this can be used defined or generated in Pythia - heavy_quark_pt.root is the default pT spectrum and is provided here, however it only contains a ccbar pT spectrum but ideally a bbbar will be added soon). Lastly, it requires a wounded nucleon distribution (the initial energy density should be fine as a root file) as well as the name of the distribtion as a TH2D.

## submit_jobs.job ##
This is the condor submission script. If you want to change how some run is processing, you can edit the arguments passed along in here.


## Additional important files ##

heavy_quark_pt.root contains the heavy quark initial pT spectra. This was generated by forcing ccbar pair creation in Pythia, and plotting their pT values.
quarkdist.root is the file copied to each run directory when setting up the events. This contains the initial energy density as a TH2 object, which is interpreted as a heavy quark spatial density function, and is used in diffusion to simulate heavy quark pair production locations.
params.txt is the primary parameters file used by sonic. There are several values in this file that need to be modified for sonic to run correctly, as listed above. The most important values are the number of bins and the bin width, but controlling the SNAPUPDATE parameter also helps keep the disk size in check.
log is the directory containing all of the output from condor_run.csh. Note that the output from the SONIC run is put in logdir within the event itself since it is more verbose.
