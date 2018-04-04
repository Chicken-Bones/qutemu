# Chaste Project
### Installation
1. Install Chaste 3.4, ideally using https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage
    * Chaste can be installed seamlessly in windows 10 using WSL and the Ubuntu app
2. Ensure chaste compiles using cmake, become accustomed to the configure/make cycle and run some of the builtin cardiac tests
    * You can run `make help` from the build directory to see a list of targets
    * Chaste will take several hours to build all the tests, so if you want a quick setup and verify, just run  
```
make -j4 TestMonodomain3dRabbitHeartTutorial
ctest -R TestMonodomain3dRabbitHearTutorial
```
3. symlink this folder to /chaste-src/projects/qutemu
    * If you're using WSL, make a directory link to have it show on both OS
```
mklink /D <chaste-src-dir>\projects\qutemu <repo-root>\chaste\qutemu
```

### Building
1. Run cmake to reconfigure, whenever files are added/removed
2. Run make, with the desired target
```
cd <chaste-build-dir>
cmake <chaste-src-dir>
make -j4 project_qutemu
```
The binaries can be found in `<build-dir>/projects/qutemu`

### Running
Simply invoke `AtrialFibrosis` from the command line with at least a `-meshfile` switch
Alternatively, run it via mpi with `mpirun AtrialFibrosis ...`

If you want to copy the output to another directory, you can use `ldd AtrialFibrosis` to get a list of dependencies and copy all the chaste libraries and `libchaste_project_qutemu.so`


### Command Line Arguments
| Switch | Params | Default | Description |
| - | - | - | - |
| `-meshfile` | `<path>` | `!!required!!` | path to atrial mesh (wthout the .node extension)
| `-outdir` | `<dir>` | `ChasteResults` | sets output directory to `testoutput/<dir>`
| `-loaddir` | `<dir>` |  | Simulation will be resumed* from a state in `testoutput/<dir>`. `-meshfile` will be ignored
| `-savedir` | `<dir>` |  | Simulation will be saved in `testoutput/<dir>`
| `-novis` ||| Suppress vtk output
| `-duration` | `<length>` | `5` | length of simlation (ms) |
| `-interval` | `<period>` | `5` | Data output/logging interval for results.h5 and results.[p]vtk (ms) |
| `-odet` | `<step>` | `0.02` | maximum ODE integration step (ms) |
| `-pdet` | `<step>` | `<odet>` | maximum PDE integration step (ms) |
| `-cell` | `maleckar`<br>`courtemanche_sr`<br>`courtemanche_caf` | `maleckar` | cell model to use |
| `-dsinus` | `<delay>` | `0` | Delay before the first sinoatrial node trigger (ms) |
| `-psinus` | `<period>` | `500` | Period of sinoatrial trigger (ms) |
| `-nsinus` | `<num>` | `4` | Number of sinoatrial triggers  |
| `-dextra` | `<delay>` | `300` | Delay between last sinoatrial trigger and first ectopic trigger (ms) |
| `-pextra` | `<period>` | `150` | Period of ectopic trigger (ms) |
| `-nextra` | `<num>` | `6` | Number of ectopic triggers  |
| `-nextra` | `<num>` | `6` | Number of ectopic triggers  |
| `-nodes` | `<nodelist>`<br>`<nodefile>` || Restrict output nodes (by number in .node file). A comma separated list of nodes to output or a file where each entry is a single line containing a node number. 
| `-fibrosis` | `<file>` || A file with one line per element, indicating how fibrotic the tissue is. Conductivity is scaled by `1 - fibrosis`

\* Duration is measured from the start of the whole simulation, not from the end of the loaded simulation, so a longer value must be provided for continuation
