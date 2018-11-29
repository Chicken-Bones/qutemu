
# Chaste Project
### Installation
1. Install Chaste 3.4, ideally using https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage
    * Chaste can be installed seamlessly in windows 10 using WSL and the Ubuntu app
    * **Replace `git clone` with the develop branch of fork, https://github.com/Chicken-Bones/Chaste**
    ```
    git clone -b develop https://github.com/Chicken-Bones/Chaste.git Chaste
    ```
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
| --- | --- | --- | --- |
| `-meshfile` | `<path>` | `!!required!!` | path to atrial mesh (wthout the .node extension)
| `-outdir` | `<dir>` | `ChasteResults` | sets output directory to `testoutput/<dir>`
| `-loaddir` | `<dir>` |  | Simulation will be resumed* from a state in `testoutput/<dir>`. `-meshfile` will be ignored
| `-savedir` | `<dir>` |  | Simulation will be saved in `testoutput/<dir>`
| `-nodes` | `<nodelist>`<br>`<nodefile>` || Restrict output nodes (by number in .node file). A comma separated list of nodes to output or a file where each entry is a single line containing a node number. |
| `-vtk` ||| Enable vtk output |
| `-duration` | `<length>` | `5` | length of simlation (ms) |
| `-interval` | `<period>` | `5` | Data output/logging interval for results.h5 and results.[p]vtk (ms) |
| `-odet` | `<step>` | `0.02` | maximum ODE integration step (ms) |
| `-pdet` | `<step>` | `<odet>` | maximum PDE integration step (ms) |
| `-cell` | `maleckar`<br>`courtemanche_sr`<br>`courtemanche_caf` | `courtemanche_sr` | cell model to use |
| `-sinus` | `<timelist>`<br>`<timefile>` || A comma separated list or newline separated file containing the stimulus times. Specifying this option will ignore `-psinus` and `-nsinus`. `-dsinus` can be used to add a constant to time values in this option. |
| `-dsinus` | `<delay>` | `0` | Delay before the first sinoatrial node trigger (ms) |
| `-psinus` | `<period>` | `500` | Period of sinoatrial trigger (ms) |
| `-nsinus` | `<num>` | `4` | Number of sinoatrial triggers  |
| `-extra` | `<timelist>`<br>`<timefile>` || `-sinus` for ectopic stimulus. Times are relative to the end of the sinus stimulus.  |
| `-dextra` | `<delay>` | `300` | Delay between last sinoatrial trigger and first ectopic trigger (ms) |
| `-pextra` | `<period>` | `150` | Period of ectopic trigger (ms) |
| `-nextra` | `<num>` | `6` | Number of ectopic triggers  |
| `-base_cond` | `<num>` | `1.75` | Base conductivity value (in Chaste's units). Not used if anatomical locations specificed in the .ele file |
| `-ar` | `<num>` | `9.21` | Ratio of conductivity values between longitudinal and transverse (default value here and for -base_cond corresponds to Chaste's traditional (1.75, 0.19, 0.19) conductivity) |
| `-condmod` | `<file>` || A file with one line per element containing conductivity multipliers
| `-svi` ||| Enables state-variable interpolation https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/StateVariableInterpolation
| `-activation` | `<threshold>` | `-40` | Activation threshold used for generating snapshots (mV). |

\* Duration is measured from the start of the whole simulation, not from the end of the loaded simulation, so a longer value must be provided for continuation
