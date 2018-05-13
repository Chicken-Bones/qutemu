# QUT Atrial Emulation
Research code for exploring the effects of fibrosis on atrial electrophysiology using [CHASTE](http://www.cs.ox.ac.uk/chaste/)

### Current Results 
[![Fast Pacing Re-entry](https://i.imgur.com/oltTiSp.jpg)](https://youtu.be/NyiyFBW7jeI?t=1m10s "Fast Pacing Re-entry")
![Fibrotic pattern](https://i.imgur.com/hU4CiFO.png)
![Electrical wave overlaid with fibrosis](https://i.imgur.com/QzXCRhi.png)
### Architecture
The software is split into 3 directories, based on language.
* [chaste](./chaste/README.md) contains the C++ user project and setup instructions for running simulations
* [MATLAB](./MATALB/README.md) contains a perlin noise implementation with functions for using noise to generate fibrotic patterns in atrial tissue
* pyscripts contains scripts and macros to aid in visualisation using [ParaView](https://www.paraview.org/)

*Unfortunately the heart model used cannot currently be provided due to IP reasons. Please contact the repository owner with a request if you want to extend this research.
