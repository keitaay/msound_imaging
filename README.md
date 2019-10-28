# Imaging with mSOUND
This repository contains an implementation of the nonlinear ultrasound simulator mSOUND for imaging applications.

mSOUND is a transient-mixed- and frequency-specific-mixed-domain (i.e. nonlinear) acoustic wave propagation simulator. This open-source, MATLAB-specific toolset was originally developed by Jiangian Gu and Yun Jing. More information, including preprint copies of the corresponding journal publication  can be found at: https://github.com/m-SOUND/mSOUND

The tutorials and examples packaged with mSOUND were primarily intended for beamforming and applied imaging scenarios, including photoacoustics and shear shock wave modeling. This toolset extends the capabilities of mSOUND by defining and imposing additional, medical imaging-related parameters into mSOUND variables and data formats. This addition may allow this toolbox to act as a nonlinear imaging simulator with capabilities similar to more popular software such as Field-II or k-wave.

## General Structure
This repository is structured as follows:

`mSOUND` - not included; should be downloaded either in the same folder level as others in this repository, or elsewhere.
`shortcut` - shortcut functions to quickly define parameters relevant for both mSOUND and ultrasonic imaging.
`visualization` - shortcut functions to generate simulation-relevant images, useful for debugging

`linear_image.m` - runs a simulation of a plane wave transmit on a designated field of scatterers

## Shortcut
`msound_medium` - using a `set_grid` custom object, this function creates a `medium` structure that incorporates information about the material properties being simulated in mSOUND.

`msound_xdcr` - using the `medium` structure among others, this function creates a structure `xdcr` that contains parameters about the transducer plane in a single, compact object.

`msound_excite` - using the `medium` structure among others, this function creates a structure `exci` that contains matrices depicting pressure waves transmitted from a transducer plane.

`msound_beamform` - using the outputs from mSOUND-native function `Forward2D`, this function applies dynamic-receive beamforming on the collected, channel-based pressure waveforms.

## Visualization
`testMedium.m` - using the `medium` structure required by mSOUND, generate 2-D images that display heterogeneous features that will be simulated.
