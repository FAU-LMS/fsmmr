# Frequency-Selective Mesh-to-Mesh Resampling for Color Upsampling of Point Clouds
This repository provides Matlab-code for color point cloud upsampling using Frequency-Selective Mesh-to-Mesh Resampling (FSMMR). The program is free of charge for personal and scientific use. Any commercial use is prohibited. If you use our software, don't forget to cite our work [Heimann2021].

## Usage

Usage is very simple. Just open Matlab and execute:

	main.m


At the moment, one exemplary point cloud from the 3D Color Mesh Dataset [3DColMesh] is incorporated as a minimial working example. The original data is stored in the 'data' folder. During computation, a 'result' folder is generated where the final point clouds are stored. Be aware, that the final stored point cloud shows only the part of the point cloud that did not have color information assigned in the beginning of the processing. Thus, it is only a subset of the original point cloud. 


## Requirements

Written and tested for MATLAB R2017b.

## Literature
[Heimann2021] V. Heimann, A. Spruck, Kaup A., "Frequency-Selective Mesh-to-Mesh Resampling for Color Upsampling of Point Clouds," Proceedings IEEE 23rd International Workshop on Multimedia Signal Processing Tampere,2021  
DOI: https://doi.org/10.1109/MMSP53017.2021.9733445  
[3DColMesh] Anass Nouri, Christophe Charrier, Olivier Lézoray - Greyc 3D Colored Mesh Database, Technical report, 2017. https://hal.archives-ouvertes.fr/hal-01441721
