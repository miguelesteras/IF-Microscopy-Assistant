# IF-Microscopy-Assistant
For the analysis and modelling of Immunofluorescence (IF) cell microscopy.<br>
This project aims to construct a cloud-based service to help with the analysis and morphological modeling in Immunofluorescence cell microscopy.
<br><br>
## MATLAB (R2017a) code
<br>
-> create_masks.m <br>
Converts a .tiff file containing a set of rgb images in which [red channel = nucleus] and [green channel = cytoplasm, into a cell array of rgb images or video frames (uint8), a binary mask of single cells, a binary mask of cell clusters and a binary mask of nucleus. <br><br>

-> create_sequences.m <br>
Returns cell array [s x f] containg the Pixel Idx List of each individual cell from the single cell mask. Each row contains a unique sequence of frames of a single cell. No.columns = No.video frames.<br><br>

-> create_InfoArrays.m <br>
Returns for each cell in cell sequence array; the center of mass, max radius (max rho in polar coord.), angle to rotate max rho pixel north, cotour polar coordinates and cell mask polar coordinates.<br><br>

-> transform_Boundary.m <br>
Transforms cell mask coordinates into a small set of [x,y] cell contour coordenates of arbitrary size, using the matlab function _boundary_.<br><br>

-> transform_DynamicImage.m
Transforms cell mask coordenates into a grayscale image containing the contour of the cell from all input frames. The intensity of the pixel contour is homogenious for one frame an increases in time from frame to frame. <br><br>

-> transform_fourierDescriptors.m
Transform cell mask coordinates into a feature vector containing fourier descriptors (based on code by Tobias Pohlen). <br><br>

-> transform_GFDforKNN.m <br>
Transform cell contour into Generic Fourier Descriptors (GFD) (based on code by Frederik Kratzert). The GFD form a dictionary of sequences to which the query sequence is matched up using _KNNsearch_.<br><br>

-> transform_RhoDescriptors.m <br>
Transforms cell contour coordenates into a vector containing the max rho coordinate for the set pixel at specific evenly spaced angles.<br><br>

<br><br><br><br>
## Reference
<br>
This code is part of the project: <br>
'Tracking of temporally occluded or overlapping structures in live cell microscopy' by Miguel Esteras. <br>
Department of Computer Science at City, University of London
