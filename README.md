# IF-Microscopy-Assistant
## For the analysis and modelling of Immunofluorescence (IF) cell microscopy
This project aims to construct a cloud-based service to help with the analysis and morphological modeling in Immunofluorescence cell microscopy.
<br><br>
## Code Files
<br>
-> create_masks.m (MATLAB R2017a) <br>
Inputs: a .tiff file containing a set of rgb images in which [red channel = nucleus] and [green channel = cytoplasm. <br>
Output: cell array of rgb images or video frames (uint8), binary mask of single cells, binary mask of cell clusters and binary mask of nucleus. <br><br>
-> create_sequences.m  (MATLAB R2017a) <br>
Inputs: cell array of rgb images or video frames (uint8) and binary mask of single cells. <br>
Outputs: cell array [s x f] containg the Pixel Idx List of cell sequences. s = no. unique sequences, f = no. video frames. 
<br><br><br><br>
**Reference**
<br>
This code is part of the project: <br>
'Tracking of temporally occluded or overlapping structures in live cell microscopy' by Miguel Esteras. <br>
Department of Computer Science at City, University of London
