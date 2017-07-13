# IF-Microscopy-Assistant
## For the analysis and modelling of Immunofluorescence (IF) cell microscopy
This project aims to construct a cloud-based service to help with the analysis and morphological modeling in Immunofluorescence cell microscopy.
<br><br>
## Code Files
<br>
-> create_masks.mat (MATLAB R2017a) <br>
Inputs: a .tiff file containing a set of rgb images in which [red channel = nucleus] and [green channel = cytoplasm. <br>
Output: cell array of rgb images or video frames, binary mask of single cells, binary mask of cell clusters and binary mask of nucleus. <br><br>
-> create_dataSet.mat  (MATLAB R2017a) <br>
Inputs: cell array of rgb images or video frames and binary mask of single cells. <br>
Outputs: --- work in progress ---
