# Polar_Dispr_index
Matlab code that calculates  PI (polarization index) and DI (dispersion index)

as described in: 
http://www.cell.com/action/showMethods?pii=S2211-1247%2811%2900019-2


PI or polarization index is calculated as:
sqrt( centroid_of_molecule - centroid_of_cytoplasm )^2 / Radius_of_gyration

Rg or Radius of gyration is caclculated by root-mean-square distance of all
pixels of the cell from the centroids (centroid as axis) to normalize cell shape

DI or dispersion index is calculated as:
variance_of_molecule / second_moment_of_cytoplasm
(i.e., variance of molecule / variance of cytoplasm)
