Surface Representation and Geometric Modelling
Exercise 4 - Surface Smoothing

3.1 Uniform Laplace Curvature and Smoothing
a) Implemented calc_uniform_mean_curvature() 

Calculated uniform Laplace operator over all the vertices by iterating over the vertices and then 1-ring neighbourhood for each of the vertices. Calculated the centroid for the neighbouring vertices which is used to calculated Laplacian.
Stored the mean curvature in vunicurvature_ vertex property.

b) Implemented uniform_smooth()

Used the same procedure as in (a) to caluculate the Laplacian and then moved each of the vertex.
Updated the normals.

3.2 Triangle Shapes

Implemented calc_triangle_quality()
Calculated the ratio of circumradius r and smallest side to calculate the triangle quality and stored it in tshape_ property.

3.3 Laplace-Beltrami curvature and smoothing
Study calc_weights to understand the data structures.

a)Implemented calc_mean_curvature()
Iterated over the 1-ring neighbourhood using half edges and calculated the Laplace Betrami operator on each of the vertices.
Stored the vcurvature_ property.

b)Implemented SmoothingViewer::smooth()
Calculated Laplace Beltrami operator in the same way as (a) but used sum of weights for normalization.

Difficulty faced in understanding the half edge iteration.

3.4 Gaussian curvature

Iterated over each of the vertices and calculated the angle around each of them. Approximated Gaussian curvature and stored it in vgausscurvature_ property.


SUBMITTED BY :
REZAUL AKRAM BARBHUIYA
USC ID: 4811-2972-49