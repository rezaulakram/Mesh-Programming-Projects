Surface Representation and Geometric Modelling
Exercise 3 - Implicit Surface Reconstruction

1. Gone through the readme and class notes for the theory. Followed the related papers.

2.1 Hoppe’s distance from tangent planes :
 a) Iterate over all the points and find the closest point.
 b) Assuming the equation of the plane is given by ax + by +cz + d = 0, we found the value of d as :
		d = -(ax + by + cz) = - (n dot p)
		where n = [a b c]'
 c) Using the value of d, we calculated the distance of the point from the plane as :
		(ax0 + by0 + cz0 + d) / sqrt(a*a + b*b + c*c)
 d) Returned the distance which is used for surface reconstruction using marching cubes.
 
2.2 Signed distance function via triharmonic RBFs :
 a) Calculated the epsilon as 1% of bounding box diagonal
 b) Stored the 2n centers
 c) Now we created a gmmMatrix A of size 2n * 2n and stored the kernel value between each pair of points and centers.
 d) Solve the equation using ImplicitRBF::solve_linear_system function
 
 Difficulty faced in finding the value of epsilon. Followed the paper - "Surface Reconstruction from Scattered Point via RBF Interpolation on GPU" which says 1 percent of bounding box diagonal.
 
 
2.3  Implement compactly supported Wendland basis function:
 a) Implemented the wendlend kernel using equation given in class.
 
 Difficulty faced : Evaluating the value of sigma. Putting small values of sigma were giving erroneous results. Checked the following papers :
 "Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree" by Holger Wendland
 "An Introduction to Scattered Data Approximation" by Greg Coombe
 "Surface Reconstruction from Scattered Point via RBF Interpolation on GPU"

 CONFUSION : The first paper doesn't use any sigma. Hence I used sigma = 1. The second paper uses sigma value as support radius. But is silent about evaluating it. Tried with the sigma calculation with third paper but it doesnot give proper result. It looks like it is evaluating the sigma for gaussian kernel rather than Wendland kernel.

SUBMITTED BY :
REZAUL AKRAM BARBHUIYA
USC ID: 4811-2972-49