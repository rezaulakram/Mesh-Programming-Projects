2.1 Installation and getting started:
Followed the steps mentioned in Exercise2.pdf. Downloaded, Configured and build the project. Loaded the model with the exe.

2.2 Subsampling
Subsampled the points uniformly. Iterate over the src points and calculated the distance from selected points. If a point is not near(below threshold) to any of the selected points, the point is included to the selected list.

2.3 Bad Pairs Rejection
Calculated the distance between src and destination and rejected the pairs whose distance is more the 3 times the median length.
Calculated the cosine angle between normals and rejected the pairs which have angle more than 60 degree.

2.4 Point-2-Point optimization
Deduced the value of A and b matrix and filled the proper values for Cholesky decomposition.
DIFFICULTY FACED : in deriving the values for matrix A as it was not given.

2.5 Point-2-Surface optimization
Filled the values in matrix A as suggested in the slides.

Checked the program by running multiple times with different models.


SUBMITTED BY :
REZAUL AKRAM BARBHUIYA
USC ID: 4811-2972-49

