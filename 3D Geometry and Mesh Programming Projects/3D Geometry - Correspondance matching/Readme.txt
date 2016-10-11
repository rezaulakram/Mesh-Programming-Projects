Digital Geometry Processing
SEMESTER PROJECT
IMPLEMENTATION OF CORRESPONDANCE PART IN PAPER: 

DEFORMATION TRANSFER FOR TRIANGLE MESHES

BY :

ROBERT W. SUMNER 

and

JOVAN POPOVIC´

Computer Science and Artificial Intelligence Laboratory
Massachusetts Institute of Technology

(Note : The project uses cholmod and umfpack libraries to solve sparse matrics equations. THey are not included in project folder)

Command : Deformation.exe horse_ref.obj camel_ref.obj horse_camel.cons

Deformation transfer applies the deformation exhibited by a source triangle mesh onto a different target triangle mesh. In this paper the author proposes a novel technique to deform target triangle mesh from the deformation in source triangle mesh.
The user builds a correspondence map between the triangles of the source and those of the target by specifying a small set of vertex markers.
The technique described in this paper creates many to many correspondence map between triangles in source mesh and target mesh, and then transfers the deformation from one mesh to other.
In this project I have implemented the many to many correspondence calculation part.

The project involves following parts :

Computing deformation gradient for each triangle using the equation (4) in the paper.
Q = V'V^-1

We then minimize the following terms :

Smoothness term :

Es(v1 ... vn) = Ti - Tj
for all j triangles adjacent to i.

Deformation identity :

Ei(v1 ... vn) = Ti - Tj
This term prevent the smoothness term to create drastic change in the whole mesh

Closest Point Term: 
We calculate the closest point in the target mesh for each vertex in the source mesh and we tend to move the source vertices towards them
Ec(v1 ... vn, c1 ... ci) = Vi - Ci

We calculate a weighted equation for each of the terms and add them. We then minimize the whole equation :

wsEs + wiEi + WcEc

In the first iteration Wc = 0 and Ws = 1 and Wi = 0.01. From the next iteration we increase Wc from 0.5 to 1.0. This is different from the paper as I am finding good resultin this range.

We then calculate the threshold for source and destination meshes, take the maximum one and add all the traingles which matches following criteria :

If the centroid of the source and destination triangle is less then threshold and angle between the face normal is less than 90 degree I add it to the correspondance list.

For each triangle of the deformed source we compute the closest compatible triangle in the correspondance list. We repeat the same procedure for the target mesh and find a compatible triangle in source mesh.

This ensures that all the corresponding pairs are in our list.


SUBMITTED BY :
REZAUL AKRAM BARBHUIYA
USC ID: 4811-2972-49