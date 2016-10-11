#include "DeformationViewer.hh"
#include "deformationTransfer.hh"
#include "correspondance.hh"
#include <unordered_map>


cholmod_common global;
Matrix9x4* triangleMatrices;
Matrix9x1* columnVectors;
//int* freeVertIndices;
int phase2mode = NORMALMODE;
int state = 0;

void _cholmod_error_handler (int status, const char *file, int line, const char *message)
{
	printf("Error occurred in :\t File: %s\n\t Line : %d\n\t Status : %d\n\t Message :%s\n\n", file, line, status, message);
}

void Initialize()
{
	cholmod_start(&global);
	columnVectors = NULL;
	triangleMatrices = NULL;
	srcWindow->inverseSurfaceMap = NULL;
	deformedSrcWindow->inverseSurfaceMap = NULL;
	targetWindow->inverseSurfaceMap = NULL;
	global.error_handler = _cholmod_error_handler;
}
void Deinitialize()
{
	cholmod_finish(&global);
}



/* Get the coord index of v1 v2 v3*/
long DeformationViewer::GetVertexCoordIndex(unsigned int vertexId, unsigned int coord)
{
	if(isConstraintVertex(mesh_.vertex_handle(vertexId)))
		return -1;
	return (freeVertIndices[vertexId] * 3 + coord);
}

/* Get the coord index of v4*/
long DeformationViewer::GetNormalCoordIndex(unsigned int triangleId, unsigned int coord)
{
	unsigned int nfree = mesh_.n_vertices() - constraintList.size();

	return (((nfree + triangleId) * 3) + coord);

}



void DeformationViewer::GetTriangleVerticesIndex(int faceId, int vertexIndexTable[][4])
{
	Mesh::FHandle fhandle = mesh_.face_handle(faceId);
	Mesh::FVIter fvIter;
	Mesh::VHandle vertices[3];
	int k = 0;

	for(fvIter = mesh_.fv_begin(fhandle); fvIter; ++fvIter, k++)
	{
		vertices[k] = fvIter.handle();
	}

	for(int i =0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			vertexIndexTable[i][j] = GetVertexCoordIndex(vertices[j].idx(), i);
		}
		vertexIndexTable[i][3] = GetNormalCoordIndex(fhandle.idx(), i);
	}
}

int Cholmod_SparseMatrix_Entry(cholmod_triplet *T, int i, int j, double x)
{
	if (T->nnz >= T->nzmax && 
        !cholmod_reallocate_triplet(2*(T->nzmax), T, &global))
    {
        return 0;  /* reallocation failed */
    }

    /* append subscript */
    if (T->itype == CHOLMOD_LONG)
    {
        ((UF_long*)T->i)[T->nnz] = i;
        ((UF_long*)T->j)[T->nnz] = j;
    }
    else
    {
        ((int*)T->i)[T->nnz] = i;
        ((int*)T->j)[T->nnz] = j;
    }
    
    ((double *)(T->x))[T->nnz] = x;


#define __DT_MAX(a,b) (((a)>(b))?(a):(b))
    T->nrow = __DT_MAX(T->nrow, (size_t)(i+1));
    T->ncol = __DT_MAX(T->ncol, (size_t)(j+1));
#undef __DT_MAX

    T->nnz += 1;

    return 1;
}

int DeformationViewer::AddClosestPointTerm(cholmod_triplet* M, cholmod_dense* C,double weight, int row)
{
	int tgt_vertex;
    int i_x, i_y, i_z;
	Mesh::VertexIter vIter, v_end = mesh_.vertices_end();

	for(vIter = mesh_.vertices_begin(); vIter != v_end; ++vIter)
    {
        if (!isConstraintVertex(vIter.handle()))
        {
			Mesh::Point p = mesh_.point(vIter);
            i_x = GetVertexCoordIndex(vIter.handle().idx(), 0);
            i_y = GetVertexCoordIndex(vIter.handle().idx(), 1);
            i_z = GetVertexCoordIndex(vIter.handle().idx(), 2);

			float location[3];
			tgt_vertex = mesh_.property(ClosestTargetVertexIdx, vIter.handle());
			targetWindow->getVertexLocation(tgt_vertex, location);

            Cholmod_SparseMatrix_Entry(M, row, i_x, weight);
			((double *)(C->x))[row] = (weight * location[0]);
            row++;

            Cholmod_SparseMatrix_Entry(M, row, i_y, weight);
			((double *)(C->x))[row] = (weight * location[1]);
            row++;

            Cholmod_SparseMatrix_Entry(M, row, i_z, weight);
			((double *)(C->x))[row] = (weight * location[2]);	
            row++;
        }
    }
	return row;
}


void AddTermToMatrix(cholmod_triplet* M, cholmod_dense* C, int triangleMatrixId, double weight, int row, int vertexIndexTable[][4])
{
	int MatrixRow = 0;
	Matrix9x4 *m = &triangleMatrices[triangleMatrixId];
	Matrix9x1 *c = &columnVectors[triangleMatrixId];

	for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++, row++, MatrixRow++)
        {
            for (int col = 0; col < 4; col++)
            {
                int i_var = vertexIndexTable[i][col];

                /* append this term if it is a free vertex */
                if (i_var != -1) {
                    Cholmod_SparseMatrix_Entry(M, row, i_var, weight * ((*m)[MatrixRow][col]));
                }
            }
            /* element of right hand side vector */
			double value = ((double *)(C->x))[row];
			((double *)(C->x))[row] = value + (weight * ((*c)[MatrixRow]));
        }
    }

}

void AddColToDenseMatrixDT(cholmod_dense *C, int row, double c[][3])
{
	int MatrixRow =0;
	double weight = 1.0;
	for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++, row++, MatrixRow++)
        {
			double value = ((double *)(C->x))[row];
			((double *)(C->x))[row] = value + (weight * (c[i][j]));
		}
	}
}

void AddTermToMatrixDT(cholmod_triplet* M, int triangleMatrixId, int row, int vertexIndexTable[][4])
{
	int MatrixRow = 0;
	double weight =1.0;
	Matrix9x4 *m = &triangleMatrices[triangleMatrixId];
	Matrix9x1 *c = &columnVectors[triangleMatrixId];

	for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++, row++, MatrixRow++)
        {
            for (int col = 0; col < 4; col++)
            {
                int i_var = vertexIndexTable[i][col];

                /* append this term if it is a free vertex */
                if (i_var != -1) {
                    Cholmod_SparseMatrix_Entry(M, row, i_var, weight * ((*m)[MatrixRow][col]));
                }
            }
        }
    }

}

void AddTermToMatrix(cholmod_triplet* M, cholmod_dense* C, Matrix9x4 m, Matrix9x1 c ,double weight, int row, int vertexIndexTable[][4])
{
	int MatrixRow = 0;
	
	for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++, row++, MatrixRow++)
        {
            for (int col = 0; col < 4; col++)
            {
                int i_var = vertexIndexTable[i][col];

                /* append this term if it is a free vertex */
                if (i_var != -1) {
                    Cholmod_SparseMatrix_Entry(M, row, i_var, (weight * m[MatrixRow][col]));
                }
            }
			
            /* element of right hand side vector */
			double value = ((double *)(C->x))[row];
			((double *)(C->x))[row] = value + (weight * c[MatrixRow]);
        }
    }

}


int AddIdentityInTriangleMatrix(cholmod_triplet* M, cholmod_dense* C, int triangleMatrixId, double weight, int row)
{
	int triangleVertIdx[3][4];
	Matrix9x1 *c = &columnVectors[triangleMatrixId];
	Matrix9x1 CIdentity = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	for(int i = 0; i < 9; i++)
		CIdentity[i] += (*c)[i];

	srcWindow->GetTriangleVerticesIndex(triangleMatrixId, triangleVertIdx);
	Matrix9x4 *m = &triangleMatrices[triangleMatrixId];
	AddTermToMatrix(M, C,(*m), CIdentity , weight, row, triangleVertIdx);

	//global.print = 5;
	//cholmod_print_triplet(M, "MATRIX", &global);
	//cholmod_print_dense(C, "MATRIX", &global);

	return (row + 9);
}


int AddSmoothnessInTriangleCouple(cholmod_triplet* M, cholmod_dense* C, int triangleMatrixId, int AdjtriangleMatrixId, double weight, int row)
{
	int triangleVertIdx[3][4];
	int AdjTriangleVertIdx[3][4];

	srcWindow->GetTriangleVerticesIndex(triangleMatrixId, triangleVertIdx);
	srcWindow->GetTriangleVerticesIndex(AdjtriangleMatrixId, AdjTriangleVertIdx);
	AddTermToMatrix(M, C, triangleMatrixId, weight, row, triangleVertIdx);
	AddTermToMatrix(M, C, AdjtriangleMatrixId, -weight, row, AdjTriangleVertIdx);
	global.print = 5;
	//cholmod_print_triplet(M, "MATRIX", &global);
	return (row + 9);
}

int DeformationViewer::AddSmoothnessCoefficients(cholmod_triplet* M, cholmod_dense* C,double weight, int row)
{
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	Mesh::FFIter ffIter;


	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{

		for(ffIter = mesh_.ff_begin(fIter.handle());ffIter; ++ffIter)
		{
			row = AddSmoothnessInTriangleCouple(M, C, fIter.handle().idx(), ffIter.handle().idx(), weight, row);
		}
	}
	return row;
}

int DeformationViewer::AddIdentityCoefficients(cholmod_triplet* M, cholmod_dense* C,double weight, int row)
{
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{
		row = AddIdentityInTriangleMatrix(M, C, fIter.handle().idx(), weight, row);
	}
	return row;
}

void LoadCorrespondanceList(char* filename, std::vector<CorrespondanceEntry>& corres_list)
{
	FILE* file = fopen(filename,"r");
	int ncorres = 0;
	fscanf(file,"%ld", &ncorres);
	
	for(int i = 0; i < ncorres; i++)
	{
		CorrespondanceEntry corres_entry;
		int src_id, tgt_id, distance;
		fscanf(file, "%ld, %ld, %lf", &corres_entry.src_triangle_id, &corres_entry.target_triangle_id, &corres_entry.distance);
		corres_list.push_back(corres_entry);
	}
	fclose(file);
}

void SaveCorrespondanceList(char* filename, std::vector<CorrespondanceEntry> corres_list)
{
	FILE* file = fopen(filename,"w");
	long length = corres_list.size();
	fprintf(file,"%ld\n",length);
	for(int i = 0; i < length; i++)
	{
		fprintf(file,"%ld, %ld, %lf\n",corres_list[i].src_triangle_id, corres_list[i].target_triangle_id,corres_list[i].distance);
	}
	fclose(file);
}

void SolveLinearEquations(cholmod_sparse *A, cholmod_dense *b, cholmod_dense *x)
{
	double *null = (double *) NULL ;
	int i ;
	void *Symbolic, *Numeric ;
	(void) umfpack_di_symbolic(A->nrow, A->ncol, (int*)A->p, (int *)A->i, (double*)A->x, &Symbolic, null, null) ;
	(void) umfpack_di_numeric ((int*)A->p, (int *)A->i, (double*)A->x, Symbolic, &Numeric, null, null) ;
	umfpack_di_free_symbolic (&Symbolic) ;
	(void) umfpack_di_solve (UMFPACK_A, (int*)A->p, (int *)A->i, (double*)A->x, (double *)x->x, (double *)b->x, Numeric, null, null) ;
	umfpack_di_free_numeric (&Numeric) ;
}

/*
This minimizes the equation 3.24 using equation 3.27 of the papers
http://people.csail.mit.edu/sumner/thesis/Sumner2005MMU.pdf
which minimizes ||Ax - f||^2 using  (At.A)x = At.f
We are using UMFPack to solve the system of linear equation with sparse matrices
The matrices are in CHOLMOD structures
*/
cholmod_dense* MinimizeUsingLeastSquares(cholmod_sparse *A, cholmod_dense *f)
{
	/* PREPARATION */
	cholmod_sparse *At, *AtA;
	cholmod_dense *b = cholmod_zeros(A->ncol, 1, CHOLMOD_REAL, &global);
	double alpha[2] = {1,0};
	double beta[2] = {0,0};

	At = cholmod_transpose(A, 2, &global);
	AtA = cholmod_aat(At, NULL, 0, 1, &global);
    cholmod_sort(AtA, &global);
	cholmod_dense *x = cholmod_allocate_dense(AtA->ncol, 1, AtA->ncol, CHOLMOD_REAL, &global);

	cholmod_sdmult(At, 0, alpha, beta, f, b, &global);

	/*
	FILE* file1 = fopen("DATANEW3.txt","w");
	FILE* file2 = fopen("DATANEW4.txt","w");
	cholmod_write_sparse(file1, AtA, NULL,"",&global);
	cholmod_write_dense(file2, b, "", &global);
	fclose(file1);
	fclose(file2);
	*/
	SolveLinearEquations(AtA, b, x);

	cholmod_free_sparse(&At, &global);
	cholmod_free_sparse(&AtA, &global);
	cholmod_free_dense(&b, &global);
	return x;
}

void ComputeCorrespondance_Phase1()
{
	printf("----PHASE I----\n");
	cholmod_dense * deformation;
	srcWindow->CalculateSurfaceInverseList();
	srcWindow->ComputeTriangleMatrix();
	srcWindow->createFreeVertexIndices();
	cholmod_triplet* M;
	cholmod_dense* C;
	

	int nrow = 9 * (srcWindow->getAdjListSize() + srcWindow->getNTriangles());
    int ncol = 3 * (srcWindow->getNFreeVertices() + srcWindow->getNTriangles());

	M = cholmod_allocate_triplet((size_t)nrow, (size_t)ncol, 0, 0, CHOLMOD_REAL, &global);
	C = cholmod_zeros((size_t)nrow, 1, CHOLMOD_REAL, &global);

	int row = 0;
	row = srcWindow->AddSmoothnessCoefficients(M, C, sqrt(weight_smoothness), row);
	row = srcWindow->AddIdentityCoefficients(M,C, sqrt(weight_identity), row);

	cholmod_sparse* A = cholmod_triplet_to_sparse(M, M->nnz, &global);
	deformation = MinimizeUsingLeastSquares(A, C);
       
	srcWindow->ApplyTransformationToModel(targetWindow, (double *)deformation->x);
	cholmod_free_triplet(&M, &global);
	cholmod_free_dense(&C, &global);
	cholmod_free_sparse(&A, &global);
	cholmod_free_dense(&deformation, &global);
	srcWindow->UpdateMovedConstraints();
	srcWindow->deinitializeFreeVertexIndices();
	//delete(inverseSurfaceMap);
	//delete(triangleMatrices);
	printf("COMPLETE\n");
}

void ComputeCorrespondance_Phase2()
{
	cholmod_dense * deformation;
	cholmod_triplet* M;
	cholmod_dense* C;
	printf("----PHASE II----\n");

	
	for (double weight_closest = weight_closest_point_min;
         weight_closest < weight_closest_point_max;
         weight_closest += weight_closest_point_step)
		 

	//float weight_closest = weight_closest_point_min + (state - 1) * weight_closest_point_step;
    {
		printf("Weight : %f\n", weight_closest);
		srcWindow->CalculateSurfaceInverseList();
		srcWindow->ComputeTriangleMatrix();
		srcWindow->createFreeVertexIndices();
	
		int nrow = 9 * (srcWindow->getAdjListSize() + srcWindow->getNTriangles()) + 3 * srcWindow->getNFreeVertices();
		int ncol = 3 * (srcWindow->getNFreeVertices() + srcWindow->getNTriangles());

		M = cholmod_allocate_triplet((size_t)nrow, (size_t)ncol, 0, 0, CHOLMOD_REAL, &global);
		C = cholmod_zeros((size_t)nrow, 1, CHOLMOD_REAL, &global);

		int ret = -1;
		if(phase2mode == LOADMODE)
		{
			char filename[40];
			sprintf(filename, "CP_%d_%d_%d_%d.txt", weight_closest_point_min, weight_closest_point_step, weight_closest_point_max, weight_closest);
			ret = srcWindow->LoadCPsFromFile(filename);
		}
		if(ret == -1)
			srcWindow->FindClosestPointOfMesh(targetWindow);

		if(phase2mode == SAVEMODE)
		{
			char filename[40];
			sprintf(filename, "CP_%d_%d_%d_%d.txt", weight_closest_point_min, weight_closest_point_step, weight_closest_point_max, weight_closest);
			srcWindow->SaveCPsToFile(filename);
		}

		int row = 0;
		row = srcWindow->AddSmoothnessCoefficients(M, C, sqrt(weight_smoothness), row);
		row = srcWindow->AddIdentityCoefficients(M,C, sqrt(weight_identity), row);
		row = srcWindow->AddClosestPointTerm(M,C, weight_closest, row);

		cholmod_sparse* A = cholmod_triplet_to_sparse(M, M->nnz, &global);
		deformation = MinimizeUsingLeastSquares(A, C);
       
		srcWindow->ApplyTransformationToModel(targetWindow, (double *)deformation->x);
		cholmod_free_triplet(&M, &global);
		cholmod_free_dense(&C, &global);
		cholmod_free_sparse(&A, &global);
		cholmod_free_dense(&deformation, &global);
		srcWindow->UpdateMovedConstraints();
		srcWindow->deinitializeFreeVertexIndices();
	}

	printf("COMPLETE\n");
}



void ComputeAndSaveCorrespondances()
{
	srcWindow->ComputeTriangleCentroid();
	targetWindow->ComputeTriangleCentroid();

	double srcthreshold = srcWindow->CalculateThreshold(); 
	double targetThreshold = targetWindow->CalculateThreshold();
	double threshold = (srcthreshold > targetThreshold)? srcthreshold : targetThreshold;

	int size = srcWindow->getNTriangles();
	int* targetTriangle = new int[size];
	double* targetTriangleDist =  new double[size];

	memset(targetTriangle, 0xFF, sizeof(int)*size);
	memset(targetTriangleDist, 0, sizeof(double)*size);

	srcWindow->ListSrcCorrespondances(targetWindow, threshold, targetTriangle, targetTriangleDist, size);
	targetWindow->ListTgtCorrespondances(srcWindow, threshold, targetTriangle, targetTriangleDist, size);


	FILE* file = fopen("CorrespondanceList.txt", "w");

	for(int i = 0; i < size; i++)
	{
		if(targetTriangle[i] >= 0)	
		{
			fprintf(file, "%d %d\n",i, targetTriangle[i]);
		}
	}
	fclose(file);
	
	delete(targetTriangle);
	delete(targetTriangleDist);
}

void DeformationViewer::ComputeTriangleMatrixForDT()
{
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	Mesh::FVIter fvIter;
	long nfaces = mesh_.n_faces();
	if(triangleMatrices)
		delete(triangleMatrices);
	if(columnVectors)
		delete(columnVectors);
	triangleMatrices = new Matrix9x4[nfaces];
	columnVectors = new Matrix9x1[nfaces];


	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{
		int idx = fIter.handle().idx();
		Matrix3x3* surfaceMatrix = &((Matrix3x3*)inverseSurfaceMap)[idx];
		Mesh::VertexHandle vertices[3];
		int k = 0;
		Matrix9x4 m = {0};
		Matrix9x1 c = {0};

		for(fvIter = mesh_.fv_begin(fIter.handle()); fvIter; ++fvIter, ++k)
		{
			vertices[k] = fvIter.handle();
		}

		for (int i = 0, row = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++, row++)
			{
				c[row] = 0;

				double coefficient = -((*surfaceMatrix)[0][j] + (*surfaceMatrix)[1][j] + (*surfaceMatrix)[2][j]);
				m[row][0] = coefficient;

				/* coefficient of v2: a[0,i]*/
				coefficient = (*surfaceMatrix)[0][j];
				m[row][1] = coefficient;

				/* coefficient of v3: a[1,i] */
				coefficient = (*surfaceMatrix)[1][j];
				m[row][2] = coefficient;

				/* coefficient of v4: a[2,i] */
				coefficient = (*surfaceMatrix)[2][j];
				m[row][3] = coefficient;
			}
		}

		int face_idx = fIter.handle().idx();
		memcpy(&triangleMatrices[face_idx], m,sizeof(Matrix9x4));
		memcpy(&columnVectors[face_idx], c,sizeof(Matrix9x1));
	}
}

void DeformationViewer::ComputeTriangleMatrix()
{
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	Mesh::FVIter fvIter;
	long nfaces = mesh_.n_faces();
	if(triangleMatrices)
		delete(triangleMatrices);
	if(columnVectors)
		delete(columnVectors);
	triangleMatrices = new Matrix9x4[nfaces];
	columnVectors = new Matrix9x1[nfaces];


	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{
		int idx = fIter.handle().idx();
		Matrix3x3* surfaceMatrix = &((Matrix3x3*)inverseSurfaceMap)[idx];
		Mesh::VertexHandle vertices[3];
		int k = 0;
		Matrix9x4 m = {0};
		Matrix9x1 c = {0};

		for(fvIter = mesh_.fv_begin(fIter.handle()); fvIter; ++fvIter, ++k)
		{
			vertices[k] = fvIter.handle();
		}

		for (int i = 0, row = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++, row++)
			{
				c[row] = 0;

				double coefficient = -((*surfaceMatrix)[0][j] + (*surfaceMatrix)[1][j] + (*surfaceMatrix)[2][j]);
				if(!isConstraintVertex(vertices[0]))
					m[row][0] = coefficient;
				else
				{
					float location[3];
					int id = getConstraintIndex(vertices[0]);
					targetWindow->getConstrainedVertexLocation(id, location);
					c[row]   -= (coefficient *	location[i]);
				}

				/* coefficient of v2: a[0,i]*/
				coefficient = (*surfaceMatrix)[0][j];
				if(!isConstraintVertex(vertices[1]))
					m[row][1] = coefficient;
				else
				{
					float location[3];
					int id = getConstraintIndex(vertices[1]);
					targetWindow->getConstrainedVertexLocation(id, location);
					c[row]   -= (coefficient *	location[i]);
				}


				/* coefficient of v3: a[1,i] */
				coefficient = (*surfaceMatrix)[1][j];
				if(!isConstraintVertex(vertices[2]))
					m[row][2] = coefficient;
				else
				{
					float location[3];
					int id = getConstraintIndex(vertices[2]);
					targetWindow->getConstrainedVertexLocation(id, location);
					c[row]   -= (coefficient *	location[i]);
				}


				/* coefficient of v4: a[2,i] */
				coefficient = (*surfaceMatrix)[2][j];
				m[row][3] = coefficient;
			}
		}

		int face_idx = fIter.handle().idx();
		memcpy(&triangleMatrices[face_idx], m,sizeof(Matrix9x4));
		memcpy(&columnVectors[face_idx], c,sizeof(Matrix9x1));
	}
}

void DeformationViewer::CalculateTriangleMatrix(int triangleIdx, double matrixOp[3][3])
{
	Mesh::FVIter fvIter;
	int i = 0;
	Mesh::Point p[3];
	Mesh::FaceHandle fH = mesh_.face_handle(triangleIdx);
	for(fvIter = mesh_.fv_begin(fH); fvIter; ++fvIter, ++i)
	{
		p[i] = mesh_.point(fvIter.handle());
	}

	/* Create the surface of the face that is v2-v1, v3-v1, v4 = (v2-v1)x(v3-v1)*/

	/////////////////////////////////////////////////////////////
	Mesh::Point v2_v1 = p[1] - p[0];
	Mesh::Point v3_v1 = p[2] - p[0];
	Mesh::Point v4 = v2_v1 % v3_v1;
	v4 /= (sqrt(v4.norm()));
	Matrix3x3 m;
	m[0][0] = v2_v1[0];
	m[1][0] = v2_v1[1];
	m[2][0] = v2_v1[2];

	m[0][1] = v3_v1[0];
	m[1][1] = v3_v1[1];
	m[2][1] = v3_v1[2];

	m[0][2] = v4[0];
	m[1][2] = v4[1];
	m[2][2] = v4[2];
	/////////////////////////////////////////////////////////////
		
	memcpy(matrixOp, m,sizeof(Matrix3x3));
}

void DeformationViewer::CalculateSurfaceInverseList()
{
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	Mesh::FVIter fvIter;
	long nfaces = mesh_.n_faces();
	if(inverseSurfaceMap)
		delete(inverseSurfaceMap);
	inverseSurfaceMap = new Matrix3x3[nfaces];

	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{
		int i = 0;
		Mesh::Point p[3];
		for(fvIter = mesh_.fv_begin(fIter.handle()); fvIter; ++fvIter, ++i)
		{
			p[i] = mesh_.point(fvIter.handle());
		}

		/* Create the surface of the face that is v2-v1, v3-v1, v4 = (v2-v1)x(v3-v1)*/

		/////////////////////////////////////////////////////////////
		Mesh::Point v2_v1 = p[1] - p[0];
		Mesh::Point v3_v1 = p[2] - p[0];
		Mesh::Point v4 = v2_v1 % v3_v1;
		v4 /= (sqrt(v4.norm()));
		Matrix3x3 m;
		m[0][0] = v2_v1[0];
		m[1][0] = v2_v1[1];
		m[2][0] = v2_v1[2];

		m[0][1] = v3_v1[0];
		m[1][1] = v3_v1[1];
		m[2][1] = v3_v1[2];

	    m[0][2] = v4[0];
		m[1][2] = v4[1];
		m[2][2] = v4[2];
		/////////////////////////////////////////////////////////////
		
		Matrix3x3 MInverse;
		MatrixInverse(m, MInverse);

		int face_idx = fIter.handle().idx();
		memcpy(&((Matrix3x3*)inverseSurfaceMap)[face_idx], MInverse,sizeof(Matrix3x3));
	}
}

void DeformationViewer::DTransfer()
{
	targetWindow->CalculateSurfaceInverseList();
	targetWindow->ComputeTriangleMatrixForDT();
	
	Mesh::FaceIter fIter, fend = mesh_.faces_end();
	Matrix3x3 identityMatrix = {0};
	identityMatrix[0][0] = identityMatrix[1][1] = identityMatrix[2][2] = 1.0;
	int nTriangles = getNTriangles();
	long nrows = 9 * nTriangles;
	long ncols = 3*(getTotalVertices() + nTriangles);
	srcWindow->CalculateSurfaceInverseList();
	cholmod_triplet* M = cholmod_allocate_triplet((size_t)nrows, (size_t)ncols, 0, 0, CHOLMOD_REAL, &global);
	cholmod_dense *C     = cholmod_zeros((size_t)nrows, (size_t)1, CHOLMOD_REAL, &global);
	createFreeVertexIndices();
	int row = 0;
	for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
	{
		int triangleVertIdx[3][4];
		int triangleMatrixId = fIter.handle().idx();
		
		unordered_map<int,list<int>>::iterator it = corresMap.find(triangleMatrixId);
		
		if(it == corresMap.end())
		{
			GetTriangleVerticesIndex(triangleMatrixId, triangleVertIdx);
			AddTermToMatrixDT(M, triangleMatrixId, row, triangleVertIdx);
			AddColToDenseMatrixDT(C, row, identityMatrix);
			row += 9;
		}
		else
		{
			list<int>  l = it->second;
			for(list<int>::iterator lIt = l.begin(); lIt != l.end(); l++)
			{
				int targetIdx = *lIt;
				GetTriangleVerticesIndex(triangleMatrixId, triangleVertIdx);
				AddTermToMatrixDT(M, triangleMatrixId, row, triangleVertIdx);
				Matrix3x3 triangleMat;
				Matrix3x3 FinalMatJ;
				deformedSrcWindow->CalculateTriangleMatrix(targetIdx,triangleMat);
				MatrixProduct(triangleMat, ((Matrix3x3*)(srcWindow->inverseSurfaceMap))[triangleMatrixId], FinalMatJ);
				AddColToDenseMatrixDT(C, row, FinalMatJ);
				row +=9;
			}
		}
		
	}
	

	cholmod_sparse* A = cholmod_triplet_to_sparse(M, M->nnz, &global);
	cholmod_dense *deformation = MinimizeUsingLeastSquares(A, C);
	ApplyDeformationToModel((double *)deformation->x);
	cholmod_free_triplet(&M, &global);
	cholmod_free_dense(&C, &global);
	cholmod_free_sparse(&A, &global);
	cholmod_free_dense(&deformation, &global);
}

void ComputeCorrespondance()
{
	std::vector<CorrespondanceEntry> corres_list;
	
	/*
	if(state == 0)
	{
		Initialize();
		ComputeCorrespondance_Phase1();
	}
	else
		ComputeCorrespondance_Phase2();

	if(state == 7)
	{
		ComputeAndSaveCorrespondances();
		Deinitialize();
	}
	state++;
	*/
	Initialize();
	ComputeCorrespondance_Phase1();
	ComputeCorrespondance_Phase2();
	ComputeAndSaveCorrespondances();
	Deinitialize();
	srcWindow->writeMesh("Horse_deformed_to_Camel.obj");
}