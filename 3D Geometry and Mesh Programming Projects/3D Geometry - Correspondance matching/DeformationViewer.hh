//=============================================================================
//
//  CLASS DeformationViewer
//
//=============================================================================


#ifndef DEFORMATION_VIEWER_HH
#define DEFORMATION_VIEWER_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <unordered_map>
#include "MeshViewer.hh"
#include "cholmod.h"
#include "umfpack.h"
#include <vector>





//== CLASS DEFINITION =========================================================

	      

class DeformationViewer : public MeshViewer
{
public:
   
	/// default constructor
	DeformationViewer(const char* _title, int _width, int _height, bool typeMesh);

	// destructor
	~DeformationViewer();
	
	/// open mesh
	virtual bool open_mesh(const char* _filename);
	void UpdateMovedConstraints();
	void SaveCPsToFile(char* filename);
	int LoadCPsFromFile(char* filename);
	void MouseMovement(int button, int state, int x, int y)
	{
		MeshViewer::mouse(button, state, x, y);
	}


	bool writeMesh(char* filename)
	{
		bool res = OpenMesh::IO::write_mesh(mesh_, filename);
		return res;
	}

	bool insertIntoConstraintList()
	{
		if(!corres_src.is_valid())
			return false;
		constraintList.push_back(corres_src);
		constraintMap.insert(std::make_pair<int, int>(corres_src.idx(), constraintList.size()-1));
		mesh_.property(ConstraintVertex, corres_src) = true;
		corres_src.invalidate();
	}
	
	void insertIntoConstraintList(int id)
	{
		Mesh::VertexHandle v = mesh_.vertex_handle(id);
		constraintList.push_back(v);
		constraintMap.insert(std::make_pair<int, int>(v.idx(), constraintList.size()-1));
		mesh_.property(ConstraintVertex, v) = true;
	}

	int getConstraintIndex(unsigned int vertexHandle)
	{
		Mesh::VertexHandle v = mesh_.vertex_handle(vertexHandle);
		std::unordered_map<int, int>::iterator it = constraintMap.find(v.idx());
		return it->second;
	}

	int getConstraintIndex(Mesh::VertexHandle v)
	{
		std::unordered_map<int, int>::iterator it = constraintMap.find(v.idx());
		return it->second;
	}

	void getConstrainedVertexLocation(unsigned int id, float* location)
	{
		Mesh::VertexHandle v = constraintList[id];
		Mesh::Point p = mesh_.point(v);
		memcpy(location, p.data(), 3* sizeof(float));
	}

	void getVertexLocation(unsigned int id, float* location)
	{
		Mesh::VertexHandle v = mesh_.vertex_handle(id);
		Mesh::Point p = mesh_.point(v);
		memcpy(location, p.data(), 3* sizeof(float));
	}

	double getValueFromDeformationGradient(double* deformation_gradient, int idx, int coordinate)
	{
		return deformation_gradient[GetVertexCoordIndex(idx,coordinate)];
	}

	void CalculateTriangleMatrix(int triangleIdx, double matrix[3][3]);
	void DTransfer();
	void DrawVertex(unsigned int id);
	void ApplyDeformationToModel(double* deformation_gradient);
	void ApplyTransformationToModel(DeformationViewer* tgt_model, double* deformation_gradient);
	void CalculateSurfaceInverseList();
	void ComputeTriangleMatrix();
	void ComputeTriangleMatrixForDT();
	long GetVertexCoordIndex(unsigned int vertexId, unsigned int coord);
	long GetNormalCoordIndex(unsigned int triangleId, unsigned int coord);
	void GetTriangleVerticesIndex(int faceId, int vertexIndexTable[][4]);
	int AddSmoothnessCoefficients(cholmod_triplet* M, cholmod_dense* C,double weight, int row);
	int AddIdentityCoefficients(cholmod_triplet* M, cholmod_dense* C,double weight, int row);
	void FindClosestPointOfMesh(DeformationViewer* target);
	int getAdjListSize();
	void createFreeVertexIndices();
	void deinitializeFreeVertexIndices();
	int getTotalVertices();
	int getNFreeVertices();
	int getNTriangles();
	int AddClosestPointTerm(cholmod_triplet* M, cholmod_dense* C,double weight, int row);
	void ComputeTriangleCentroid();
	std::vector<unsigned int> GetCorrespondingTriangles(Mesh::Point centroid, Mesh::Point normal, double threshold, double* distance);
	void ListCorrespondances(DeformationViewer* target, double threshold, char* filename);
	void ListSrcCorrespondances(DeformationViewer* target, double threshold, int* tgts, double* dist, int n = -1);
	void ListTgtCorrespondances(DeformationViewer* target, double threshold, int* tgts, double* dist, int n = -1);
	
	double CalculateThreshold();

	int findClosestPoint(Mesh::Point location, Mesh::Point Normal)
	{
		Mesh::VertexIter vIter, v_end = mesh_.vertices_end();
		int closest_tgt = -1;
		double min_dist = DBL_MAX;

		for(vIter = mesh_.vertices_begin(); vIter != v_end; ++vIter)
		{
			Mesh::Point targetPoint = mesh_.point(vIter.handle());
			Mesh::Point targetNormal = mesh_.normal(vIter.handle());

			if((Normal | targetNormal) > 0)
			{
				double dist = (targetPoint - location).length();
				if((closest_tgt == -1) || dist < min_dist)
				{
					min_dist = dist;
					closest_tgt = vIter.handle().idx();
				}
			}
		}

		return closest_tgt;
	}


private:

	virtual void keyboard(int key, int x, int y);
	virtual void draw(const std::string& _draw_mode);
	//virtual void motion(int x, int y);
	virtual void mouse(int button, int state, int x, int y);

	
	enum Mode { MOVE, PICK, DRAG } mode_;
	void set_mode( Mode _mode );

	void glText(int x, int y, const std::string& _text);
	
	Mesh::Point& orig_point(Mesh::VertexHandle _vh) 
	{ return mesh_.property(orig_point_, _vh); }

	bool isConstraintVertex(Mesh::VertexHandle _vh) 
	{
		return mesh_.property(ConstraintVertex, _vh);	
	}

	Mesh::VertexHandle findClosestVertex(Mesh::Point p)
	{
		Mesh::VertexIter v_it, v_end;
		v_it = mesh_.vertices_begin();
		v_end = mesh_.vertices_end();
		Mesh::VertexHandle resPoint = v_it.handle();
		float minLength=FLT_MAX, newLength=0;
		for(; v_it != v_end; ++v_it)
		{
			newLength = (mesh_.point(v_it.handle()) - p).length();
			if(minLength > newLength)
			{
				minLength = newLength;
				resPoint = v_it.handle();
			}
		}
		return resPoint;
	}

	friend void ComputeCorrespondance();

private:


	std::vector<Mesh::VertexHandle> constraintList;
	std::unordered_map<int, int> constraintMap;

	Mesh::VertexHandle corres_src;

	OpenMesh::VPropHandleT<Mesh::Point>  orig_point_;
	OpenMesh::VPropHandleT<bool>  ConstraintVertex;
	OpenMesh::VPropHandleT<int>  ClosestTargetVertexIdx;
	OpenMesh::FPropHandleT<Mesh::Point>  Centroid;
	std::vector<Vec3f>  orig_constraints_, moved_constraints_;
	int                 active_sphere_;
	bool sourceMesh;
	Mesh::Vertex v;
	int* freeVertIndices;
public :
	void* inverseSurfaceMap;

};


//=============================================================================
#endif //
//=============================================================================

extern int numPointsSelected;
extern bool isFirstSelected;
extern std::vector<Vec3f> mesh1Points,mesh2Points;

extern DeformationViewer *srcWindow;
extern DeformationViewer *targetWindow;
extern DeformationViewer *deformedSrcWindow;

extern double weight_smoothness;
extern double weight_identity;
extern double weight_closest_point_min;
extern double weight_closest_point_max;
extern double weight_closest_point_step;
extern int mode;
//bool isSelectionDone;