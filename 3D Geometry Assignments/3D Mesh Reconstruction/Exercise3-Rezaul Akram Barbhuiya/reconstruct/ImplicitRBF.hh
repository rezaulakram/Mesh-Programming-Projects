//=============================================================================


#ifndef RBF_HH
#define RBF_HH

#define USE_WENDLAND  0
//=============================================================================


#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>
#include <gmm.h>


//=============================================================================


class ImplicitRBF
{
public:

	typedef OpenMesh::Vec3f            Vec3f;
	typedef OpenMesh::Vec3d            Vec3d;
  typedef gmm::dense_matrix<double>  gmmMatrix;
	typedef std::vector<double>        gmmVector;

	
  // fit RBF to given constraints
  ImplicitRBF( const std::vector<Vec3f>& _points, 
							 const std::vector<Vec3f>& _normals );


  // evaluate RBF at position _p
  double operator()(const Vec3f& _p) const;

  

private:

	// evaluate basis function of RBF
	static double kernel(const Vec3d& _c, const Vec3d& _x)
  {
    double r = (_x-_c).norm();
#if USE_WENDLAND
	// Cannot find the reference of sigma in original paper. Setting sigma  = 1 gives a nice output
	double h = 1; // sigma
	double frac = r / h;
	if((1-frac) < 0) 
		return 0;
	else
		return (pow((1 - frac) ,4))*(4 * frac + 1);
#else
    return r*r*r;
#endif
  }  

	// solve linear system _A * _x = _b
	void solve_linear_system( gmmMatrix& _A, 
														gmmVector& _b, 
														gmmVector& _x );
  
private:

  /*
  // Cannot find the reference of sigma in original paper. Setting sigma  = 1 gives a nice output
  static double separationDistance;
  */
  std::vector<Vec3d>   centers_;
  std::vector<double>  weights_;
};


//=============================================================================
#endif // RBF_HH defined
//=============================================================================

