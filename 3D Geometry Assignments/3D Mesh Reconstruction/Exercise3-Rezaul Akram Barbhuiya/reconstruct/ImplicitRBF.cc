//=============================================================================


#include "ImplicitRBF.hh"

//double ImplicitRBF::separationDistance = 0;

//== IMPLEMENTATION ==========================================================


ImplicitRBF::
ImplicitRBF( const std::vector<Vec3f>& _points, 
						 const std::vector<Vec3f>& _normals )
{
	//////////////////////////////////////////////////////////////////////
	// INSERT CODE:
	// 1) collect constraints (on-surface and off-surface)
	// 2) setup matrix
	// 3) solve linear system for weights_
#if USE_WENDLAND
	
	/*
	double minDist = (_points[0] - _points[1]) .length();
	for(int i = 0; i < _points.size(); i++ )
	{	
		double tDist = 0;
		for(int j = i+1; j < _points.size(); j++ )
		{			
			tDist = (_points[i] - _points[j]) .length();
			if(tDist < minDist)
				minDist = tDist;
		}
	}
	separationDistance = minDist;
	*/
#endif

	Vec3f bb_min( _points[0]), bb_max( _points[0]);
	for (unsigned int i=1; i<_points.size(); ++i)
	{
		bb_min.minimize( _points[i] );
		bb_max.maximize( _points[i] );
	}
	double epsilon = 0.01 * (bb_max - bb_min).length(); 

	int nPoints = _points.size();
	gmmMatrix A(2*nPoints, 2*nPoints);
	gmmVector dist(2*nPoints);
	centers_.resize(2 * nPoints);
	weights_.resize(2 * nPoints);

	for(int i = 0; i < nPoints; i++)
	{
		Vec3d p(_points[i]);
		Vec3d dp(_points[i] + (_normals[i] * epsilon));
		centers_[i] = p;
		centers_[i+nPoints] = dp;
	}

	for(int i = 0; i < nPoints; i++)
	{
		std::vector<Vec3d>::const_iterator c_it;
		int j = 0;
		for(c_it = centers_.begin(); c_it != centers_.end(); ++c_it, ++j)
		{	
			Vec3d p(_points[i]);
			Vec3d dp(_points[i] + (_normals[i] * epsilon));
			A(i,j) = kernel(*c_it, p);
			A(i+nPoints,j) = kernel(*c_it, dp);
		}
		dist[i] = 0;
		dist[i+nPoints] = epsilon;
	}

	solve_linear_system(A, dist, weights_);

	//////////////////////////////////////////////////////////////////////
}


//-----------------------------------------------------------------------------


void
ImplicitRBF::solve_linear_system( gmmMatrix& _M, 
																	gmmVector& _b, 
																	gmmVector& _x)
{
	// solve linear system by gmm's LU factorization
	unsigned int N = _b.size();
	_x.resize(N);
  std::vector< size_t >  ipvt(N);
  gmm::lu_factor( _M, ipvt );
  gmm::lu_solve( _M, ipvt, _x, _b );
}


//-----------------------------------------------------------------------------


double 
ImplicitRBF::operator()(const Vec3f& _p) const
{
  std::vector<Vec3d>::const_iterator  
    c_it(centers_.begin()),
    c_end(centers_.end());

  std::vector<double>::const_iterator   
    w_it(weights_.begin());

	const Vec3d p(_p);
  double f(0);

  for (; c_it!=c_end; ++c_it, ++w_it)
    f += *w_it * kernel(*c_it, p);

  return f;
}


//=============================================================================
