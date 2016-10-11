//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS SmoothingViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "SmoothingViewer.hh"



//== IMPLEMENTATION ========================================================== 


SmoothingViewer::
SmoothingViewer(const char* _title, int _width, int _height)
  : QualityViewer(_title, _width, _height)
{ 
  mesh_.add_property(vpos_);
}


//-----------------------------------------------------------------------------


void
SmoothingViewer::
keyboard(int key, int x, int y)
{
	switch (toupper(key))
	{
	case 'N':
		{
			std::cout << "10 Laplace-Beltrami smoothing iterations: " << std::flush;
			smooth(10);
			calc_weights();
			calc_mean_curvature();
			calc_uniform_mean_curvature();
			calc_gauss_curvature();
			calc_triangle_quality();
			face_color_coding();

			glutPostRedisplay();
			std::cout << "done\n";
			break;
		}
	case 'U':
		{
			std::cout << "10 uniform smoothing iterations: " << std::flush;
			uniform_smooth(10);
			calc_weights();
			calc_mean_curvature();
			calc_uniform_mean_curvature();
			calc_gauss_curvature();
			calc_triangle_quality();
			face_color_coding();

			glutPostRedisplay();
			std::cout << "done\n";
			break;
		}


	default:
		{
			QualityViewer::keyboard(key, x, y);
			break;
		}
	}
}


//-----------------------------------------------------------------------------


void 
SmoothingViewer::
smooth(unsigned int _iters)
{
	Mesh::VertexIter        v_it, v_end(mesh_.vertices_end());
	Mesh::HalfedgeHandle    h;
	Mesh::EdgeHandle        e;
	Mesh::VertexVertexIter  vv_it;
	Mesh::Point             laplace(0.0, 0.0, 0.0);
	Mesh::Scalar            w, ww;

	// ------------- IMPLEMENT HERE ---------
	// TASK 3.3.b Smoothing using the Laplace-Beltrami.
	// Use eweight_ properties for the individual edge weights
	// and their sum for the normalization term.
	// ------------- IMPLEMENT HERE ---------
	for(uint i = 0; i < _iters; ++i)
	{
		for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			Mesh::Point m(0.0,0.0,0.0);
			ww = 0;
			for( Mesh::VertexOHalfedgeIter vh_it = mesh_.voh_begin(v_it); vh_it; ++vh_it)
			{	
				Mesh::Point n(0.0,0.0,0.0);
				VertexHandle p0 = mesh_.to_vertex_handle(vh_it.handle());
				EdgeHandle eh = mesh_.edge_handle(vh_it.handle());
				n = mesh_.point(p0) - mesh_.point(v_it);
				w = mesh_.property(eweight_,eh);
				n *= w;
				m += n;
				ww += w;
			}
			m /= ww;
			mesh_.property(vpos_, v_it) = m;
		}
		for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			Mesh::Point p = mesh_.point(v_it);
			Mesh::Point shift = mesh_.property(vpos_, v_it);
			mesh_.set_point(v_it.handle(), (p + (shift * 0.5)));
		}
		calc_weights();
		mesh_.update_normals();
	}
}

//-----------------------------------------------------------------------------


void 
SmoothingViewer::
uniform_smooth(unsigned int _iters)
{
	Mesh::VertexIter        v_it, v_end(mesh_.vertices_end());
	Mesh::VertexVertexIter  vv_it;
	Mesh::Point             centroid(0.0, 0.0, 0.0);
	unsigned				  w;

	// ------------- IMPLEMENT HERE ---------
	// TASK 3.1.b Smoothing using the uniform Laplacian approximation
	// ------------- IMPLEMENT HERE ---------

	for(w = 0; w < _iters; ++w)
	{
		for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			Mesh::Point             centroid(0.0, 0.0, 0.0);
			long n_neighbours = 0;
			for(vv_it = mesh_.vv_iter(v_it); vv_it; ++vv_it)
			{
				centroid += mesh_.point(vv_it);
				n_neighbours++;
			}
			centroid /= n_neighbours;
			Mesh::Point p = mesh_.point(v_it);
			Mesh::Point laplacian =  (centroid - p);
			mesh_.property(vpos_, v_it) = laplacian;
		}
		for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			Mesh::Point p = mesh_.point(v_it);
			Mesh::Point shift = mesh_.property(vpos_, v_it);
			mesh_.set_point(v_it.handle(), (p + (shift * 0.5)));
		}
		mesh_.update_normals();
	}
}

//=============================================================================
