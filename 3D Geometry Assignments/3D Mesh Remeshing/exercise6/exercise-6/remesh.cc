//== INCLUDES =================================================================
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <set>
#include <cfloat>



//== IMPLEMENTATION ===========================================================


typedef OpenMesh::TriMesh_ArrayKernelT<>        Mesh;
Mesh                                            mesh;
OpenMesh::VPropHandleT<Mesh::Point>             update;

//desired edge length
float                                           target_length;

// Curvatures
OpenMesh::VPropHandleT<Mesh::Scalar>            vweight_, vunicurvature_, vcurvature_, vgausscurvature_, vtargetlength_, vnewtargetlength_;
OpenMesh::EPropHandleT<Mesh::Scalar>            etargetlength_;
OpenMesh::EPropHandleT<Mesh::Scalar>            eweight_;

int remeshingType;

//-----------------------------------------------------------------------------


void  split_long_edges();
void  collapse_short_edges();
void  equalize_valences();
void  tangential_relaxation();

void calc_weights();
void calc_mean_curvature();
void calc_uniform_mean_curvature();
void calc_gauss_curvature();
void calc_target_length();


//-----------------------------------------------------------------------------


int main(int argc, char **argv)
{
  if (argc < 4) 
  {
    std::cerr << "Usage: \n" 
	      << argv[0] << " <edge-length>  <input_mesh>  <output_mesh> <Uniform Remeshing = 0|Adaptive Remeshing = 1|Edge Curvature = 2>\n\n";
    exit(1);
  }
	target_length = (float)atof(argv[1]);

	if(argv[4])
		remeshingType = (int)atoi(argv[4]);
	else
		remeshingType = 1;
	
  // add required properties
  mesh.request_vertex_status();
  mesh.request_vertex_normals();
  mesh.request_edge_status();
  mesh.request_face_status();
  mesh.request_face_normals();
  mesh.add_property(update);


  // read mesh
  OpenMesh::IO::read_mesh(mesh, argv[2]);
  std::cout << "#vertices: " << mesh.n_vertices() << std::endl;
	mesh.update_normals();


  // compute curvatures
	mesh.add_property(vcurvature_);
	mesh.add_property(vunicurvature_);
	mesh.add_property(vweight_);
	mesh.add_property(vgausscurvature_);
	mesh.add_property(eweight_);
	mesh.add_property(etargetlength_);
	mesh.add_property(vtargetlength_);
	mesh.add_property(vnewtargetlength_);

	calc_weights();
	calc_mean_curvature();
	calc_uniform_mean_curvature();
	calc_gauss_curvature();
	calc_target_length();

	
  // main remeshing loop
  for (int i=0; i<5; ++i)
	{
		split_long_edges();
		collapse_short_edges();
		equalize_valences();
		tangential_relaxation();
	}

	
	// write mesh
  OpenMesh::IO::write_mesh(mesh, argv[3]);
}


//-----------------------------------------------------------------------------


void  split_long_edges()
{
	Mesh::EIter     e_it, e_end;
	Mesh::VHandle   v0, v1, vh, v2, v3;
	Mesh::EHandle   eh, e0, e1;
	Mesh::FHandle   f0, f1, f2, f3;
	bool            finished;
	int             i;
	float           desired_length;


	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;
		
		for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
		{
			// Exercise 6.1 ----------------------------------------------
			// INSERT CODE:
			//  1) compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
			//  If the edge is longer than 4/3 * desired length
			//  2) add the midpoint to the mesh
			//  3) set the interpolated normal and interpolated vtargetlength_ property to the vertex
			//  4) split the edge with this vertex (use openMesh function split)
			// Leave the loop running until no splits are done (use the finished variable)
			// -----------------------------------------------------------

			eh = e_it.handle();
			v0 = mesh.to_vertex_handle(mesh.halfedge_handle(eh,0));
			v1 = mesh.from_vertex_handle(mesh.halfedge_handle(eh,0));
			v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh,0)));
			v3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh,1)));


			desired_length = (mesh.property(vtargetlength_, v0) + mesh.property(vtargetlength_, v1)) / 2;

			float length =  mesh.calc_edge_length(eh);
			desired_length = 2.00;
			if(length > (desired_length * 4.0 / 3))
			{
				vh = mesh.add_vertex((mesh.point(v0) + mesh.point(v1))/2);
				mesh.set_normal(vh, (mesh.normal(v0) + mesh.normal(v1))/ 2);
				mesh.property(vtargetlength_, vh) = desired_length;
				
				mesh.split(eh, vh);

				for(Mesh::VertexFaceIter vf_it = mesh.vf_begin(vh); vf_it; ++vf_it) 
				{
					mesh.update_normal(vf_it.handle());
				}
				
				finished = false;
			}
		}
	}
}


//-----------------------------------------------------------------------------


void  collapse_short_edges()
{
	Mesh::EIter     e_it, e_end;
	Mesh::CVVIter   vv_it;
	Mesh::VHandle   v0, v1;
	Mesh::HHandle   h0, h1, h01, h10;
	bool            finished, /*skip,*/ b0, b1/*, l0, l1, f0, f1*/;
	int             i;
	bool            hcol01, hcol10;

	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
		{
			if (!mesh.status(e_it).deleted()) // might already be deleted
			{

				// Exercise 6.2 ----------------------------------------------
				// INSERT CODE:
				//  1) Compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
				// If the edge is shorter than 4/5 of the desired length
				//  2) Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse. Otherwise
				//  3) Check if halfedges collapsible
				//  4) Select the halfedge to be collapsed if at least one halfedge can be collapsed
				//  5) Collapse the halfedge
				// Leave the loop running until no collapse has been done (use the finished variable)
				// -----------------------------------------------------------

				Mesh::EHandle eh = e_it.handle();
				h0 = mesh.halfedge_handle(eh,0);
				h1 = mesh.halfedge_handle(eh,1);
				v0 = mesh.from_vertex_handle(h0);
				v1 = mesh.to_vertex_handle(h0);
				Mesh::Point p0 = mesh.point(v0);
				Mesh::Point p1 = mesh.point(v1);


				double desired_length;

				desired_length = (mesh.property(vtargetlength_, v0) + mesh.property(vtargetlength_, v1)) / 2;

				float length =  mesh.calc_edge_length(eh);
				if(length < (desired_length * 4.0 / 5) && (mesh.is_boundary(v0) == mesh.is_boundary(v1)))
				{
					bool collapseb0 = true, collapseb1 = true;
					
					for(vv_it = mesh.vv_begin(v0); vv_it; ++vv_it)
					{
						double maxLength = (mesh.property(vtargetlength_,vv_it) + mesh.property(vtargetlength_,v1)) * 2 / 3;
						Mesh::Point p = mesh.point(vv_it.handle());
						if((p - p1).length() > maxLength)
						{
							collapseb0 = false;
							break;
						}
					}

					for(vv_it = mesh.vv_begin(v1); vv_it; ++vv_it)
					{
						double maxLength = (mesh.property(vtargetlength_,vv_it) + mesh.property(vtargetlength_,v0)) * 2 / 3;
						Mesh::Point p = mesh.point(vv_it.handle());
						if((p - p0).length() > maxLength)
						{
							collapseb1 = false;
							break;
						}
					}

					b0 = mesh.is_collapse_ok(h0) && collapseb0;
					b1 = mesh.is_collapse_ok(h1) && collapseb1;

					if(b0 && b1)
					{
						if(mesh.valence(v0) < mesh.valence(v1))
						{
							mesh.collapse(h0);
						}
						else
						{
							mesh.collapse(h1);
						}
						finished = false;
					}
					else if(b0)
					{
						mesh.collapse(h0);
						finished = false;
					}
					else if(b1)
					{
						mesh.collapse(h1);
						finished = false;
					}
				}
			}
		}
	}

	mesh.garbage_collection();

	if (i==100) std::cerr << "collapse break\n";
}


//-----------------------------------------------------------------------------


void  equalize_valences()
{
  Mesh::EIter     e_it, e_end;
  Mesh::VHandle   v0, v1, v2, v3;
  Mesh::HHandle   hh, h0;
  int             val0, val1, val2, val3;
  int             val_opt0, val_opt1, val_opt2, val_opt3;
  int             ve0, ve1, ve2, ve3, ve_before, ve_after;
  bool            finished;
  int             i;
  
	


  // flip all edges
  for (finished=false, i=0; !finished && i<100; ++i)
  {
	  finished = true;

	  for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
	  {
		  if (!mesh.is_boundary(e_it))
		  {

				// Exercise 6.3 ----------------------------------------------
				// INSERT CODE:
				//  1) Extract valences of the four vertices involved to an eventual flip.
				//  2) Compute the sum of the squared valence deviances before flip
				//  3) Compute the sum of the squared valence deviances after and eventual flip
				//  4) If valence deviance is decreased and flip is possible, flip the vertex
				// Leave the loop running until no collapse has been done (use the finished variable)
				// -----------------------------------------------------------
				Mesh::EHandle eh = e_it.handle();
				hh = mesh.halfedge_handle(eh,0);
				h0 = mesh.halfedge_handle(eh,1);
				
				Mesh::EHandle e1 = mesh.edge_handle(mesh.next_halfedge_handle(hh));
				Mesh::EHandle e2 = mesh.edge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(hh)));
				Mesh::EHandle e3 = mesh.edge_handle(mesh.next_halfedge_handle(h0));
				Mesh::EHandle e4 = mesh.edge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(h0)));

				Mesh::Normal e1_vec, e2_vec, e3_vec, e4_vec;

				mesh.calc_edge_vector(e1, e1_vec);
				mesh.calc_edge_vector(e2, e2_vec);
				mesh.calc_edge_vector(e3, e3_vec);
				mesh.calc_edge_vector(e4, e4_vec);

				v0 = mesh.from_vertex_handle(hh);
				v1 = mesh.to_vertex_handle(hh);
				v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
				v3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(hh)));

				val0 =  mesh.valence(v0);
				val1 =  mesh.valence(v1);
				val2 =  mesh.valence(v2);
				val3 =  mesh.valence(v3);

				ve0 = val0 - 1;
				ve1 = val1 - 1;
				ve2 = val2 + 1;
				ve3 = val3 + 1;

				val_opt0 = (mesh.is_boundary(v0))? 4 :6;
				val_opt1 = (mesh.is_boundary(v1))? 4 :6;
				val_opt2 = (mesh.is_boundary(v2))? 4 :6;
				val_opt3 = (mesh.is_boundary(v3))? 4 :6;

				ve_before = ((val0 - val_opt0) * (val0 - val_opt0)) +
							((val1 - val_opt1) * (val1 - val_opt1)) +
							((val2 - val_opt2) * (val2 - val_opt2)) +
							((val3 - val_opt3) * (val3 - val_opt3));

				ve_after  = ((ve0 - val_opt0) * (ve0 - val_opt0)) +
							((ve1 - val_opt1) * (ve1 - val_opt1)) +
							((ve2 - val_opt2) * (ve2 - val_opt2)) +
							((ve3 - val_opt3) * (ve3 - val_opt3));

				float cosAngle1 = (e1_vec | e4_vec)/(e1_vec.length() * e4_vec.length());
				float cosAngle2 = (e2_vec | e3_vec)/(e2_vec.length() * e3_vec.length());


				if((ve_after < ve_before) && mesh.is_flip_ok(eh) &&  cosAngle1 < 0.9999 && cosAngle1 > -0.9999 && cosAngle2 < 0.9999 && cosAngle2 > -0.9999)
				{
					mesh.flip(eh);
					finished = false;
				}
		  }
	  }
  }

  if (i==100) std::cerr << "flip break\n";
}


//-----------------------------------------------------------------------------


void  tangential_relaxation()
{
	Mesh::VIter     v_it, v_end(mesh.vertices_end());
	Mesh::CVVIter   vv_it;
	Mesh::Scalar    valence;
	Mesh::Point     u, n;

	
  // smooth
  for (int iters=0; iters<10; ++iters)
  {
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
    {
      if (!mesh.is_boundary(v_it))
      {
				// Exercise 6.4 ----------------------------------------------
				// INSERT CODE:
				//  1) Compute uniform laplacian approximation vector
				//  2) Compute the tangential component of the laplacian vector
				//  3) Store smoothed vertex location in the update vertex property.
				//     (you don't have to use 1/2 attenuation in this case, it's fine without attenuation)
				// -----------------------------------------------------------

			valence = 0;
			Mesh::Point m(0.0,0.0,0.0);
			for(vv_it = mesh.vv_iter(v_it); vv_it; ++vv_it)
			{
				m += mesh.point(vv_it);
				valence++;
			}
			m /= valence;
			Mesh::Point laplace = m - mesh.point(v_it);

			Mesh::Normal normal = mesh.normal(v_it);
			normal.normalize();

			Mesh::Point normalizedComp =  (laplace | normal) * normal;
			Mesh::Point tangentialComp =  laplace - normalizedComp;

			mesh.property(update, v_it) = tangentialComp;
      }
    }
		
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      if (!mesh.is_boundary(v_it))
				mesh.point(v_it) += mesh.property(update, v_it);
  }
}


//-----------------------------------------------------------------------------

void calc_target_length()
{
	Mesh::VertexIter        v_it, v_end;
	Mesh::HalfedgeHandle    h;
	Mesh::EdgeHandle        e;
	Mesh::VertexVertexIter  vv_it;
	float                   length;
	float                   meanlength;
	float                   H;
	float                   K;
	int                     i;

		// Exercise 6.5 ----------------------------------------------
		// INSERT CODE:
		//  1) Compute the maximal curvature at each vertex (use the precomputed mean (property: vcurvature_) and gaussian curvature (property: vgausscurvature_))
		//  2) Calculate the desired edge length as the target_length divided by the maximal curvature at each vertex, and assign it to the property vtargetlength_
		//  2) Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
		//  3) Rescale the property vtargetlength_ such that it's mean equals the user specified target_length
		// -----------------------------------------------------------

	if(remeshingType == 0)
	{
		for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
		{
			length = target_length;
			mesh.property(vtargetlength_, v_it) = length;
		}
	}
	else 
	{		
		if(remeshingType == 1)
		{
			//calculate desired length 
			for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
			{
				length = target_length;
				if (!mesh.is_boundary(v_it.handle()))
				{
					//1)

					H = mesh.property(vcurvature_, v_it);
					K = mesh.property(vgausscurvature_, v_it);

					float kmax =  H + sqrt(H*H - K);

					length /= kmax;
				}
				mesh.property(vtargetlength_, v_it) = length;
			}
		}
		else
		{
			for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
			{
				length = target_length;
				if (!mesh.is_boundary(v_it.handle()))
				{
					Mesh::Point p = mesh.point(v_it.handle());
					Mesh::Normal N =  mesh.normal(v_it.handle());
					float maxCurvature = FLT_MIN;
					for(Mesh::VertexVertexIter vv_it = mesh.vv_begin(v_it); vv_it; ++vv_it)
					{
						Mesh::Point pi = mesh.point(vv_it.handle());
						Mesh::Normal Ni =  mesh.normal(vv_it.handle());
			
						Mesh::Point deltaP =  pi - p;
						Mesh::Point deltaN =  Ni -N;

						Mesh::Scalar curvature =  -(deltaP | deltaN)/(deltaP | deltaP);
						if(abs(curvature) > maxCurvature) maxCurvature = abs(curvature);

					}
					length /= maxCurvature;
				}
				mesh.property(vtargetlength_, v_it) = length;
			}
		}
		
		//smooth desired length a few times
		for( i=0;i<5;i++)
		{
			//2)
			for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
			{
				meanlength = mesh.property(vtargetlength_, v_it);
				i = 1;
				for(vv_it = mesh.vv_iter(v_it);vv_it;++vv_it, ++i)
				{
					meanlength += mesh.property(vtargetlength_, vv_it);
				}
				meanlength /=  i;

				mesh.property(vnewtargetlength_, v_it) = meanlength;
			}

			for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
			{
				mesh.property(vtargetlength_, v_it) =  mesh.property(vnewtargetlength_, v_it);
			}
		}

		//rescale desired length:
		//3)

		meanlength = 0;
		i = 0;
		for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it, ++i)
		{
			meanlength += mesh.property(vtargetlength_, v_it);
		}

		meanlength /= i;

		float scale = target_length / meanlength;

		for (v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
		{
			length = mesh.property(vtargetlength_, v_it);
			mesh.property(vtargetlength_, v_it) = length * scale;
			//mesh.property(vtargetlength_, v_it) = target_length;
		}
	}
}


//-----------------------------------------------------------------------------

void calc_weights()
{
	Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	Mesh::EdgeIter          e_it, e_end(mesh.edges_end());
	Mesh::VertexFaceIter    vf_it;
	Mesh::FaceVertexIter    fv_it;
	Mesh::HalfedgeHandle    h0, h1, h2;
	Mesh::VertexHandle      v0, v1;
	Mesh::Point             p0, p1, p2, d0, d1;
	Mesh::Scalar            w, area;



	for (e_it=mesh.edges_begin(); e_it!=e_end; ++e_it)
	{
		w  = 0.0;

		h0 = mesh.halfedge_handle(e_it.handle(), 0);
		v0 = mesh.to_vertex_handle(h0);
		p0 = mesh.point(v0);

		h1 = mesh.halfedge_handle(e_it.handle(), 1);
		v1 = mesh.to_vertex_handle(h1);
		p1 = mesh.point(v1);

		h2 = mesh.next_halfedge_handle(h0);
		p2 = mesh.point(mesh.to_vertex_handle(h2));
		d0 = (p0 - p2).normalize();
		d1 = (p1 - p2).normalize();
		w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1)))));

		h2 = mesh.next_halfedge_handle(h1);
		p2 = mesh.point(mesh.to_vertex_handle(h2));
		d0 = (p0 - p2).normalize();
		d1 = (p1 - p2).normalize();
		w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1)))));

		w = std::max(0.0f, w);
		mesh.property(eweight_,e_it) = w;
	}


	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	{
		area = 0.0;

		for (vf_it=mesh.vf_iter(v_it); vf_it; ++vf_it)
		{
			fv_it = mesh.fv_iter(vf_it);
 
			const Mesh::Point& P = mesh.point(fv_it);  ++fv_it;
			const Mesh::Point& Q = mesh.point(fv_it);  ++fv_it;
			const Mesh::Point& R = mesh.point(fv_it);

			area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;
		}

		mesh.property(vweight_,v_it) = 1.0 / (2.0 * area);
	}
}


//-----------------------------------------------------------------------------


void calc_mean_curvature()
{
	Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	Mesh::HalfedgeHandle    h;
	Mesh::EdgeHandle        e;
	Mesh::VertexVertexIter  vv_it;
	Mesh::Point             laplace(0.0, 0.0, 0.0);


	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	{
		mesh.property(vcurvature_, v_it) = 0.0;
		laplace = Mesh::Point(0,0,0);

		if (!mesh.is_boundary(v_it.handle()))
		{
			for (vv_it=mesh.vv_iter(v_it); vv_it; ++vv_it)
			{
				h = vv_it.current_halfedge_handle();
				e = mesh.edge_handle(h);

				laplace += mesh.property(eweight_,e) * (mesh.point(vv_it) - mesh.point(v_it));
			}
			laplace *= mesh.property(vweight_, v_it);

			mesh.property(vcurvature_, v_it) = laplace.norm();
		}
	}
}

void calc_uniform_mean_curvature()
{
	Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	Mesh::VertexVertexIter  vv_it;
	Mesh::Point             laplace(0.0, 0.0, 0.0);


	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	{
		mesh.property(vunicurvature_, v_it) = 0.0;
		laplace = Mesh::Point(0,0,0);

		if (!mesh.is_boundary(v_it.handle()))
		{
			for (vv_it=mesh.vv_iter(v_it); vv_it; ++vv_it)
			{
				laplace += (mesh.point(vv_it) - mesh.point(v_it));
			}
			laplace *= mesh.property(vweight_, v_it);

			mesh.property(vunicurvature_, v_it) = laplace.norm();
		}
	}
}

void calc_gauss_curvature()
{
	Mesh::VertexIter        v_it, v_end(mesh.vertices_end());
	Mesh::VertexVertexIter  vv_it, vv_it2;
	Mesh::Point             d0, d1;
	Mesh::Scalar            angles, cos_angle;
	Mesh::Scalar            lb(-1.0), ub(1.0);


	// compute for all non-boundary vertices
	for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
	{
		if (!mesh.is_boundary(v_it))
		{
			angles = 0.0;

			for (vv_it=mesh.vv_iter(v_it); vv_it; ++vv_it)
			{
				vv_it2 = vv_it; ++vv_it2;
				d0 = (mesh.point(vv_it)  - mesh.point(v_it)).normalize();
				d1 = (mesh.point(vv_it2) - mesh.point(v_it)).normalize();
				cos_angle = std::max(lb, std::min(ub, (d0 | d1)));
				angles += acos(cos_angle);
			}

			mesh.property(vgausscurvature_, v_it) =
				(2*M_PI-angles) * 2.0f * mesh.property(vweight_,v_it);
		}
	}

}
