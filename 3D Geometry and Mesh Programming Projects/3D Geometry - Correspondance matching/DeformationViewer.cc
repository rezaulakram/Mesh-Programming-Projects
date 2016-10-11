//=============================================================================
//
//  CLASS DeformationViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "DeformationViewer.hh"
#include "deformationTransfer.hh"
#include "correspondance.hh"
#include <vector>
#include <float.h>




//== IMPLEMENTATION ========================================================== 


DeformationViewer::
DeformationViewer(const char* _title, int _width, int _height, bool typeMesh)
: MeshViewer(_title, _width, _height)
{ 
  	mesh_.add_property(orig_point_);
	mesh_.add_property(ConstraintVertex);
	mesh_.add_property(ClosestTargetVertexIdx);
	mesh_.add_property(Centroid);
  	init();
	set_mode(MOVE);
	sourceMesh = typeMesh;
}


//-----------------------------------------------------------------------------


DeformationViewer::
~DeformationViewer()
{
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
keyboard(int key, int x, int y)
{
	switch (key)
	{
    	case ' ':
    	{
      		std::cout << "Printing hard constraints..." << std::endl;
      		//glutPostRedisplay();
			for(int i=0; i < numPointsSelected; i++)
			{
				std::cout<<mesh1Points[i]<<"\t"<<mesh2Points[i]<<std::endl;
			}
      		std::cout << "done\n";
      		break;
    	}

		case 'l':
		{
			phase2mode = LOADMODE;
			break;	
		}
		case 'k':
		{
			phase2mode = SAVEMODE;
			break;	
		}

		case 'n':
		{
			phase2mode = NORMALMODE;
			break;	
		}
			
		case 'm':
		{
			set_mode(MOVE);
			glutPostRedisplay();
			break;
		}
			
		case 'p':
		{
			set_mode(PICK);
			glutPostRedisplay();
			break;
		}
		case 'c':
		{
			if(mode == 0)
				ComputeCorrespondance();
			else 
				PerformDeformationTransfer();
			break;
		}
			
		/*case 'd':
		{
			set_mode(DRAG);
			glutPostRedisplay();
			break;
		}*/
			
		case 'r': // reset
		{
			set_mode(MOVE);

			Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
				mesh_.set_point(v_it, orig_point(v_it));

			mesh_.update_normals();

			moved_constraints_ = orig_constraints_;
			
			srcWindow->constraintList.empty();
			targetWindow->constraintList.empty();
			
			glutPostRedisplay();
			break;
		}
		case 's': // reset
		{
			FILE* file = fopen("correspondance.cor","w");
			int length_ = srcWindow->constraintList.size();
			fprintf(file,"%d\n",length_);
			for(int i = 0; i < length_; i++)
			{
				fprintf(file,"%d, %d\n",srcWindow->constraintList[i], targetWindow->constraintList[i]);
			}
			fclose(file);
			break;
		}
			
		default:
		{
			MeshViewer::keyboard(key, x, y);
			break;
		}
	}
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
set_mode( Mode _mode )
{
	switch(mode_ = _mode)
	{
		case MOVE:
			glutSetCursor( GLUT_CURSOR_LEFT_ARROW );
			break;

		case PICK:
			glutSetCursor( GLUT_CURSOR_CROSSHAIR );
			break;

		/*case DRAG:
			glutSetCursor( GLUT_CURSOR_INFO );
			break;*/
	}	
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
glText(int x, int y, const std::string& _text)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// set raster pos
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, (GLfloat) viewport[2], 0.0, (GLfloat) viewport[3]);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);

	// draw characters
	std::string::const_iterator s_it(_text.begin()), s_end(_text.end());
	for (; s_it!=s_end; ++s_it)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *s_it);
	
	// restore matrices
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}


//-----------------------------------------------------------------------------


bool
DeformationViewer::
open_mesh(const char* _filename1)
{
	// load mesh
	bool isValidFile = MeshViewer::open_mesh(_filename1);

	Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
	if (isValidFile)
	{
		// store original vertex positions
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
		{
			orig_point(v_it) = mesh_.point(v_it);
			mesh_.property(ConstraintVertex, v_it) = false;
		}
		//return true;
	}
  return false;
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
draw(const std::string& _draw_mode)
{
	// draw mesh
  	MeshViewer::draw(_draw_mode);
	
	
	// draw spheres
	GLfloat mat_sphere[4] = {1.0, 0.0, 0.0, 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_sphere);
	
	for (int i=0; i<moved_constraints_.size(); ++i)
	{
		glPushMatrix();
		glTranslatef( moved_constraints_[i][0],
					  moved_constraints_[i][1],
					  moved_constraints_[i][2] );
		glutSolidSphere(0.007*radius_, 20, 20);
		glPopMatrix();
	}

	/*glBegin( GL_LINES );
	for (int i=0; i<moved_constraints_.size(); ++i)
	{
		GL::glVertex( orig_constraints_[i] );
		GL::glVertex( moved_constraints_[i] );
	}
	glEnd();
	*/
	GLfloat mat_mesh[4] = {0.4, 0.4, 0.4, 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_mesh);

	
	// draw text
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	switch (mode_)
	{
		case MOVE:
			glText(10, 10, "Move");
			break;

		case PICK:
			glText(10, 10, "Pick");
			break;

		/*case DRAG:
			glText(10, 10, "Drag");
			break;*/
	}

	switch (phase2mode)
	{
		case NORMALMODE:
			glText(100, 10, "Normal Mode");
			break;
		case LOADMODE:
			glText(100, 10, "Load Mode");
			break;
		case SAVEMODE:
			glText(100, 10, "Save Mode");
			break;
	
	}
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
mouse(int button, int state, int x, int y)
{
	switch (mode_)
	{
		// move the mesh
		case MOVE:
		{
			MeshViewer::mouse(button, state, x, y);
			break;
		}
		
		// create a new handle point 
		case PICK:
		{
			if (state == GLUT_DOWN)
			{	
				GLdouble  modelview[16], projection[16];
				GLint     viewport[4];
				glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
				glGetDoublev(GL_PROJECTION_MATRIX, projection);
				glGetIntegerv(GL_VIEWPORT, viewport);
				
				
				// read depth buffer value at (x, y_new)
				float  z;
				double y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
				glReadPixels(x, y_new, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );
				
				
				// reverse projection to get 3D point
				double pos[3];
				gluUnProject(x, y_new, z, 
							 modelview, 
							 projection, 
							 viewport, 
							 &pos[0], &pos[1], &pos[2]);
				
				if (z != 1.0f)
				{
					if(sourceMesh)
					{
						if(isFirstSelected)
							std::cout << "First select the point in the second mesh" << std::endl;
						else
						{
							Mesh::Point p((float)pos[0],(float)pos[1],(float)pos[2]);
							orig_constraints_.push_back(p);
							Mesh::VertexHandle v = findClosestVertex(p);
							moved_constraints_.push_back(mesh_.point(v));
							corres_src = v;
							active_sphere_ = moved_constraints_.size()-1;
							std::cout << "Source Mesh : Create new handle sphere at " << p << std::endl;
							isFirstSelected = true;
							numPointsSelected++;
							mesh1Points.push_back(mesh_.point(v));
						}
					}
					else
					{
						if(!isFirstSelected)
							std::cout << "First select the point in the first mesh" << std::endl;
						else
						{
							Mesh::Point p((float)pos[0],(float)pos[1],(float)pos[2]);
							orig_constraints_.push_back(p);
							Mesh::VertexHandle v = findClosestVertex(p);
							moved_constraints_.push_back(mesh_.point(v));
							if (v.is_valid())
							{
								if(!(srcWindow->insertIntoConstraintList()))
									printf("error");
								corres_src = v;
								insertIntoConstraintList();
								
							}
							corres_src.reset();
							active_sphere_ = moved_constraints_.size()-1;
							
							std::cout << "Target Mesh : Create new handle sphere at " << p << std::endl;
							isFirstSelected = false;
							mesh2Points.push_back(mesh_.point(v));
						}
					}
				}
				else
					std::cout << "Cannot create new sphere handle" << std::endl;
				glutPostRedisplay();
			}

			break;
		}
	}
}

void DeformationViewer::UpdateMovedConstraints()
{
	moved_constraints_.erase(moved_constraints_.begin(), moved_constraints_.end());

	int length = constraintList.size();
	for(int i = 0; i < length; i++)
	{
		Mesh::VHandle vh = constraintList[i];
		Mesh::Point p = mesh_.point(vh);
		moved_constraints_.push_back(p);
	}
}

void 
DeformationViewer::DrawVertex(unsigned int id)
{
	Mesh::Point p =  mesh_.point(mesh_.vertex_handle(id));
	moved_constraints_.push_back(p);
	orig_constraints_.push_back(p);
}

void DeformationViewer::ApplyTransformationToModel(DeformationViewer* tgt_model,
								double* deformation_gradient)
{
	Mesh::VertexIter viter, vend = mesh_.vertices_end();
	Mesh::VertexHandle vh;

	for(viter = mesh_.vertices_begin(); viter != vend; ++viter)
	{
		float loc[3];
		Mesh::Point prev;
		Mesh::Point newL;
		vh = viter.handle();
		if(isConstraintVertex(vh))
		{
			
			int id = getConstraintIndex(vh);
			targetWindow->getConstrainedVertexLocation(id, loc);
			Mesh::Point p(loc);
			mesh_.point(vh) = p;
			
		}
		else
		{
			loc[0] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 0);
			loc[1] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 1);
			loc[2] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 2);
			Mesh::Point p(loc);
			prev = mesh_.point(vh);
			mesh_.point(vh) = p;
			newL = mesh_.point(vh);
			mesh_.point(vh) = p;
		}

	}
	mesh_.update_normals();
}


void DeformationViewer::ApplyDeformationToModel(double* deformation_gradient)
{
	Mesh::VertexIter viter, vend = mesh_.vertices_end();
	Mesh::VertexHandle vh;

	for(viter = mesh_.vertices_begin(); viter != vend; ++viter)
	{
		float loc[3];
		Mesh::Point prev;
		Mesh::Point newL;
		vh = viter.handle();
		loc[0] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 0);
		loc[1] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 1);
		loc[2] = (float)getValueFromDeformationGradient(deformation_gradient, vh.idx(), 2);
		Mesh::Point p(loc);
		prev = mesh_.point(vh);
		mesh_.point(vh) = p;
		newL = mesh_.point(vh);
		mesh_.point(vh) = p;
	}
	mesh_.update_normals();
}

int DeformationViewer::getNTriangles()
	{
		return mesh_.n_faces();
	}

	int DeformationViewer::getNFreeVertices()
	{
		return (mesh_.n_vertices() - constraintList.size());
	}

	int DeformationViewer::getTotalVertices()
	{
		return (mesh_.n_vertices());
	}

	void DeformationViewer::deinitializeFreeVertexIndices()
	{
		delete(freeVertIndices);
	}

	void DeformationViewer::createFreeVertexIndices()
	{
		freeVertIndices = new int[getTotalVertices()];
		Mesh::VertexIter vIter, vEnd = mesh_.vertices_end();
		int k = 0;
		int index = 0;
		for(vIter = mesh_.vertices_begin(); vIter != vEnd; ++vIter, k++)
		{
			if((mode == 1) || !isConstraintVertex(vIter.handle()))
			{
				freeVertIndices[k] = index;
				index++;
			}
			else
				freeVertIndices[k] = -1;
		}
	}

	int DeformationViewer::getAdjListSize()
	{
		Mesh::FaceIter fIter, fend = mesh_.faces_end();
		Mesh::FFIter ffIter;

		int nadj = 0;

		for(fIter = mesh_.faces_begin(); fIter != fend; ++fIter)
		{
			for(ffIter = mesh_.ff_begin(fIter);ffIter; ++ffIter)
			{
				nadj++;
			}
		}
		return nadj;
	}

	
	

	void DeformationViewer::FindClosestPointOfMesh(DeformationViewer* target)
	{

		Mesh::VertexIter vIter, v_end = mesh_.vertices_end();
		for(vIter = mesh_.vertices_begin(); vIter != v_end; ++vIter)
		{
			Mesh::Point srcPoint = mesh_.point(vIter.handle());
			Mesh::Point srcNormal = mesh_.normal(vIter.handle());

			int id = target->findClosestPoint(srcPoint, srcNormal);

			mesh_.property(ClosestTargetVertexIdx, vIter) = id;
		}
	}

	void DeformationViewer::ComputeTriangleCentroid()
	{
		Mesh::FIter fIter, f_end = mesh_.faces_end();
		Mesh::FVIter fvIter;
		for(fIter = mesh_.faces_begin(); fIter != f_end; ++fIter)
		{
			Mesh::VertexHandle vertices[3];
			int k =0;
			for(fvIter = mesh_.fv_begin(fIter); fvIter; ++fvIter, k++)
			{
				vertices[k] = fvIter.handle();
			}

			Mesh::Point centroid = (mesh_.point(vertices[0]) +
									mesh_.point(vertices[1]) +
									mesh_.point(vertices[2]))/3;

			mesh_.property(Centroid, fIter.handle()) = centroid;
		}
	}

	std::vector<unsigned int> DeformationViewer::GetCorrespondingTriangles(Mesh::Point centroid, Mesh::Point normal, double threshold, double* distance)
	{
		Mesh::FIter fIter, f_end = mesh_.faces_end();
		
		std::vector<unsigned int> result;
		int min_id = -1;
		double min_length = DBL_MAX;
		for(fIter = mesh_.faces_begin(); fIter != f_end; ++fIter)
		{
			Mesh::Point cen = mesh_.property(Centroid, fIter.handle());
			Mesh::Point nor = mesh_.normal(fIter.handle());

			if(((cen -centroid).length() < threshold) && ((nor | normal) > 0))
			{
				if((cen - centroid).length() < min_length)
				{
					min_id = fIter.handle().idx();
					min_length = (cen - centroid).length();
				}
			}
		}

		*distance = min_length;
		result.push_back(min_id);
		return result;
	}

	void DeformationViewer::ListCorrespondances(DeformationViewer* target, double threshold, char* filename)
	{
		FILE* file = fopen(filename, "a+");
		Mesh::FIter fIter, f_end = mesh_.faces_end();
		for(fIter = mesh_.faces_begin(); fIter != f_end; ++fIter)
		{
			Mesh::Point cen = mesh_.property(Centroid, fIter.handle());
			Mesh::Point nor = mesh_.normal(fIter.handle());
			double distance;
			std::vector<unsigned int> result = target->GetCorrespondingTriangles(cen, nor, threshold, &distance);

			int size = result.size();
			for(int i = 0; i < size; i++)
			{
				fprintf(file,"%d %d\n", fIter.handle().idx(), result[i]);
			}
		}

		fclose(file);
	}

	void DeformationViewer::ListSrcCorrespondances(DeformationViewer* target, double threshold, int* tgts, double* dist, int n)
	{
		Mesh::FIter fIter, f_end = mesh_.faces_end();
		for(fIter = mesh_.faces_begin(); fIter != f_end; ++fIter)
		{
			Mesh::Point cen = mesh_.property(Centroid, fIter.handle());
			Mesh::Point nor = mesh_.normal(fIter.handle());

			double distance;
			std::vector<unsigned int> result = target->GetCorrespondingTriangles(cen, nor, threshold, &distance);
			int srcIndex = fIter.handle().idx();
			if(result.size() > 0)
			{
				int tgtIndex = result[0];	
				if(tgts[srcIndex] < 0 || distance < dist[srcIndex])
				{
					//if(n >= 0 && srcIndex >= n)
						//printf("Error\n");
					tgts[srcIndex] = tgtIndex;
					dist[srcIndex] = distance;
				}
			}

		}
	}

	void DeformationViewer::ListTgtCorrespondances(DeformationViewer* target, double threshold, int* tgts, double* dist, int n)
	{
		Mesh::FIter fIter, f_end = mesh_.faces_end();
		for(fIter = mesh_.faces_begin(); fIter != f_end; ++fIter)
		{
			Mesh::Point cen = mesh_.property(Centroid, fIter.handle());
			Mesh::Point nor = mesh_.normal(fIter.handle());
			
			double distance;
			std::vector<unsigned int> result = target->GetCorrespondingTriangles(cen, nor, threshold, &distance);
			int tgtIndex = fIter.handle().idx();
			if(result.size() > 0)
			{
				int srcIndex = result[0];	
				if(tgts[srcIndex] < 0 || distance < dist[srcIndex])
				{
					//if(n >= 0 && srcIndex >= n)
						//printf("Error\n");
					tgts[srcIndex] = tgtIndex;
					dist[srcIndex] = distance;
				}
			}
		}
	}

	double DeformationViewer::CalculateThreshold()
	{
		double x_min, y_min, z_min;
		x_min = y_min = z_min = DBL_MAX;
		double x_max, y_max, z_max;
		x_max = y_max = z_max = DBL_MIN;

		Mesh::VertexIter vIter, v_end = mesh_.vertices_end();
		for(vIter = mesh_.vertices_begin(); vIter != v_end; ++vIter)
		{
			Mesh::Point p = mesh_.point(vIter.handle());

			if(p[0] < x_min) x_min = p[0];
			if(p[1] < y_min) y_min = p[1];
			if(p[2] < z_min) z_min = p[2];
			if(p[0] > x_max) x_max = p[0];
			if(p[1] > y_max) y_max = p[1];
			if(p[2] > z_max) z_max = p[2];
		}

		double dx = x_max - x_min;
		double dy = y_max - y_min;
		double dz = z_max - z_min;

		double threshold = sqrt(4 * (dx*dy + dy*dz + dx*dz) / mesh_.n_faces());
		return threshold;
	}

	void DeformationViewer::SaveCPsToFile(char* filename)
	{
		FILE* file = fopen(filename, "w");
		if(file)
		{
			fprintf(file, "%d\n", mesh_.n_vertices());
			Mesh::VertexIter vIter, v_end = mesh_.vertices_end();
			for(vIter = mesh_.vertices_begin(); vIter != v_end; ++vIter)
			{
				int src_id = vIter.handle().idx();
				int tgt_id = mesh_.property(ClosestTargetVertexIdx, vIter);
				fprintf(file, "%d  %d\n", src_id, tgt_id);
			}

			fclose(file);
		}
	}

	int DeformationViewer::LoadCPsFromFile(char* filename)
	{
		FILE* file = fopen(filename, "r");
		if(file)
		{
			int nvertices;
			fscanf(file, "%d", &nvertices);
			if(nvertices != mesh_.n_vertices())
				return -1;
			int k =0;
			for(int i = 0; i < nvertices; i++, k++)
			{

				int src_id;
				int tgt_id;
				fscanf(file, "%d %d", &src_id, &tgt_id);
				mesh_.property(ClosestTargetVertexIdx, mesh_.vertex_handle(src_id)) = tgt_id;
			}
			if(k != mesh_.n_vertices())
				return -1;

			fclose(file);
			return 0;
		}
		return -1;
	}

//-----------------------------------------------------------------------------


/*void 
DeformationViewer::
motion(int x, int y)
{
	switch (mode_)
	{
		case MOVE:
		{
			MeshViewer::motion(x, y);
			break;
		}

			
		case DRAG:
		{
			GLdouble  modelview[16], projection[16];
			GLint     viewport[4];
			glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
			glGetDoublev(GL_PROJECTION_MATRIX, projection);
			glGetIntegerv(GL_VIEWPORT, viewport);
			double y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'


			GLdouble ox, oy, oz, wx, wy, wz;

			gluProject( moved_constraints_[active_sphere_][0],
						moved_constraints_[active_sphere_][1],
						moved_constraints_[active_sphere_][2],
						modelview,
						projection,
						viewport,
						&wx, &wy, &wz );
			
			gluUnProject( x, y_new, wz,
						  modelview, 
						  projection, 
						  viewport, 
						  &ox, &oy, &oz );
			
			moved_constraints_[active_sphere_] = Mesh::Point(ox, oy, oz);

			// compute deformation on the fly...
			deform_mesh();

			glutPostRedisplay();
			break;
		}
	}
}*/


//-----------------------------------------------------------------------------


/*void 
DeformationViewer::deform_mesh()
{
	// TODO:
	// add some deformation codes here
	// you can add whatever extra files you want
}*/


//=============================================================================