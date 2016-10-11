#include "DeformationViewer.hh"
#include "deformationTransfer.hh"

int numPointsSelected;
bool isFirstSelected;
std::vector<Vec3f> mesh1Points,mesh2Points;

DeformationViewer *srcWindow;
DeformationViewer *targetWindow;
DeformationViewer *deformedSrcWindow;

double weight_smoothness;
double weight_identity;
double weight_closest_point_min;
double weight_closest_point_max;
double weight_closest_point_step;
int mode;



int main(int argc, char **argv)
{
	if(argc < 5 || (((atoi(argv[1]) != 0) && (atoi(argv[1]) != 1))))
	{
		std::cout<<"Command use : deformation.exe <mode:0=Compute Correspondance:1=Deformation Transfer><source mesh> <target mesh> <Correspondance List File>"<<std::endl;
		exit(0);
	}
	mode = atoi(argv[1]);

	glutInit(&argc, argv);

	numPointsSelected = 0;
	isFirstSelected = false;

	mesh1Points.clear();
	mesh2Points.clear();

	
	DeformationViewer window1("Source Window", 512, 512, true);
	DeformationViewer window2("Target Window", 512, 512, false);

	srcWindow = &window1;
	targetWindow = &window2;


	window1.open_mesh(argv[3]);
	window2.open_mesh(argv[2]);
	
	if(mode == 0)
	{
		if(argc >3)
		{
			FILE* file = fopen(argv[4],"r");
			int ncorres = 0;;
			fscanf(file,"%d", &ncorres);
			for(int i = 0; i < ncorres; i++)
			{
				int src_id, tgt_id;
				fscanf(file, "%d, %d", &src_id, &tgt_id);
				srcWindow->insertIntoConstraintList(src_id);
				srcWindow->DrawVertex(src_id);
				targetWindow->insertIntoConstraintList(tgt_id);
				targetWindow->DrawVertex(tgt_id);
			}
			fclose(file);
		}

		weight_smoothness = 1.0;
		weight_identity = 0.01;
		weight_closest_point_min = 1;
		weight_closest_point_max = 5000;
		weight_closest_point_step = 1000;


		glutMainLoop();
	}
	else
	{
		DeformationViewer window3("Source Window 2", 512, 512, true);
		window3.open_mesh(argv[5]);
		deformedSrcWindow = &window3;
		LoadCorrespondanceList(argv[4]);


		glutMainLoop();
	}
}
