#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;
	int projectionType;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Model* > models;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	//Here starts our helper functions
	Matrix4 CalculateCameraTransformationMatrix(Camera*);
	void applyCameraTransformationToVertices(Matrix4);
	void applyCameraTransformationToCameras(Matrix4);
	void applyTransformationsToModels(void);
	void applyTransformationToModelsVertices(Model*,Matrix4);
	Matrix4 createTranslationMatrix(double tx,double ty,double tz);
	Matrix4 createScalingMatrix(double sx, double sy, double sz);
	Matrix4 createRotationMatrix(double angle, double ux, double uy, double uz);
};

#endif
