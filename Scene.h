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

typedef struct Line {
	Vec4 *starting_point;
	Vec4 *ending_point;
};

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
	vector<Vec4> copyVertices(vector<bool> &vertex_visibility);
	void applyTransformationsToModels(vector<Vec4> &);
	Matrix4 createTranslationMatrix(double tx,double ty,double tz);
	Matrix4 createScalingMatrix(double sx, double sy, double sz);
	Matrix4 createRotationMatrix(double angle, double ux, double uy, double uz);
	void applyTransformationToModelsVertices(vector<Vec4> &,Model*,Matrix4);
	Matrix4 CalculateCameraTransformationMatrix(Camera*);
	void applyCameraTransformationToVertices(vector<Vec4> &,Matrix4);
	Matrix4 createProjectionMatrix(Camera *);
	void applyProjectionMatrix(vector<Vec4> &,Matrix4);
	bool isVisible(int den, int num, double &te, double &tl);
	void clip(vector<Vec4> &copied_vertices, vector<bool> &vertex_visibility, vector<Line> &lines);
	Matrix4 createViewportMatrix(Camera *camera);
	void applyViewportTransformation(vector<Line> &lines, Matrix4 viewport_matrix);
	void createLines(vector<Line> &lines,vector<Vec4> &copied_vertices);
};

#endif
