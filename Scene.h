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
	void clip(vector<Vec4> &copied_vertices);
	Matrix4 createViewportMatrix(Camera *camera);
	void applyViewportTransformation(vector<Line> &lines, Matrix4 viewport_matrix);
	void applyViewportTransformation(vector<Vec4> &copied_vertices, Matrix4 viewport_matrix, Camera* camera);
	void createLines(vector<Line> &lines,vector<Vec4> &copied_vertices);
	void applyProjectionDivide(vector<Line> &lines);
	void applyProjectionDivide(vector<Vec4> &copied_vertices);
	Vec3 calculateVectorV(Vec3 u);
	void applyMidPointAlgorithm(vector<Vec4> &copied_vertices,Model*);
	vector<Line*> getLinesOfTriangle(Triangle triangle, vector<Vec4> &copied_vertices);
	vector<Line*> getClippedLinesOfTriangle(Triangle *triangle, vector<Vec4> &copied_vertices);
	void midPoint(Line *line);
	double calculateSlope(Vec4 *starting_point, Vec4 *ending_point);
	void midPointStandard(Line *line, bool isReflected);
	void midPointSwapped(Line *line, bool isReflected);
	void drawStandard(int x, int y, Color color);
	void drawReflected(int x, int y, int reflection_coefficient, Color color);
	void drawVerticalLine(Line *line);
	void drawHorizontalLine(Line *line);
	bool isStandardLine(Line *line);
	void swapLinePoints(Line *line);
	Color computeDc(Line *line,bool);
	Color addColor(Color,Color);
	void doBackfaceCulling(vector<Vec4> &copied_vertices);
	Vec3 calculateNormal(vector<Line*> lines);
	void render(vector<Vec4> &);
	void rasterizeTriangles(vector<Vec4>&,Model*);
	vector<Vec3> getLineEquations(vector<Line*>);
	double calculateLineEquations(double,double,Vec3);
	Color calculateAverageColor(Color* c0, Color* c1, Color* c2, double alpha, double beta, double gamma);
};

#endif
