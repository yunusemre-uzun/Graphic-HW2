#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You can define helper functions inside Scene class implementation.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	Matrix4 cameraTransformationMatrix = this->CalculateCameraTransformationMatrix(camera);
	this->applyCameraTransformationToVertices(cameraTransformationMatrix); // Transform vertices
	this->applyCameraTransformationToCameras(cameraTransformationMatrix); // Transform cameras
	this->applyTransformationsToModels();

}

Matrix4 Scene::CalculateCameraTransformationMatrix(Camera *camera)
{
	Matrix4 cameraTransformationMatrix = Matrix4();
	cameraTransformationMatrix.val[0][0] = camera->u.x;
	cameraTransformationMatrix.val[0][1] = camera->u.y;
	cameraTransformationMatrix.val[0][2] = camera->u.z;
	cameraTransformationMatrix.val[1][0] = camera->v.x;
	cameraTransformationMatrix.val[1][1] = camera->v.y;
	cameraTransformationMatrix.val[1][2] = camera->v.z;
	cameraTransformationMatrix.val[2][0] = camera->w.x;
	cameraTransformationMatrix.val[2][1] = camera->w.y;
	cameraTransformationMatrix.val[2][2] = camera->w.z;
	cameraTransformationMatrix.val[3][0] = 0;
	cameraTransformationMatrix.val[3][1] = 0;
	cameraTransformationMatrix.val[3][2] = 0;
	cameraTransformationMatrix.val[0][3] = -1*(camera->u.x*camera->pos.x+camera->u.y*camera->pos.y+camera->u.z*camera->pos.z);
	cameraTransformationMatrix.val[1][3] = -1*(camera->v.x*camera->pos.x+camera->v.y*camera->pos.y+camera->v.z*camera->pos.z);
	cameraTransformationMatrix.val[2][3] = -1*(camera->w.x*camera->pos.x+camera->w.y*camera->pos.y+camera->w.z*camera->pos.z);
	cameraTransformationMatrix.val[3][3] = 1;
	return cameraTransformationMatrix;
}

void Scene::applyCameraTransformationToVertices(Matrix4 cameraTransformationMatrix)
{
	for(int i=0;i<this->vertices.size();i++) {
		Vec4 homogenous_coordinates = Vec4(this->vertices[i]->x,this->vertices[i]->y,this->vertices[i]->z,1,this->vertices[i]->colorId);
		Vec4 camera_transformed_vertice = multiplyMatrixWithVec4(cameraTransformationMatrix,homogenous_coordinates);
		this->vertices[i]->x = camera_transformed_vertice.x;
		this->vertices[i]->y = camera_transformed_vertice.y;
		this->vertices[i]->z = camera_transformed_vertice.z;
	}
}

// Note camera transformaion appllied to all cameras under the assumption a camera never visited twice
void Scene::applyCameraTransformationToCameras(Matrix4 cameraTransformationMatrix)
{
	for(int i=0;i<this->cameras.size();i++) {
		Vec4 homogenous_coordinates = Vec4(this->cameras[i]->pos.x,this->cameras[i]->pos.y,this->cameras[i]->pos.z,1,this->cameras[i]->pos.colorId);
		Vec4 camera_transformed_vertice = multiplyMatrixWithVec4(cameraTransformationMatrix,homogenous_coordinates);
		this->cameras[i]->pos.x = camera_transformed_vertice.x;
		this->cameras[i]->pos.y = camera_transformed_vertice.y;
		this->cameras[i]->pos.z = camera_transformed_vertice.z;
	}
}

void Scene::applyTransformationsToModels(void)
{
	for(int model_iterator=0;model_iterator<this->models.size();model_iterator++) {
		Model *model = this->models[model_iterator];
		for(int transformation_iterator=0;transformation_iterator<model->numberOfTransformations;transformation_iterator++) {
			switch (model->transformationTypes[transformation_iterator])
			{
			case 't':
				Translation *translation = this->translations[model->transformationIds[transformation_iterator]];
				Matrix4 translation_matrix = this->createTranslationMatrix(translation->tx,translation->ty,translation->tz);
				this->applyTransformationToModelsVertices(model, translation_matrix);
				break;
			case 's':
				Scaling *scaling = this->scalings[model->transformationIds[transformation_iterator]];
				Matrix4 scaling_matrix = this->createScalingMatrix(scaling->sx,scaling->sy,scaling->sz);
				this->applyTransformationToModelsVertices(model,scaling_matrix);
				break;
			case 'r':
				Rotation *rotation = this->rotations[model->transformationIds[transformation_iterator]];
				Matrix4 rotation_matrix = this->createRotationMatrix(rotation->angle,rotation->ux,rotation->uy,rotation->uz);
				this->applyTransformationToModelsVertices(model,rotation_matrix);
				break;
			default:
				break;
			}
		}
	}
	return;
}

Matrix4 Scene::createScalingMatrix(double sx, double sy, double sz)
{
	Matrix4 scaling_matrix = Matrix4();
	scaling_matrix.val[0][0] = sx;
	scaling_matrix.val[1][1] = sy;
	scaling_matrix.val[2][2] = sz;
	scaling_matrix.val[3][3] = 1;
	return scaling_matrix;
}

Matrix4 Scene::createTranslationMatrix(double tx, double ty, double tz)
{
	Matrix4 translation_matrix = Matrix4();
	translation_matrix.val[0][0] = 1;
	translation_matrix.val[1][1] = 1;
	translation_matrix.val[2][2] = 1;
	translation_matrix.val[3][3] = 1;
	translation_matrix.val[0][3] = tx;
	translation_matrix.val[1][3] = ty;
	translation_matrix.val[2][3] = tz;
	return translation_matrix;
}

Matrix4 Scene::createRotationMatrix(double angle, double ux, double uy, double uz)
{
	Vec3 u = Vec3(ux,uy,ux,0);
	Vec3 v = normalizeVec3(Vec3(-1*uy,ux,0,0));
	Vec3 w = normalizeVec3(crossProductVec3(u,v));
	Matrix4 m = Matrix4();
	Matrix4 m_inverse = Matrix4();
	m.val[0][0] = m_inverse.val[0][0] =ux;
	m.val[1][0] = m_inverse.val[0][1] =uy;
	m.val[2][0] = m_inverse.val[0][2] =uz;
	m.val[0][1] = m_inverse.val[1][0] =v.x;
	m.val[1][1] = m_inverse.val[1][1] =v.y;
	m.val[2][1] = m_inverse.val[1][2] =v.z;
	m.val[0][2] = m_inverse.val[2][0] =w.x;
	m.val[1][2] = m_inverse.val[2][1] =w.y;
	m.val[2][2] = m_inverse.val[2][2] =w.z;
	Matrix4 rotate_x = Matrix4();
	rotate_x.val[0][0] = 1; 
	rotate_x.val[1][1] = cos(angle);
	rotate_x.val[1][2] = -1*sin(angle);
	rotate_x.val[2][1] = sin(angle);
	rotate_x.val[2][2] = cos(angle);
	rotate_x.val[3][3] = 1;     
	Matrix4 rotation_matrix = multiplyMatrixWithMatrix(m_inverse,multiplyMatrixWithMatrix(rotate_x,m));
	return rotation_matrix;
}

void Scene::applyTransformationToModelsVertices(Model* model, Matrix4 transformation_matrix)
{
	for(int triangle_iterator=0;triangle_iterator<model->numberOfTriangles;triangle_iterator++) {
		for(int vertice=0;vertice<3;vertice++) {
			int vertex_id = model->triangles[triangle_iterator].getFirstVertexId();
			Vec3 *vertex = this->vertices[vertex_id];
			Vec4 homogenous_coordinates = Vec4(vertex->x,vertex->y,vertex->z,1,vertex->colorId);
			Vec4 transformed_vertex = multiplyMatrixWithVec4(transformation_matrix, homogenous_coordinates);
			vertex->x = transformed_vertex.x;
			vertex->y = transformed_vertex.y;
			vertex->z = transformed_vertex.z;
			vertex->colorId = transformed_vertex.colorId;
		}
	}
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
		pElement->QueryBoolText(&cullingEnabled);

	// read projection type
	pElement = pRoot->FirstChildElement("ProjectionType");
	if (pElement != NULL)
		pElement->QueryIntText(&projectionType);

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read models
	pElement = pRoot->FirstChildElement("Models");

	XMLElement *pModel = pElement->FirstChildElement("Model");
	XMLElement *modelElement;
	while (pModel != NULL)
	{
		Model *model = new Model();

		pModel->QueryIntAttribute("id", &model->modelId);
		pModel->QueryIntAttribute("type", &model->type);

		// read model transformations
		XMLElement *pTransformations = pModel->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		pTransformations->QueryIntAttribute("count", &model->numberOfTransformations);

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			model->transformationTypes.push_back(transformationType);
			model->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		// read model triangles
		XMLElement *pTriangles = pModel->FirstChildElement("Triangles");
		XMLElement *pTriangle = pTriangles->FirstChildElement("Triangle");

		pTriangles->QueryIntAttribute("count", &model->numberOfTriangles);

		while (pTriangle != NULL)
		{
			int v1, v2, v3;

			str = pTriangle->GetText();
			sscanf(str, "%d %d %d", &v1, &v2, &v3);

			model->triangles.push_back(Triangle(v1, v2, v3));

			pTriangle = pTriangle->NextSiblingElement("Triangle");
		}

		models.push_back(model);

		pModel = pModel->NextSiblingElement("Model");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	// if image is filled before, just change color rgb values with the background color
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}