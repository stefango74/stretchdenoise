//============================================================================
// Author      : Stefan Ohrhallinger
// Version     :
// Copyright   : GPL v3
// Description : Test driver for: Reconstruct a curve from noisy 2D points
//============================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <list>
#include <png++/png.hpp>
#include <GL/glew.h>
#include <GL/glut.h>

#include "Reconstruct2D.h"

using namespace std;

#define INIT_WIDTH 704	// for screen height of 768 not to squash square window
#define INIT_HEIGHT 704
int viewWidth = INIT_WIDTH;
int viewHeight = INIT_HEIGHT;

template<typename T> inline T CUB(T t) { return (SQR(t)*(t)); };

enum class PointState: char { EPSHALF, EPSONE, EPSLARGE, EPSHALF_NM, EPSONE_NM, EPSLARGE_NM };

//vector<float> noise;

class TestReconstruct2D
{
private:
	int viewWidth, viewHeight;
	float scale, diagonal;
	float offset[2], lfsMin, lfsMax;
	bool isClosed;
	int outputSize, iterations, handledFits, handledPoints, squaredFits;
	float runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum;
	vector<Point> points, projPoints, curvePoints, extremaPoints, mPoints, outputPoints, denoisedPoints;
	vector<float> lfsCurvePoints, kCurvePoints, noise;
	vector<PointState> pointStates;
	vector<Point> normals;
	vector<PointsEnum> pClasses;
	map<pair<int, int>, EdgeEnum> edgeMap;
	vector<Circle> circles;
	vector<pair<int, int> > arcs;

	void drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap);
	void drawPoints(vector<Point> &points);
	void drawPoints(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawPointsWithStates(vector<Point> &points, vector<PointState> &pointStates);
	void drawPointsWithClasses(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawCurvePoints(vector<Point> &points, vector<float> &lfsCurvePoints, float minLfs, float maxLfs);
	void drawCircle(Circle &circle);
	void drawCircles(vector<Circle> &circles);
	void drawArc(Circle &c, pair<int, int> arc);
	void drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arc);
	void drawDisc(Circle &circle);
	void drawPie(Circle &c, pair<int, int> arc);
	void drawCover(Circle &c, pair<int, int> arc);
	void drawCovers(vector<Circle> &circles, vector<pair<int, int> > &arc);
	void drawLabels(vector<Point> &points);
	void drawNormals(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawNoiseExtent(vector<Point> &points, vector<float> &noise);
	Point transform(Point);
	float scaleOfPoints(float *translate, vector<Point> &points);

public:
	TestReconstruct2D();
	virtual ~TestReconstruct2D();
	void reshape(int width, int height);
	void draw(int items, float pointSize);
	void draw(string filename, int items, float pointSize);
	void loadPointSet(string name);
	void setNoise(vector<float> p_noise);
	void invertY();
	void loadCurvePointsFromPNG(char *name);
	void generateCurveFromBezierString(string str);
	void generateCircle();
	void generateZigZag();
	void generatePointSetCosinesOnCircleWithVaryingAmplitude(float minAmp, float maxAmp);
	void generatePointSetCosinesOnCircleWithVaryingFrequency(float minFreq, float maxFreq);
	void sampleCurveByEps(float maxEps, float error, bool isClosed, float perturb);
	void sampleCurveByReach(float maxReach, float error, bool isClosed, float perturb);
	void generatePointSetBySamplingCurveEpsAndManifold(float error);
	void generateCurveCosinesOnCircleWithVaryingFrequency(int count, float minFreq, float maxFreq, float amplitude);
	bool reconstruct(int mode, int maxIter);
	void reconstruct(float minNoise, int mode, int maxIter);
	void reconstruct(vector<float> &noise, int mode, int maxIter);
	void setMaxIter(int p_maxIter);
	void computeScale();
	float getDiagonal();
	vector<Point> getPoints();
	vector<Point> getProjPoints();
	map<pair<int, int>, EdgeEnum> getEdges();
	void getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtimeC, float &p_runtimeD, float &p_oldAvgSD, float &p_newAvgSD, float &p_oldAngleSum, float &p_newAngleSum);
	void perturbPoints(float stddev);
	void perturbPointsMax(float max);
	void addCornerOutliers();
	void addNoiseToPoints(int multiplier, float maxDistance);
	void determinePointStates(float maxEps);
	float scaleOfCurvePoints();
	void generatePointsFromCircle();
	void generatePointsFrom2Circles(float noise, int seed);
	void generateNoisySharpCorner(float aspectRatio, float noise, int seed);
	void generateNoisyLeafLine(float noise, int seed);
	void generateNoisyLeafFeature(float noise, int seed);
	void generateNoisySpiral(float noise, int seed);
	void generateVaryingNoisyCircle(float noise, int seed);
	void generateVaryingDensityNoisyRect();
	void generateNoisyCircle(float noise, int seed);
	void generateZigZagRectangle();
	void generateZigZagLine();
	void generateZigZagCircle();
	void generateNoisyHole(float noise, float holeSize, int seed);
	void generateNoisyTCrossing(float noise, int seed);
	void generateNoisyFeatures(float noise, int seed);
	void generateOpenCurve();
	void generateTCurve();
	void generateSharpCorner();
	void generateHighNoise();
	void generateHighNoiseClosed();
	void generateRegularNoisyCircle();
	void generateRegularVaryingNoisyCircle();
	void generateEqualizedNoisyCircle();
	void generateVaryingNoiseCircle();
};

void TestReconstruct2D::setNoise(vector<float> p_noise)
{
	noise = p_noise;
}

/*
 * determine scale of curve points for unit square
 */
float TestReconstruct2D::scaleOfCurvePoints()
{
	float offset[2];

	return scaleOfPoints(offset, curvePoints);
}

void TestReconstruct2D::addCornerOutliers()
{
	points.push_back(Point(-0.9, -2.9));
	points.push_back(Point(4.9, 2.9));
}

vector<Point> TestReconstruct2D::getPoints()
{
	return points;
}

vector<Point> TestReconstruct2D::getProjPoints()
{
	return projPoints;
}

map<pair<int, int>, EdgeEnum> TestReconstruct2D::getEdges()
{
	return edgeMap;
}

void TestReconstruct2D::getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtimeC, float &p_runtimeD, float &p_oldAvgSD, float &p_newAvgSD, float &p_oldAngleSum, float &p_newAngleSum)
{
	p_output = outputSize;
	p_iter = iterations;
	p_fit = handledFits;
	p_point = handledPoints;
	p_squared = squaredFits;
	p_runtimeC = runtimeC;
	p_runtimeD = runtimeD;
	p_oldAvgSD = oldAvgSD;
	p_newAvgSD = newAvgSD;
	p_oldAngleSum = oldAngleSum;
	p_newAngleSum = newAngleSum;
}

/*
 * determine scale of points for unit square
 */
float TestReconstruct2D::scaleOfPoints(float *offset, vector<Point> &points)
{
	int i;
	float min[2] = { numeric_limits<float>::max(), numeric_limits<float>::max() };
	float max[2] = { -numeric_limits<float>::max(), -numeric_limits<float>::max() };

	for (auto p:points)
	{
		for (i = 0; i < 2; i++)
		{
			if (p[i] < min[i])
				min[i] = p[i];

			if (p[i] > max[i])
				max[i] = p[i];
		}
	}

	float dim[2] = { max[0] - min[0], max[1] - min[1] };

	for (i = 0; i < 2; i++)
		offset[i] = min[i];

	i = (dim[0] > dim[1]) ? 0 : 1;

	offset[1 - i] -= (dim[i] - dim[1 - i])/2;

	// set diagonal
	diagonal = sqrt(SQR(dim[0]) + SQR(dim[1]));

	return dim[i];
}

float TestReconstruct2D::getDiagonal()
{
	return diagonal;
}

bool TestReconstruct2D::reconstruct(int mode, int maxIter)
{
	Reconstruct2D *instance = new Reconstruct2D(points, noise, mode);
	instance->setMaxIter(maxIter);
	bool isClosed = instance->reconstructNoisy();
	instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum);
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	outputPoints = instance->getOutputPoints();
	denoisedPoints = instance->getDenoisedPoints();

	if (noise.size() == 0)
		noise = instance->getEstimatedNoiseExtents();

	delete instance;

	return isClosed;
}

void TestReconstruct2D::reconstruct(float minNoise, int mode, int maxIter)
{
	Reconstruct2D *instance = new Reconstruct2D(points, minNoise, mode);
	instance->setMaxIter(maxIter);
	instance->reconstructNoisy();
	instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum);
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	outputPoints = instance->getOutputPoints();
	denoisedPoints = instance->getDenoisedPoints();

	if (noise.size() == 0)
		noise = instance->getEstimatedNoiseExtents();

	delete instance;
}

void TestReconstruct2D::reconstruct(vector<float> &p_noise, int mode, int maxIter)
{
	noise = p_noise;
	Reconstruct2D *instance = new Reconstruct2D(points, noise, mode);
	instance->setMaxIter(maxIter);
	instance->reconstructNoisy();
	instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum);
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	outputPoints = instance->getOutputPoints();
	denoisedPoints = instance->getDenoisedPoints();

	delete instance;
}

void TestReconstruct2D::computeScale()
{
	scale = scaleOfPoints(offset, points);
}

/*
 * create a random float in the unit range [0..1]
 */
double unitrand()
{
	return (double)rand()/RAND_MAX;
}

void TestReconstruct2D::perturbPoints(float range)
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
	{
		float sigma = range*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += sigma*cos(angle);
		points[i][1] += sigma*sin(angle);
	}
}

/*
 * load point set from file with name
 */
void TestReconstruct2D::loadPointSet(string name)
{
	int i;

	points.clear();
	ifstream file(name);

	if (file)
	{
		while (file)
		{
			float x, y;
			file >> x >> y;

			if (file)
				points.push_back(Point(x, y));
		}
	}
	else
	{
		cerr << "ERROR: input file " << name << " could not be read." << endl;
		exit(2);
	}

	if (points.size() < 3)
	{
		cerr << "ERROR: input file " << name << " contains less than 3 points." << endl;
		exit(3);
	}

	for (i = 0; i < (int)points.size(); i++)
	{
		pClasses.push_back(CONFORMING);
		pointStates.push_back(PointState::EPSLARGE_NM);
	}

	computeScale();
}

/*
 * invert y-coordinate
 */
void TestReconstruct2D::invertY()
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
		points[i][1] = -points[i][1];
}

/*
 * converts a string two comma-separated floats into a Point
 */
Point str2point(string str)
{

	int ofs = str.find(",");
//	float t = stof(str);

//	string test = str.substr(0, ofs);
//	string test2 = str.substr(ofs + 1, str.length() - ofs - 1);
//	float t = stof(test);

	return Point(stof(str.substr(0, ofs)), stof(str.substr(ofs + 1, str.length() - ofs - 1)));
}

/*
 * parses SVG curveto string and extracts bezier control points
 */
void parseSVGCurveToString(string str, vector<Point> &bezierVec)
{
	int i;
	Point p[2];
	assert((str[0] == 'M') || (str[0] == 'm'));
	string::size_type ofs = str.find(' ', 2);
	string pStr = str.substr(2, ofs - 1);
	p[0] = str2point(pStr);
	char cc = str[ofs + 1];
	bool isRelative = (cc == 'c');
	assert(isRelative || (cc == 'C'));
	bool isEnd = false;
	ofs += 2;

	do
	{
		bezierVec.push_back(p[0]);

		for (i = 0; i < 3; i++)
		{
			int prevOfs = ofs;
			ofs = str.find(' ', ofs + 1);

			if (ofs == string::npos)
			{
				isEnd = true;
				ofs = str.length();
			}

			pStr.clear();
			pStr = str.substr(prevOfs, ofs - prevOfs);
			p[1] = str2point(pStr);

			if (isRelative)
				p[1] = p[1] + p[0];

			bezierVec.push_back(p[1]);
		}

		if (!isEnd)
			p[0] = p[1];

	} while (!isEnd);
}

/*
 * generate curve point set from SVG curveto string
 */
void TestReconstruct2D::generateCurveFromBezierString(string str)
{
	int i, j;
	const int SAMPLE_COUNT = 300;
	string bezierStr(str);
	vector<Point> b;
	parseSVGCurveToString(bezierStr, b);

	// iterate all cubic bezier curves
	for (j = 0; j < (int)b.size()/4; j++)
	{
		// sample the cubic bezier curve by evaluating with parameter t [0..1]
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			float t = (float)i/SAMPLE_COUNT;
			Point p = b[j*4]*CUB(1 - t) + b[j*4 + 1]*3*t*SQR(1 - t) + b[j*4 + 2]*3*SQR(t)*(1 - t) + b[j*4 + 3]*CUB(t);
			curvePoints.push_back(p);
		}
	}
}

/*
 * generate point set from sharp corner, perturb points by noise
 */
void TestReconstruct2D::generateNoisySharpCorner(float aspectRatio, float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float x = (float)i/SAMPLE_COUNT;
		float y = (float)i/SAMPLE_COUNT*aspectRatio;
		points.push_back(Point(x - 1.0, y));
		points.push_back(Point(x - 1.0, -y));
		points.push_back(Point(x, aspectRatio - y));
		points.push_back(Point(x, -aspectRatio + y));
	}

	perturbPoints(noise);
}

/*
 * generate point set from line, perturb points by noise
 */
void TestReconstruct2D::generateNoisyLeafLine(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));

	perturbPoints(noise);
}

/*
 * generate point set from line with feature, perturb points by noise
 */
void TestReconstruct2D::generateNoisyLeafFeature(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25, FEATURE_COUNT = 5;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));

	for (i = 0; i < FEATURE_COUNT; i++)
		points.push_back(Point(0.0, (float)i/SAMPLE_COUNT));


	perturbPoints(noise);
}

/*
 * generate point set from spiral, perturb points by noise
 * distance between boundary = 1
 * density by angle, therefore features denser sampled in interior
 */
void TestReconstruct2D::generateNoisySpiral(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;
	const int rotations = 2;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float i2 = 1.0 - SQR(1.0 - (float)i/(SAMPLE_COUNT*1.2));	// same range as i/SAMPLE_COUNT but growing superlinearly
		float angle = rotations*2*PI*i2;
		float x = (0.0 + (float)i2*rotations)*cos(angle);
		float y = (0.0 + (float)i2*rotations)*sin(angle);
		points.push_back(Point(x, y));
	}

	perturbPoints(noise);
}

/*
 * generate point set from circle, perturb points by varying noise
 */
void TestReconstruct2D::generateVaryingNoisyCircle(float noiseExtent, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;
	const int CLUSTER_COUNT = 2;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);
		points.push_back(Point(x, y));
	}

	// perturb with varying noise up to stddev 'noise'
	noise.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
	{
		float level = abs(1.0 - fmod(CLUSTER_COUNT*2*(float)i/points.size(), 2.0));
		float sigma = noiseExtent*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += level*sigma*cos(angle);
		points[i][1] += level*sigma*sin(angle);
		noise[i] = level*noiseExtent*scale;

		if (noise[i] < 0.001)
			noise[i] = 0.001;
	}
}

/*
 * generate point set from rectangle with varying density, perturb points by fixed noise
 */
void TestReconstruct2D::generateVaryingDensityNoisyRect()
{
	int i;
	const int NOISY_COUNT = 10;
	const float noise1 = 0.2, noise2 = 0.05;
	points.push_back(Point(0.0, 0.0));
	points.push_back(Point(0.0, 1.0));
	points.push_back(Point(0.0, 2.0));
	points.push_back(Point(1.0, 2.0));
	points.push_back(Point(2.0, 2.0));
	points.push_back(Point(2.0, 1.0));
	points.push_back(Point(2.0, 0.0));
	points.push_back(Point(1.0, 0.0));

	for (i = 0; i < NOISY_COUNT - 1; i++)
		points.push_back(Point((i + 1)*0.1, -noise1*scale*unitrand()));

	points.push_back(Point(1.5, 0.2));

	for (i = 0; i < (int)points.size(); i++)
	{
		points[i][0] += (noise2*unitrand() - noise2/2);
		points[i][1] += (noise2*unitrand() - noise2/2);
	}

	noise.resize(points.size());

	for (i = 0; i < (int)noise.size(); i++)
	{
		if (i < 8)
			noise[i] = noise2;
		else
			noise[i] = noise1;
	}
}

/*
 * generate point set from circle, perturb points by fixed noise
 */
void TestReconstruct2D::generateNoisyCircle(float noiseExtent, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;
//	const int CLUSTER_COUNT = 2;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);
		points.push_back(Point(x, y));
	}

	// perturb with varying noise up to stddev 'noise'
	noise.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
	{
		float sigma = noiseExtent*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += sigma*cos(angle);
		points[i][1] += sigma*sin(angle);
		noise[i] = noiseExtent*scale;
	}
}

void TestReconstruct2D::generateZigZagRectangle()
{
	int i;
	const int SIDE_COUNT = 4;
	const int SAMPLE_COUNT = 4*SIDE_COUNT;	// 4 sides
	const float STEP = 1.0;

	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SIDE_COUNT; i++)
	{
		int odd = i & 1;
		points[i] = Point(STEP/2 + i*STEP, STEP/2 - odd);
		points[SIDE_COUNT + i] = Point(STEP/2 + SIDE_COUNT*STEP + odd, STEP/2 + i*STEP);
		points[2*SIDE_COUNT + i] = Point(STEP/2 + (SIDE_COUNT - i)*STEP, STEP/2 + SIDE_COUNT*STEP + odd);
		points[3*SIDE_COUNT + i] = Point(STEP/2 - odd, STEP/2 + (SIDE_COUNT - i)*STEP);
	}

	for (i = 0; i < SAMPLE_COUNT; i++)
		noise[i] = STEP/2;
}

#ifdef OLD
void TestReconstruct2D::generateZigZagCircle()
{
	int i;
	const int SAMPLE_COUNT = 100;
	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		points[i] = Point(sin(i*2.0*PI/SAMPLE_COUNT), cos(i*2.0*PI/SAMPLE_COUNT));
		noise[i] = 0.1;
	}
}
#endif

void TestReconstruct2D::generateZigZagLine()
{
	int i;
	const int SAMPLE_COUNT = 10;
	const float STEP = 1.0;

	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points[i] = Point(STEP/2 + i*STEP, STEP/2 - (i & 1));

	for (i = 0; i < SAMPLE_COUNT; i++)
		noise[i] = STEP/2;
}

/*
 * generate point set from circle with non-random zig zag perturbation
 */
void TestReconstruct2D::generateZigZagCircle()
{
	int i;
	const int SAMPLE_COUNT = 40;
	const float NOISE = 0.1;
	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);

		if (i & 1)
			points[i] = Point(x, y)*(1.0 + NOISE);
		else
			points[i] = Point(x, y)*(1.0 - NOISE);

		noise[i] = NOISE;
	}
}

/*
 * generate point set from spiral, perturb points by noise
 */
void TestReconstruct2D::generateNoisyTCrossing(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		if ((i & 1) == 0)
			points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));
		else
			points.push_back(Point(0.5, (float)i/SAMPLE_COUNT));
	}

	perturbPoints(noise);
}

/*
 * generate sine curve with frequency FREQ and increasing amplitude, perturb points by noise
 */
void TestReconstruct2D::generateNoisyFeatures(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 50;
	const int FREQ = 5;
	const float AMP = 0.1;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float x = ((float)i)/SAMPLE_COUNT*2 - 1.0;
		float y = AMP*i*sin((float)i/SAMPLE_COUNT*2*PI*FREQ)/SAMPLE_COUNT;
		points.push_back(Point(x, y));
	}

	for (i = 0; i < (int)points.size(); i++)
		points[i][1] += noise*(2*(double)rand()/RAND_MAX - 1.0);
}

/*
 * generate point set from 2 circles
 */
void TestReconstruct2D::generatePointsFrom2Circles(float noise, int seed)
{
	int i, j;
	const int SAMPLE_COUNT = 10;

	srand(seed);

	for (j = 0; j < 2; j++)
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			double t = (double)i/SAMPLE_COUNT*2*PI;
			Point p(4*j + sin(t), cos(t));
			points.push_back(p);
		}

	// add outliers
	points.push_back(Point(2.0, -1.1));
	points.push_back(Point(2.0, 1.1));

	computeScale();
	perturbPoints(noise);
}

void TestReconstruct2D::generateRegularNoisyCircle()
{
	int i;
	const int SAMPLE_COUNT = 16;
	const float NOISE_EXTENT = 0.1;

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		double t = (double)i/SAMPLE_COUNT*2*PI;
		double r = 1.0 + 2*NOISE_EXTENT*(i % 2 - 1);
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back(NOISE_EXTENT);
	}

	computeScale();
}

void TestReconstruct2D::generateRegularVaryingNoisyCircle()
{
	int i;
	const int SAMPLE_COUNT = 32;
	const double extent = 0.5;

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		double t = (double)i/SAMPLE_COUNT*2*PI;
		double random = (double)rand()/RAND_MAX;
		double r = 1.0 + extent*(0.5 - random);
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
//		noise.push_back(0.5*random*extent);	// error
		noise.push_back(0.5*extent);
	}

	computeScale();
}

void TestReconstruct2D::generateEqualizedNoisyCircle()
{
	int i;
	const int SAMPLE_COUNT = 16, NOISY_COUNT = 8;
	const double extent = 0.1;

	for (i = 0; i < SAMPLE_COUNT - 1; i++)
	{
		double t = (double)i/SAMPLE_COUNT*2*PI;
		double r = (i == 0) ? (1.0 - extent) : 1.0;
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back((i == 0) ? extent : 0.001);
	}

	for (i = 1; i < NOISY_COUNT; i++)
	{
		double t = (double)i/NOISY_COUNT/SAMPLE_COUNT*2*PI;
		double random = (double)rand()/RAND_MAX;
		double r = 1.0 + extent*random;
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back(extent);
	}

	computeScale();
}

void TestReconstruct2D::generateVaryingNoiseCircle()
{
	int i;
	const int SAMPLE_COUNT = 8, NOISY_COUNT = 8;
	const double extent = 0.2;

	srand(1);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		double t = (double)i/SAMPLE_COUNT*2*PI;
		double r = 1.0;
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back(0.001);
	}

	for (i = 1; i < NOISY_COUNT; i++)
	{
		double t = (double)i/NOISY_COUNT/SAMPLE_COUNT*2*PI;
		double random = (double)rand()/RAND_MAX;
		double r = 1.0 + extent*(0.5 - random);
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back(abs(0.5 - random)*extent + 0.05);
//		noise.push_back(abs(0.5 - random)*extent);
	}

	computeScale();
}

/*
 * generate closed high noise point set
 */
void TestReconstruct2D::generateHighNoiseClosed()
{
	float coords[34][2] = {
			{ 0.0, 0.0 },
			{ 1.0, 0.1 },
			{ 2.0, -0.2 },
			{ 8.0, 0.2 },
			{ 9.0, -0.3 },
			{ 10.0, -0.1 },
			{ 0.0, 2.0 },
			{ 0.01, 4.0 },
			{ 0.0, 7.0 },
			{ 4.0, 7.01 },
			{ 7.0, 7.01 },
			{ 10.0, 7.0 },
			{ 10.01, 4.0 },
			{ 10.0, 2.0 },
		};

	srand(0);

	for (int i = 0; i < 20; i++)
	{
		coords[14 + i][0] = 5.0 + 2.0*(2.0*(float)rand()/RAND_MAX - 1.0);
		coords[14 + i][1] = 0.0 + 2.0*(2.0*(float)rand()/RAND_MAX - 1.0);
	}

	for (auto coord:coords)
		points.push_back(Point(coord[0], coord[1]));

	pointStates.resize(points.size());
}

/*
 * return radius for circle through point p with normalized normal n and point q
 * q can be mirrored on the line through n, therefore the radius is the circumradius of the triangle pqq'
 */
float radiusForCircleThrough2PointsandNormal(Point p, Point n, Point q)
{
	float a = p.distance(q);
	Point n2(-n[1], n[0]);
	Point pq = p - q;
	float dist = abs(pq*n2);
	float b = 2*dist;

	if (b == 0)
		return 0.5*sqrt(pq.squared_length());	// distance pq = diameter

	float e = (2*a + b)*SQR(b)*(2*a - b);	// r=abc/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) -> isosceles a=b

	if (e <= 0.0)
		return numeric_limits<float>::max();	// triangle points are collinear, infinite radius

	float d = sqrt(e);
	return abs(SQR(a)*b/d);	// circumradius of triangle adapted to isosceles version
}

/*
 * return radius for circle through point p with normalized normal n and point q
 */
float radiusForCircleThrough2PointsandNormal2(Point p, Point n, Point q)
{
	// circle center c=p+t*n, |pc|=|qc|
	double px = p[0];
	double py = p[1];
	double qx = q[0];
	double qy = q[1];
	double nx = n[0];
	double ny = n[1];
	double vx = qx - px;
	double vy = qy - py;

	return abs((SQR(vx) + SQR(vy))/(2*(vx*nx + vy*ny)));
}

/*
 * compute LFS values for curve points
 */
void computeLFSForCurve(vector<Point> &curvePoints, vector<float> &lfsCurvePoints, vector<Point> &mPoints, bool isClosed)
{
	int i, j;

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		float minR = numeric_limits<float>::max();

		for (j = 0; j < (int)curvePoints.size(); j++)
			if (i != j)
			{
				// at point p, determine radius r for maximum empty circle with center through normal n (one neighbor q on its boundary)
				Point curr2P = curvePoints[j];
				float radius = radiusForCircleThrough2PointsandNormal(currP, normal, curr2P);
				radius = radiusForCircleThrough2PointsandNormal2(currP, normal, curr2P);

				if (radius < minR)
				{
					minR = radius;
					float direction = (normal*(curr2P - currP) < 0) ? -1.0 : 1.0;
					mPoints[i] = currP + normal*(radius*direction);
				}
			}
	}

	// compute lfs from nearest medial axis points
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points;

	ann_points = annAllocPts(mPoints.size(), 2);

	for(i = 0; i < (int)mPoints.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = mPoints[i][0];
		p[1] = mPoints[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, mPoints.size(), 2);
	ANNpointArray search_point = annAllocPts(1, 2);
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray distances = new ANNdist[1];

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		// get nearest neighbor in medial axis
		search_point[0][0] = curvePoints[i][0];
		search_point[0][1] = curvePoints[i][1];
		kdTree->annkSearch(search_point[0], 1, nnIdx, distances);
		lfsCurvePoints[i] = sqrt(distances[0]);
	}

	delete nnIdx;
	delete distances;
	annDeallocPts(ann_points);
}

/*
 * return distance of p0 from line p1-p2
 */
float distancePFromLine(Point p0, Point p1, Point p2)
{
	Point normal(p2 - p1);
	swap(normal[0], normal[1]);
	normal[0] = -normal[0];
	normal.normalize();
	return abs(normal*(p0 - p1));
}

/*
 * generate point set by sampling with condition < epsMax, max error to original curve and optionally manifold condition
 */
void TestReconstruct2D::sampleCurveByEps(float maxEps, float error, bool isClosed, float perturb)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
//	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	i = 0;

	do
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		int prevI = i;
		prevP = curvePoints[prevI];
		Point newP = prevP;

		if (perturb != 0.0)
		{
			// perturb p by up to perturb*lfs(p) along its normal in any direction
			float random = (float)rand()/RAND_MAX*2 - 1.0;
			newP = newP + normal*(random*perturb*lfsCurvePoints[i]);
		}

		points.push_back(newP);
		noise.push_back(perturb*lfsCurvePoints[i]);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];

			j = prevI + 1;
			isConforming = true;

			// test eps and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for epsilon condition (a sample within dist/lfs < maxEps)
					float lfs = lfsCurvePoints[j];

					float dist = currP.distance(prevP);

					if (dist/lfs < maxEps)
						isConforming = true;
					else
					{
						dist = currP.distance(p);
						isConforming = (dist/lfs < maxEps);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

//	determinePointStates(maxEps);
	computeScale();

//	for (i = 0; i < (int)noise.size(); i++)
//		noise[i]*=scale;
}

/*
 * generate point set by sampling with condition < maxReach, max error to original curve
 */
void TestReconstruct2D::sampleCurveByReach(float maxReach, float error, bool isClosed, float perturb)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
//	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	i = 0;

	do
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		int prevI = i;
		prevP = curvePoints[prevI];
		Point newP = prevP;

		if (perturb != 0.0)
		{
			// perturb p by up to perturb*lfs(p) along its normal in any direction
			float random = (float)rand()/RAND_MAX*2 - 1.0;
			newP = newP + normal*(random*perturb*lfsCurvePoints[i]);
		}

		points.push_back(newP);
		noise.push_back(perturb*lfsCurvePoints[i]);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];
			float reach = numeric_limits<float>::max();

			for (j = prevI; j <= i; j++)
			{
				float lfs = lfsCurvePoints[j];

				if (lfs < reach)
					reach = lfs;
			}

			j = prevI + 1;
			isConforming = true;

			// test reach and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for reach condition (a sample within dist/reach < maxReach)
					// use prev/next points so that discrete interval of curve points does not impact
					float dist = curvePoints[j + 1].distance(prevP);

					if (dist/reach < maxReach)
						isConforming = true;
					else
					{
						dist = curvePoints[j - 1].distance(p);
						isConforming = (dist/reach < maxReach);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

	if (!isClosed)
		points.push_back(curvePoints[curvePoints.size() - 1]);

//	determinePointStates(maxReach);
	computeScale();
}

TestReconstruct2D::TestReconstruct2D()
{
	viewWidth = 1;
	viewHeight = 1;
	scale = 1.0;
	offset[0] = 0.0;
	offset[1] = 0.0;
	lfsMin = 0.0;
	lfsMax = 0.0;
	isClosed = true;
}

TestReconstruct2D::~TestReconstruct2D()
{
}

void TestReconstruct2D::reshape(int width, int height)
{
	this->viewWidth = width;
	this->viewHeight = height;
}

/*
 * transform point into unit square
 */
Point TestReconstruct2D::transform(Point p)
{
	int i;
	Point tp;

	for (i = 0; i < 2; i++)
		tp[i] = 0.05 + (p[i] - offset[i])/scale*0.9;	// add 5% border

	return tp;
}

void TestReconstruct2D::drawPoints(vector<Point> &points)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawPoints(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
		if (pClasses[i] == FITTED)
		{
			Point tp = transform(points[i]);
			glVertex2f(tp[0], tp[1]);
		}

	glEnd();
}

void TestReconstruct2D::drawPointsWithClasses(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		if (pClasses[i] == FITTED)
			glColor3f(0.0, 0.0, 0.0);	// black
		else
			glColor3f(0.5, 0.5, 0.5);	// gray

		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap)
{
	int i;

	glBegin(GL_LINES);

	for (auto edgeItem:edgeMap)
	{
		for (i = 0; i < 2; i++)
		{
			int index = ((i == 0) ? edgeItem.first.first : edgeItem.first.second);
			Point tp = transform(points[index]);
			glVertex2f(tp[0], tp[1]);
		}
	}

	glEnd();
}

void TestReconstruct2D::drawArc(Circle &c, pair<int, int> arc)
{
	int i;

	glBegin(GL_LINES);
	Point oldTP;

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	for (i = arc.first; i < arc.second; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point tp = transform(p);

		if (i != arc.first)
		{
			glVertex2f(oldTP[0], oldTP[1]);
			glVertex2f(tp[0], tp[1]);
		}

		oldTP = tp;
	}

	glEnd();
}

void TestReconstruct2D::drawPie(Circle &c, pair<int, int> arc)
{
	int i;

	glBegin(GL_TRIANGLES);
	Point oldTP, tc = transform(Point(c.a, c.b));

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	float degInRad = (arc.second - 1)*3.14159/180.0;
	Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
	oldTP = transform(Point(p));

	for (i = arc.first; i < arc.second; i++)
	{
		degInRad = i*3.14159/180.0;
		p = Point(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point tp = transform(p);

		glVertex2f(tc[0], tc[1]);
		glVertex2f(oldTP[0], oldTP[1]);
		glVertex2f(tp[0], tp[1]);

		oldTP = tp;
	}

	glEnd();
}

void TestReconstruct2D::drawCircle(Circle &c)
{
	pair<int, int> arc(0, 360);
	drawArc(c, arc);
}

void TestReconstruct2D::drawDisc(Circle &c)
{
	pair<int, int> arc(0, 360);
	drawPie(c, arc);
}

void TestReconstruct2D::drawCover(Circle &c, pair<int, int> arc)
{
	int i, j;

	glBegin(GL_QUADS);
	Point oldTP[2];

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	for (i = arc.first; i < arc.second; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point n(cos(degInRad)*c.variance, sin(degInRad)*c.variance);
		Point pp[2] = { p + n, p - n };
		Point tp[2];

		for (j = 0; j < 2; j++)
			tp[j] = transform(pp[j]);

		if (i != arc.first)
		{
			glVertex2f(oldTP[0][0], oldTP[0][1]);
			glVertex2f(tp[0][0], tp[0][1]);
			glVertex2f(tp[1][0], tp[1][1]);
			glVertex2f(oldTP[1][0], oldTP[1][1]);
		}

		for (j = 0; j < 2; j++)
			oldTP[j] = tp[j];
	}

	glEnd();
}

void TestReconstruct2D::drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arcs)
{
	int i;

	if (arcs.size() > 0)
		for (i = 0; i < (int)circles.size(); i++)
			drawArc(circles[i], arcs[i]);
}

void TestReconstruct2D::drawCovers(vector<Circle> &circles, vector<pair<int, int> > &arcs)
{
	int i;

	if (arcs.size() > 0)
		for (i = 0; i < (int)circles.size(); i++)
			drawCover(circles[i], arcs[i]);
}

void drawNormal(Point tp, Point n)
{
	glVertex2f(tp[0], tp[1]);
	Point endP = tp + n*0.05;
	glVertex2f(endP[0], endP[1]);

	// draw arrow hooks
	Point hookVec;
	float cs = cos(0.75*PI);
	float sn = sin(0.75*PI);
	hookVec[0] = 0.01*(cs*n[0] - sn*n[1]);
	hookVec[1] = 0.01*(sn*n[0] + cs*n[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] + hookVec[0], endP[1] + hookVec[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] - hookVec[1], endP[1] + hookVec[0]);
}

void TestReconstruct2D::drawNormals(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_LINES);

	for (i = 0; i < (int)points.size(); i++)
		if (pClasses[i] == FITTED)
		{
			// draw arrow line
			Point tp = transform(points[i]);
			drawNormal(tp, normals[i]);
		}

	glEnd();
}

void TestReconstruct2D::drawLabels(vector<Point> &points)
{
	int i, j;

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glRasterPos2f(tp[0] + 0.01, tp[1] - 0.00);
		char str[5];
		sprintf(str, "%d", i);

		for (j = 0; j < (int)strlen(str); j++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, str[j]);
	}
}

void TestReconstruct2D::drawNoiseExtent(vector<Point> &points, vector<float> &noise)
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
	{
		Circle circle(points[i][0], points[i][1], noise[i]);
//		drawCircle(circle);
		drawDisc(circle);
	}
}

const int DRAW_POINT = 1;
const int DRAW_PROJPOINT = 2;
const int DRAW_CURVEPOINT = 4;
const int DRAW_MEDIALPOINT = 8;
const int DRAW_NORMAL = 16;
const int DRAW_LABEL = 32;
const int DRAW_EDGE = 64;
const int DRAW_ARC = 128;
const int DRAW_COVER = 256;
const int DRAW_POINTCLASS = 512;
const int DRAW_DENOISEDPOINT = 1024;
const int DRAW_NOISEEXTENT = 2048;

void TestReconstruct2D::draw(int items, float pointSize)
{
	// clear screen to white
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	if (items & DRAW_COVER)
	{
		// draw covers
		glColor3f(0.9, 0.9, 1.0);	// light blue
		glLineWidth(1.0);
		drawCovers(circles, arcs);
	}

	if (items & DRAW_ARC)
	{
		// draw arcs
		glColor3f(0.0, 0.0, 1.0);	// blue
		glLineWidth(1.0);
		drawArcs(circles, arcs);
	}

	if (items & DRAW_NOISEEXTENT)
	{
		glColor3f(0.8, 0.8, 0.8);	// grey
		drawNoiseExtent(points, noise);
	}

	if (items & DRAW_POINT)
	{
		// draw points
		glColor3f(0.0, 0.0, 0.0);	// black
		glPointSize(pointSize);
		glEnable(GL_POINT_SMOOTH);
		drawPoints(points);
	}

	if (items & DRAW_POINTCLASS)
	{
		// draw points
		glPointSize(pointSize);
		glEnable(GL_POINT_SMOOTH);
		drawPointsWithClasses(points, pClasses);
	}

	if (items & DRAW_PROJPOINT)
	{
		// draw projected points
		glPointSize(pointSize);
		glColor3f(0.0, 0.0, 1.0);	// blue
		drawPoints(projPoints, pClasses);
	}

	if (items & DRAW_DENOISEDPOINT)
	{
		// draw projected points
		glPointSize(pointSize);
		glColor3f(0.0, 0.0, 1.0);	// blue
		drawPoints(denoisedPoints, pClasses);
	}

	if (items & DRAW_NORMAL)
	{
		// draw normals
		glColor3f(0.0, 0.0, 0.0);	// black
		glLineWidth(1.0);
		drawNormals(projPoints, pClasses);
	}

	if (items & DRAW_EDGE)
	{
		// draw lines
		glColor3f(1.0, 0.0, 0.0);	// red
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
		drawEdges(projPoints, edgeMap);
	}

	if (items & DRAW_LABEL)
	{
		// draw labels
		glColor3f(0.0, 0.0, 0.0);	// black
		drawLabels(points);
	}
}

/*
 * write screen to PNG file
 */
void writePNG(int width, int height, string fileName)
{
    int x, y, npixels = width*height;
	GLfloat* pixels = new GLfloat[npixels*3];
	glReadPixels(0.0, 0.0, width, height, GL_RGB, GL_FLOAT, pixels);
	png::image<png::rgb_pixel> image(width, height);
	double R, G, B;

	for (y = height - 1; y >= 0; y--)
		for (x = 0; x < width; x++)
		{
			R = *pixels++;
			G = *pixels++;
			B = *pixels++;
			image[y][x] = png::rgb_pixel(255*R, 255*G, 255*B);     // set pixel to position (x, y)
		}

	image.write(fileName);
}

/*
 * draw and write to PNG file
 */
void TestReconstruct2D::draw(string filename, int items, float pointSize)
{
	draw(items, pointSize);
	writePNG(viewWidth, viewHeight, filename);
}

static TestReconstruct2D *instance = NULL;

///////////////////////////// GL FUNCTIONS /////////////////////////////////

// GLUT idle function
void idle()
{
}

// GLUT display function
void display(string filename, int items, float pointSize)
{
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);
	glViewport(0, 0, viewWidth, viewHeight);
	instance->draw(filename, items, pointSize);

	// restore the stored viewport dimensions
	glViewport(vp[0], vp[1], vp[2], vp[3]);
	glutSwapBuffers();
}

// GLUT reshape function
void reshape(int newWidth, int newHeight)
{
	if (newWidth == 0)
		newWidth = 1;

	if (newHeight == 0)
		newHeight = 1;

	viewWidth = newWidth;
	viewHeight = newHeight;

	glViewport(0, 0, viewWidth, viewHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// update application
	instance->reshape(viewWidth, viewHeight);
}

void handleKeypress(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 27:
			exit(0);
	}
}

void do_circle(TestReconstruct2D *instance, float noise, int mode)
{
	instance->generateVaryingNoisyCircle(noise, 1);
	instance->computeScale();
	instance->reconstruct(mode, -1);
}

void fig_circle(float noise, string subfig, int mode)
{
	instance = new TestReconstruct2D();
	do_circle(instance, noise, mode);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << mode;
	display("fig_circle_" + subfig + "_" + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 10.0);
}

void fig_circle_points(float noise, string subfig)
{
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(noise, 1);
	instance->computeScale();
	reshape(viewWidth, viewHeight);
	display("fig_circle_" + subfig + ".png", DRAW_POINT | DRAW_NOISEEXTENT, 10.0);
}

void fig_circle_a(int mode)
{
	fig_circle(0.1, "a", mode);
}

void fig_circle_b(int mode)
{
	fig_circle(0.25, "b", mode);
}

void fig_circle_c(int mode)
{
	fig_circle(0.5, "c", mode);
}

void fig_circle_d(int mode)
{
	fig_circle(0.75, "d", mode);
}

void fig_circle_e(int mode)
{
	fig_circle(1.0, "e", mode);
}

void fig_circle_f(int mode)
{
	fig_circle(1.25, "f", mode);
}

void tab_circle_reconstruct_mode(int mode, float noise, float &inputMaxErr, float &inputMeanErr, float &inputRMSErr, float &outputMaxErr, float &outputMeanErr, float &outputRMSErr)
{
	int i, j;
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(noise, 1);

	// compute mean error of noisy points
	vector<Point> points = instance->getPoints();
	float sumDist = 0.0;
	float sqrSumDist = 0.0;
	inputMaxErr = 0.0;

	for (i = 0; i < (int)points.size(); i++)
	{
		float dist = abs(sqrt(SQR(points[i][0]) + SQR(points[i][1])) - 1.0);
		sumDist += dist;
		sqrSumDist += SQR(dist);

		if (dist > inputMaxErr)
			inputMaxErr = dist;
	}

	inputMeanErr = sumDist/points.size();
	inputRMSErr = sqrt(sqrSumDist/points.size());

	instance->computeScale();
	instance->reconstruct(mode, -1);

	// compute max and root mean square error of reconstructed polygon (sampled at edges with constant density)
	points = instance->getProjPoints();
	map<pair<int, int>, EdgeEnum> edgeMap = instance->getEdges();
	sumDist = 0.0;
	sqrSumDist = 0.0;
	outputMaxErr = 0.0;
	int count = 0;

	for (auto edgeItem:edgeMap)
	{
		Point p0 = points[edgeItem.first.first], p1 = points[edgeItem.first.second];
		Point v = p1 - p0;
		float edgeLen = sqrt(v.squared_length());

		for (j = 0; j < (int)(edgeLen*100.0); j++)
		{
			float len = 0.01/edgeLen;
			Point p = p0 + v*j*len;
			float dist = abs(sqrt(SQR(p[0]) + SQR(p[1])) - 1.0);
			sumDist += dist;
			sqrSumDist += SQR(dist);

			if (dist > outputMaxErr)
				outputMaxErr = dist;

			count++;
		}
	}

	outputMeanErr = sumDist/count;
	outputRMSErr = sqrt(sqrSumDist/count);
}

void tab_circle_reconstruct(stringstream &sout, float noise)
{
	float inputMaxErr, inputMeanErr, inputRMSErr, outputMaxErr, outputMeanErr, outputRMSErr;

	tab_circle_reconstruct_mode(MODE_BLEND, noise, inputMaxErr, inputMeanErr, inputRMSErr, outputMaxErr, outputMeanErr, outputRMSErr);
	sout << noise << "\t| " << inputMaxErr << "\t| " << inputMeanErr << "\t| " << inputRMSErr << "\t| " << outputMaxErr << "\t| " << outputMeanErr << "\t| " << outputRMSErr;
	tab_circle_reconstruct_mode(MODE_DENOISE, noise, inputMaxErr, inputMeanErr, inputRMSErr, outputMaxErr, outputMeanErr, outputRMSErr);
	sout << "\t| " << outputMaxErr << "\t| " << outputMeanErr << "\t| " << outputRMSErr << endl;
}

void tab_circle()
{
	int i;
	stringstream sout;
	sout.setf(ios::fixed, ios::floatfield);
	sout.precision(3);
	float noise[] = { 0.1, 0.25, 0.5, 0.75, 1.0 };
	sout << "Table 2:" << endl;
	sout << "Noise \t| maxI \t| meanI\t| rmsI \t| maxO \t| meanO\t| rmsO \t| maxD \t| meanD\t| rmsD" << endl;
	sout << "======================================================================================" << endl;

	for (i = 0; i < 5; i++)
		tab_circle_reconstruct(sout, noise[i]);

	// output table
	cout << sout.str();
}

void do_rect(TestReconstruct2D *instance, int mode)
{
	srand(1);
	instance->generateVaryingDensityNoisyRect();
	instance->computeScale();
	instance->reconstruct(mode, -1);
}

void fig_rect(int mode)
{
	instance = new TestReconstruct2D();
	do_rect(instance, mode);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << mode;
	display("fig_rect_" + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 10.0);
}

void fig_rects()
{
	fig_rect(MODE_NONE);
	fig_rect(MODE_DENOISE);
}

void do_varying_subset(TestReconstruct2D *instance, int mode)
{
	srand(1);
	instance->generateVaryingNoiseCircle();
	instance->computeScale();
	instance->reconstruct(mode, -1);
}

void fig_varying_subset(int mode)
{
	instance = new TestReconstruct2D();
	do_varying_subset(instance, mode);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << mode;
	display("fig_varying_subset_" + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 10.0);
}

void fig_silhouette(string pointsetname, int index, float noiseVar, bool denoise, bool invert, bool curveonly, bool minNoise)
{
	vector<float> noise;
	float pointsize = 2;
	srand(1);

	instance = new TestReconstruct2D();
	instance->loadPointSet("data/" + pointsetname + ".txt");

	if (invert)
		instance->invertY();

	instance->computeScale();

	if (noiseVar == 0.0)
		instance->reconstruct(denoise ? MODE_DENOISE : MODE_NONE, -1);
	else
	if (minNoise)
		instance->reconstruct(noiseVar, denoise ? MODE_DENOISE : MODE_NONE, -1);
	else
	{
		vector<float> noise;

		for (int i = 0; i < (int)instance->getPoints().size(); i++)
			noise.push_back(noiseVar);

		instance->reconstruct(noise, denoise ? MODE_DENOISE : MODE_NONE, -1);
	}

	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << index;
	display("fig_" + pointsetname + "_" + ss.str() + "_denoised.png", (curveonly ? 0 : (DRAW_POINT | DRAW_NOISEEXTENT)) | DRAW_EDGE , pointsize);
}

void fig_silhouette_a1()
{
	fig_silhouette("Mouse0", 1, 0.0, false, true, false, true);
}

void fig_silhouette_a2()
{
	fig_silhouette("Mouse0", 2, 0.001, true, true, false, true);
}

void fig_silhouette_b1()
{
	fig_silhouette("Monitor0", 1, 0.0, false, true, false, true);
}

void fig_silhouette_b2()
{
	fig_silhouette("Monitor0", 2, 0.001, true, true, false, true);
}

void fig_silhouette_c1()
{
	fig_silhouette("Keyboard0", 1, 0.0, false, true, false, true);
}

void fig_silhouette_c2()
{
	fig_silhouette("Keyboard0", 2, 0.001, true, true, false, true);
}

void fig_silhouette_d1()
{
	fig_silhouette("Cup0", 1, 0.0, false, true, false, true);
}

void fig_silhouette_d2()
{
	fig_silhouette("Cup0", 2, 0.001, true, true, false, true);
}

void fig_silhouette_e1()
{
	fig_silhouette("drill", 1, 0.0, false, false, false, false);
}

void fig_silhouette_e2()
{
	fig_silhouette("drill", 2, 0.003, true, false, false, false);
}

void fig_silhouette_e3()
{
	fig_silhouette("drill", 3, 0.0, false, false, true, false);
}

void fig_silhouette_e4()
{
	fig_silhouette("drill", 4, 0.003, true, false, true, false);
}

void fig_silhouette_f1()
{
	fig_silhouette("thing", 1, 0.0, false, false, false, false);
}

void fig_silhouette_f2()
{
	fig_silhouette("thing", 2, 0.003, true, false, false, false);
}

void fig_silhouette_f3()
{
	fig_silhouette("thing", 3, 0.0, false, false, true, false);
}

void fig_silhouette_f4()
{
	fig_silhouette("thing", 4, 0.003, true, false, true, false);
}

void fig_silhouette_g1()
{
	fig_silhouette("vase", 1, 0.0, false, false, false, false);
}

void fig_silhouette_g2()
{
	fig_silhouette("vase", 2, 0.003, true, false, false, false);
}

void fig_silhouette_g3()
{
	fig_silhouette("vase", 3, 0.0, false, false, true, false);
}

void fig_silhouette_g4()
{
	fig_silhouette("vase", 4, 0.003, true, false, true, false);
}

void do_comparison(TestReconstruct2D *instance, string filename, float noiseVar, bool denoise)
{
	int i;
	vector<float> noise;
	instance->loadPointSet("./data/" + filename + ".txt");
	instance->invertY();

	if (noiseVar != 0.0)
	{
		for (i = 0; i < (int)instance->getPoints().size(); i++)
			noise.push_back(noiseVar);
	}

	instance->setNoise(noise);
	instance->computeScale();
	instance->reconstruct(denoise ? MODE_DENOISE : MODE_NONE, -1);
}

void fig_comparison(string filename, string figname, float noiseVar, float pointSize, bool denoise)
{
	instance = new TestReconstruct2D();
	do_comparison(instance, filename, noiseVar, denoise);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << (denoise ? 1 : 0);
	display(figname + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | ((noiseVar != 0.0) ? DRAW_NOISEEXTENT : 0), pointSize);
}

void fig_comparison_a1()
{
	fig_comparison("apple_2percent_noise", "fig_comparison_a", 0.03, 7.0, false);
}

void fig_comparison_a2()
{
	fig_comparison("apple_2percent_noise", "fig_comparison_a", 0.03, 7.0, true);
}

void fig_comparison_b1()
{
	fig_comparison("butterfly_2percent_noise", "fig_comparison_b", 0.03, 7.0, false);
}

void fig_comparison_b2()
{
	fig_comparison("butterfly_2percent_noise", "fig_comparison_b", 0.03, 7.0, true);
}

void fig_comparison_c1()
{
	fig_comparison("crab_2percent_noise", "fig_comparison_c", 0.03, 7.0, false);
}

void fig_comparison_c2()
{
	fig_comparison("crab_2percent_noise", "fig_comparison_c", 0.03, 7.0, true);
}

void fig_comparison_d1()
{
	fig_comparison("dolphin_2percent_noise", "fig_comparison_d", 0.03, 7.0, false);
}

void fig_comparison_d2()
{
	fig_comparison("dolphin_2percent_noise", "fig_comparison_d", 0.03, 7.0, true);
}

void fig_comparison2_a1()
{
	fig_comparison("fish", "fig_comparison2_a", 0.5, 5.0, false);
}

void fig_comparison2_a2()
{
	fig_comparison("fish", "fig_comparison2_a", 0.5, 5.0, true);
}

void fig_comparison2_b1()
{
	fig_comparison("bot", "fig_comparison2_b", 0.5, 5.0, false);
}

void fig_comparison2_b2()
{
	fig_comparison("bot", "fig_comparison2_b", 0.5, 5.0, true);
}

void do_highnoise(TestReconstruct2D *instance, bool denoise)
{
	vector<float> noise;
	srand(1);
	string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";
	instance->generateCurveFromBezierString(bunny); instance->sampleCurveByReach(0.01, 0.0, true, 1.0/3);
	instance->computeScale();
	instance->reconstruct(denoise ? MODE_DENOISE : MODE_NONE, -1);
}

void fig_highnoise(float noiseVar, bool denoise)
{
	instance = new TestReconstruct2D();
	do_highnoise(instance, denoise);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << (denoise ? 1 : 0);
	display("fig_highnoise_" + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 5);
}

void fig_highnoise_1()
{
	fig_highnoise(0.0, false);
}

void fig_highnoise_2()
{
	fig_highnoise(0.0, true);
}

void do_lfs(TestReconstruct2D *instance, float rho, float noise, bool denoise)
{
	string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";
	srand(1);
	instance->generateCurveFromBezierString(bunny);
	instance->sampleCurveByReach(rho, 0.0, true, noise);	// rho=0.43 -> eps>0.3
	instance->computeScale();
	instance->reconstruct(denoise ? MODE_DENOISE : MODE_NONE, -1);
}

void do_density(TestReconstruct2D *instance, float eps, bool denoise)
{
	do_lfs(instance, eps/(1.0 - eps), 1.0/3.0, denoise);
}

void fig_density(float eps, string subfig, bool denoise)
{
	instance = new TestReconstruct2D();
	do_density(instance, eps, denoise);
	reshape(viewWidth, viewHeight);
	display("fig_density_" + subfig + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 10.0);
}

void fig_density_a1()
{
	fig_density(0.4, "a1", false);
}

void fig_density_a2()
{
	fig_density(0.4, "a2", true);
}

void fig_density_b1()
{
	fig_density(0.3, "b1", false);
}

void fig_density_b2()
{
	fig_density(0.3, "b2", true);
}

void fig_density_c1()
{
	fig_density(0.2, "c1", false);
}

void fig_density_c2()
{
	fig_density(0.2, "c2", true);
}

void fig_density_d1()
{
	fig_density(0.1, "d1", false);
}

void fig_density_d2()
{
	fig_density(0.1, "d2", true);
}

void fig_circles()
{
	fig_circle_points(0.1, "a");
	fig_circle_points(0.25, "b");
	fig_circle_points(0.5, "c");
	fig_circle_points(0.75, "d");
	fig_circle_points(1.0, "e");
	fig_circle_a(MODE_NONE);
	fig_circle_b(MODE_NONE);
	fig_circle_c(MODE_NONE);
	fig_circle_d(MODE_NONE);
	fig_circle_e(MODE_NONE);
	fig_circle_a(MODE_BLEND);
	fig_circle_b(MODE_BLEND);
	fig_circle_c(MODE_BLEND);
	fig_circle_d(MODE_BLEND);
	fig_circle_e(MODE_BLEND);
	fig_circle_a(MODE_DENOISE);
	fig_circle_b(MODE_DENOISE);
	fig_circle_c(MODE_DENOISE);
	fig_circle_d(MODE_DENOISE);
	fig_circle_e(MODE_DENOISE);
}

void fig_bunny_orig()
{
	instance = new TestReconstruct2D();
	do_lfs(instance, 0.1/(1.0 - 0.1), 0.0, MODE_NONE);	// eps=0.1
	reshape(viewWidth, viewHeight);
	display("fig_bunny_orig.png", DRAW_EDGE, 10.0);
}

void fig_densities()
{
	fig_varying_subset(MODE_NONE);
	fig_varying_subset(MODE_DENOISE);
	fig_density_a1();
	fig_density_a2();
	fig_density_b1();
	fig_density_b2();
	fig_density_c1();
	fig_density_c2();
	fig_density_d1();
	fig_density_d2();
	fig_bunny_orig();
}

void fig_comparisons()
{
	fig_comparison_a1();
	fig_comparison_a2();
	fig_comparison_b1();
	fig_comparison_b2();
	fig_comparison_c1();
	fig_comparison_c2();
	fig_comparison_d1();
	fig_comparison_d2();
}

void fig_high_noise()
{
	fig_comparison2_a2();
	fig_comparison2_b2();
	fig_highnoise_2();
}

void fig_silhouettes1()
{
	fig_silhouette_a1();
	fig_silhouette_a2();
	fig_silhouette_b1();
	fig_silhouette_b2();
	fig_silhouette_c1();
	fig_silhouette_c2();
	fig_silhouette_d1();
	fig_silhouette_d2();
}

void fig_silhouettes2()
{
	fig_silhouette_e3();
	fig_silhouette_e4();
	fig_silhouette_f3();
	fig_silhouette_f4();
	fig_silhouette_g3();
	fig_silhouette_g4();
}

void TestReconstruct2D::generateZigZag()
{
	int i;
	const int SAMPLE_COUNT_1 = 16, SAMPLE_COUNT_2 = 16;

	for (i = 1; i < SAMPLE_COUNT_1; i++)
	{
		double t = -PI/2 + (double)i/SAMPLE_COUNT_1*PI;
		double r = 1.0;
		Point p(r*sin(t), r*cos(t));
		points.push_back(p);
		noise.push_back(0.001);
	}

	for (i = 1; i < SAMPLE_COUNT_2; i++)
	{
		points.push_back(Point((2.0*(double)i/SAMPLE_COUNT_2 - 1.0), (double)(i & 1)*2.0/SAMPLE_COUNT_2));
		noise.push_back((0.5 + (double)i/SAMPLE_COUNT_2)/SAMPLE_COUNT_2);
	}

	computeScale();
}

void do_zigzag(TestReconstruct2D *instance, int mode)
{
	instance->generateZigZag();
	instance->computeScale();
	instance->reconstruct(mode, -1);
}

void fig_zigzag(int mode)
{
	instance = new TestReconstruct2D();
	do_zigzag(instance, mode);
	reshape(viewWidth, viewHeight);
	stringstream ss;
	ss << mode;
	display("fig_zigzag_" + ss.str() + ".png", DRAW_POINT | DRAW_EDGE | DRAW_NOISEEXTENT, 10.0);
}

void fig_zigzags()
{
	fig_zigzag(MODE_NONE);
	fig_zigzag(MODE_DENOISE);
}

void tab_runtime()
{
	int i, outputSize, iterations, handledFits, handledPoints, squaredFits;
	float runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum;
	string name[] = { "{\\scshape Circle} 0.1r", "{\\scshape Circle} 0.25r", "{\\scshape Circle} 0.5r", "{\\scshape Circle} 0.75r", "{\\scshape Circle} r",
			"{\\scshape Bunny} $\\epsilon=0.4$", "{\\scshape Bunny} $\\epsilon=0.3$", "{\\scshape Bunny} $\\epsilon=0.2$", "{\\scshape Bunny} $\\epsilon=0.1$", "{\\scshape Keyboard}", "{\\scshape Monitor}", "{\\scshape Cup}", "{\\scshape Mouse}",
			"{\\scshape Apple}", "{\\scshape Butterfly}", "{\\scshape Crab}", "{\\scshape Dolphin}", "{\\scshape Fish}", "{\\scshape Bottle}", "{\\scshape Bunny} hi noise", "{\\scshape VarCircle}", "{\\scshape Square}", "{\\scshape Sawtooth}"
	};

	stringstream sout;
	sout.setf(ios::fixed, ios::floatfield);
	sout.precision(3);

	sout << "Table 1:" << endl;
	sout << "      Figure/Noise & \\# & old $\\overline{SD}$ & new $\\overline{SD}$ & old$\\sum\\angle$ & new$\\sum\\angle$ & Conn & Den \\\\\n    \\hline\n";

	for (i = 0; i < 23; i++)
	{
		instance = new TestReconstruct2D();

		if (i == 0)
			do_circle(instance, 0.1, MODE_DENOISE);
		else
		if (i == 1)
			do_circle(instance, 0.25, MODE_DENOISE);
		else
		if (i == 2)
			do_circle(instance, 0.5, MODE_DENOISE);
		else
		if (i == 3)
			do_circle(instance, 0.75, MODE_DENOISE);
		else
		if (i == 4)
			do_circle(instance, 1.0, MODE_DENOISE);
		else
		if (i == 5)
			do_density(instance, 0.4, true);
		else
		if (i == 6)
			do_density(instance, 0.3, true);
		else
		if (i == 7)
			do_density(instance, 0.2, true);
		else
		if (i == 8)
			do_density(instance, 0.1, true);
		else
		if (i == 9)
			do_comparison(instance, "Keyboard0", 0.0, true);
		else
		if (i == 10)
			do_comparison(instance, "Monitor0", 0.0, true);
		else
		if (i == 11)
			do_comparison(instance, "Cup0", 0.0, true);
		else
		if (i == 12)
			do_comparison(instance, "Mouse0", 0.0, true);
		else
		if (i == 13)
			do_comparison(instance, "apple_2percent_noise", 0.03, true);
		else
		if (i == 14)
			do_comparison(instance, "butterfly_2percent_noise", 0.03, true);
		else
		if (i == 15)
			do_comparison(instance, "crab_2percent_noise", 0.03, true);
		else
		if (i == 16)
			do_comparison(instance, "dolphin_2percent_noise", 0.03, true);
		else
		if (i == 17)
			do_comparison(instance, "fish", 0.5, true);
		else
		if (i == 18)
			do_comparison(instance, "bot", 0.5, true);
		else
		if (i == 19)
			do_highnoise(instance, true);
		else
		if (i == 20)
			do_varying_subset(instance, MODE_DENOISE);
		else
		if (i == 21)
			do_rect(instance, MODE_DENOISE);
		else
		if (i == 22)
			do_zigzag(instance, MODE_DENOISE);

		instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtimeC, runtimeD, oldAvgSD, newAvgSD, oldAngleSum, newAngleSum);
		vector<Point> points = instance->getPoints();
		float diagonal = instance->getDiagonal();
		oldAvgSD *= 100.0/diagonal;
		newAvgSD *= 100.0/diagonal;
		sout << "      " << name[i] << " & " << points.size() << " & " << oldAvgSD << " & " << newAvgSD << " & ";
		sout.precision(0);
		sout << oldAngleSum << " & " << newAngleSum << " & ";
		sout.precision(3);
		sout << runtimeC << " & " << runtimeD << " \\\\" << endl;
		delete instance;
	}

	// output table
	cout << sout.str();
}

void fig_bunny_points()
{
	instance = new TestReconstruct2D();
	do_density(instance, 0.1, false);
	reshape(viewWidth, viewHeight);
	display("fig_bunny_points.png", DRAW_POINT | DRAW_NOISEEXTENT, 10.0);
}

void printUsage(char *argv[])
{
	cout << "Usage: " << argv[0] << " fig1|fig4|fig5|fig6|fig7|fig8|fig9|fig10|fig11|tab1|tab2" << endl;
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowSize(viewWidth, viewHeight);
	glutCreateWindow("Reconstruct from 2D Points");

	if (argc <= 1)
		printUsage(argv);
	else
	{
		string a = argv[1];

		if (!a.compare("fig1"))
		{
			fig_bunny_orig();
			fig_bunny_points();
			fig_density_d1();
			fig_density_d2();
 		}
		else
		if (!a.compare("fig4"))
			fig_circles();
		else
		if (!a.compare("fig5"))
			fig_densities();
		else
		if (!a.compare("fig6"))
			fig_comparisons();
		else
		if (!a.compare("fig7"))
			fig_high_noise();
		else
		if (!a.compare("fig8"))
			fig_zigzags();
		else
		if (!a.compare("fig9"))
			fig_silhouettes1();
		else
		if (!a.compare("fig10"))
			fig_silhouettes2();
		else
		if (!a.compare("fig11"))
			fig_rects();
		else
		if (!a.compare("tab1"))
			tab_runtime();
		else
		if (!a.compare("tab2"))
			tab_circle();
		else
		{
			cout << "No valid parameter chosen" << endl;
			printUsage(argv);
		}
	}

	return 0;
}
