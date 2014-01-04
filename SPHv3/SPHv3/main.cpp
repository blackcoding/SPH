#define _USE_MATH_DEFINES


#include <time.h>
#include<iostream>
#include "helper.h"
#include "particle.h"
#include "sph.h"
#include<glut.h>
#include<glm\glm.hpp>
//scene information:
//block of water in sphere

//properties
int threshold = 7.0;
int partCount = 10 * 10 * 10;
double rho0 = 1000;
double mass = rho0 * 1 * 1 * 1 / partCount;
double h = 0.2;
double k = 3;
double eta = 3.5;
double sigma = 0.728;
double eps = 0.01;
double timestep = 0.03;

//array of particles holding particle data
Particle *particles = new Particle[partCount];

//radius of the bounding sphere
double sphereRadius = 0.9;

//rendering properties

//radius of the particle drawn
double radius = 0.03;
//drawing of the bounding sphere 0=off 1=on
int drawSphere = 0;

// initial condition if sumulation is running 
bool running = false;

// angle of rotation for the camera direction
float angleX = 0.0;
float angleZ = 0.0;
float rx = 0.0f, rz = 0.0f;

void processSpecialKeys(int key, int xx, int yy)
{

	float fraction = 1.5f;

	switch (key)
	{
	case GLUT_KEY_F1:
		if (drawSphere == 1)
			drawSphere = 0;
		else
			drawSphere = 1;
		break;

	case GLUT_KEY_LEFT:
		angleZ -= fraction;
		rz = angleZ;
		break;

	case GLUT_KEY_RIGHT:
		angleZ += fraction;
		rz = angleZ;
		break;

	case GLUT_KEY_UP:
		angleX -= fraction;
		rx = angleX;
		break;

	case GLUT_KEY_DOWN:
		angleX += fraction;
		rx = angleX;
		break;

	case GLUT_KEY_END:
		if (running)
			running = false;
		else
			running = true;
		break;
	}
}

void renderParticle(Particle p, double radius)
{
	glPushMatrix();
	glTranslated(p.x.x, p.x.y, p.x.z);
	float densityfactor = p.density / p.densityRest*0.3;	
	glColor3f(densityfactor, 0.0f, 1.0f - densityfactor);
	glutSolidSphere(radius, 10, 10);
	glPopMatrix();
	glColor3f(1.0f, 1.0f, 1.0f);
}

//initiate particles position, velocity, acceleration etc.
void initParticles()
{
	//setting all to zero or simulation parameters
	for (int i = 0; i < partCount; i++)
	{
		particles[i].x = glm::dvec3(0, 0, 0);
		particles[i].v = glm::dvec3(0, 0, 0);
		particles[i].a = glm::dvec3(00.0, 0, 0.0);
		particles[i].force = glm::dvec3(0, 0, 0);
		particles[i].mass = mass;
		particles[i].pressure = 0;
		particles[i].density = 0;
		particles[i].densityRest = rho0;
		particles[i].k = k;		
	}

	//set position of particles
	int ppD = (int)cbrt(partCount);
	double dx = 1;
	double dy = dx;
	double dz = dx;
	for (int i = 0; i < partCount; i++)
	{
		//equaliy spaced in a dx*dy*dz cube
		int z = i % ppD;
		int y = (i / ppD) % ppD;
		int x = i / (ppD * ppD);
		particles[i].x = (glm::dvec3(x * dx / ppD, y * dy / ppD, z * dz / ppD) - glm::dvec3(dx / 2, dy / 2, dz / 2));
		particles[i].vm = particles[i].v;
	}
}

void init()
{
	// General Graphics
	glClearColor(0.1, 0.1, 0.1, 1);	
	glEnable(GL_DEPTH_TEST);
	//particle  initiation
	initParticles();
}

static void resize(int width, int height)
{
	const float ar = (float)width / (float)height;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glFrustum(-ar, ar, 1.0, -1.0, 1.0, 500.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glPushMatrix();
	glTranslated(0., 0., -2);
	glRotatef(rx, 1.0, 0.0, 0.0);
	glRotatef(rz, 0.0, 1.0, 0.0);
	glTranslated(0., 0., 10);
	glTranslated(0., 0., -10);
	//rendering each particle
	for (int i = 0; i < partCount; i++)
	{
		renderParticle(particles[i], radius);
	}	
	//draw bounding sphere if drawSphere==1
	if (drawSphere == 1)
		glutWireSphere(sphereRadius, 15, 15);
	glPopMatrix();
	glutSwapBuffers();
}

void idle()
{
	glm::dvec3 colorFieldGradient;
	double colorFieldLaplacian;
	// Print OpenGL errors, if there are any (for debugging)
	if (GLenum err = glGetError())
	{
		std::cerr << "OpenGL ERROR: " << gluErrorString(err) << std::endl;
	}
	
	static double previous_time;
	static bool initialized = false;

	if (!initialized)
	{
		previous_time = double(clock()) / CLOCKS_PER_SEC;
		initialized = true;
	}

	double current_time = double(clock()) / CLOCKS_PER_SEC;
	double elapsed_time = current_time - previous_time;
	double dt = timestep;// elapsed_time;

	//reset force and density
	if (running)
		updateParticles(particles, partCount, h, eta, dt, sphereRadius, threshold, sigma);

	//check each particle for collision
	for (int i = 0; i < partCount; i++)
	{
		//if colliding
		if (isColliding(particles[i].x, sphereRadius)>0)
		{
			glm::dvec3 cp;
			double depth;
			glm::dvec3 normal;
			//get collision data
			collisionData(cp, depth, normal, particles[i].x, sphereRadius);

			//collision Response
			double slip = 0.0;
			particles[i].x = cp;
			particles[i].v = particles[i].v - (1.0 + slip*depth / (dt* glm::length(particles[i].v)))* glm::dot(particles[i].v, normal)*normal;		
		}
	}
	previous_time = current_time;

	if (running)
		glutSetWindowTitle("SPH running");
	else
		glutSetWindowTitle("SPH paused");
	glutPostRedisplay();
}


int main(int argc, char** argv)
{
	printVec(glm::dvec3(1, 2, 3));
	//init GLUT and create windows
	glutInit(&argc, argv);
	glutInitWindowSize(1024, 768);
	glutInitWindowPosition(500, 110);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

	glutCreateWindow("SPH");

	init();
	glutReshapeFunc(resize);
	glutDisplayFunc(display);
	glutSpecialFunc(processSpecialKeys);
	glutIdleFunc(idle);

	// enter GLUT event processing cycle
	glutMainLoop();
}