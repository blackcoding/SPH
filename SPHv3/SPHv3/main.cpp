#define _USE_MATH_DEFINES
#define GLFW_INCLUDE_GLU

#include <time.h>
#include<iostream>
#include"vector.h"
#include "helper.h"
#include "particle.h"
#include "sph.h"
#include <math.h>
#include<glew.h>
#include<glut.h>
#include<glm\glm.hpp>
//scene information:
//block of water 1x1x1m ->1000kg mass;



double totalmass =1000;

//properties
int threshold =7.0;
int partCount =8*8*8;
double rho0 =1000;
double mass = rho0*1*1*1/partCount;// totalmass / partCount;
double h =0.2;// 1.1255 + 0.2;
double k =3;
double eta = 3.5;
double sigma = 0.728;
double eps = 0.01;
//global variables
double boxSize = 6;

double radius = 0.03;
double sphereRadius =0.9;

double timePassed = 0;

//glut stuff
int drawSphere = 0;
// angle of rotation for the camera direction
float angleX = 0.0;
float angleZ = 0.0;

// actual vector representing the camera's direction
bool running = false;
float rx = 0.0f, rz = 0.0f;



Particle *particles = new Particle[partCount];




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
	float densityfactor = p.density/p.densityRest*0.8;
	//GLfloat mat_specular[] = { densityfactor, 1.0, 1.0f - densityfactor, 1.0 };
	glColor3f(densityfactor, 0.0f, 1.0f - densityfactor);

		
	
	
		
	
	

	glutSolidSphere(radius, 10, 10);
	
	glPopMatrix();
	glColor3f(1.0f, 1.0f, 1.0f);
}

//initiate particles position, velocity, acceleration
void initParticles()
{
	//setting all to zero or simulation parameters
	for (int i = 0; i < partCount;i++)
	{
		particles[i].x = glm::dvec3(0, 0, 0);
		particles[i].v = glm::dvec3(0, 0, 0);
		particles[i].a = glm::dvec3(0, 0, 0);
		particles[i].force = glm::dvec3(0, 0, 0);
		particles[i].mass = mass;
		particles[i].pressure = 0;
		particles[i].density = 0;
		particles[i].densityRest = rho0;
		particles[i].k = k;
		//particles[i].eta = eta;
	}

	//set position of particles
	int ppD = (int)cbrt(partCount);
	double dx =1;
	double dy = dx;
	double dz = dx;
	for (int i = 0; i < partCount; i++)
	{
		//equaliy spaced in a 5x5x5 cube
		
		
		int z = i % ppD;
		int y = (i / ppD) % ppD;
		int x = i / (ppD * ppD);
		particles[i].x = (glm::dvec3(x * dx / ppD, y * dy / ppD, z * dz / ppD) - glm::dvec3(dx/2,dy/2, dz/2));
		particles[i].vm = particles[i].v;
	}
	//compute missing parameters:
	//density
	
	//updateDensity(particles, partCount,h);
	//pressure
	//updatePressure(particles, partCount);
}

void init()
{
	// General Graphics
	
	glClearColor(0.1, 0.1, 0.1, 1);


	glEnable(GL_DEPTH_TEST);





	initParticles();
}

static void resize(int width, int height)
{
	const float ar = (float)width / (float)height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);


	glFrustum(-ar, ar, 1.0, -1.0, 1.0, 500.0);
	//gluPerspective(60, double(width) / height, 1.0, 10.0);
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
	glTranslated(0., 0.,-10);
	for (int i = 0; i < partCount; i++)
	{
		renderParticle(particles[i], radius);
		if (i == 0)
		{
			glColor3f(0.0f, 1.0f, 1.0f);
			glBegin(GL_LINES); 
			glVertex3d(particles[i].x.x, particles[i].x.y, particles[i].x.z);

			glVertex3d(particles[i].x.x + particles[i].a.x, particles[i].x.y + particles[i].a.y, particles[i].x.z + particles[i].a.z);
			glEnd();
			/*glTranslated(particles[i].x.x, particles[i].x.y, particles[i].x.z);
			glutWireSphere(h, 15,15);
			glTranslated(-particles[i].x.x, -particles[i].x.y, -particles[i].x.z);*/
			//glLoadIdentity();
			glFlush();
			
		}
		
	}
		
	//glTranslated(0., 0., -2.5);
	//glutWireCube(boxSize);
	if (drawSphere==1)
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
	double dt = 0;
	
	dt = 0.03;// elapsed_time;
	
	//reset force and density
	if (running)
		updateParticles(particles, partCount, h, eta, dt, sphereRadius,threshold,sigma);
	for (int i = 0; i < partCount; i++)
	{
		
		if (isColliding(particles[i].x,sphereRadius)>0)
		{
			glm::dvec3 cp;
			double depth;
			glm::dvec3 normal;
			collisionData(cp, depth, normal, particles[i].x, sphereRadius);		
			//collision Response
			double slip =0.0;
			particles[i].x = cp;
			
			particles[i].v =  particles[i].v - (1.0 + slip*depth / (dt* glm::length(particles[i].v)))* glm::dot(particles[i].v,normal)*normal;
			//particles[i].v = particles[i].v.Inverse()*0.3;// -2.0*particles[i].v.Dot(normal)*normal;
		}
	}
	

	
	//std::cout << "render" << std::endl;

	previous_time = current_time;
	if (running)
		glutSetWindowTitle("SPH running");
	else
		glutSetWindowTitle("SPH");
	glutPostRedisplay();
}


int main(int argc, char** argv)
{
	printVec(glm::dvec3(1, 2,3));
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

	//initialize particles:

}