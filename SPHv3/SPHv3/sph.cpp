
#include "particle.h"
#include "helper.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include<glm\glm.hpp>

bool first = true;

//laplacian of the poly6 kernel
double poly6LaplceKernel(glm::dvec3 vec, double h)
{
	double rsq = glm::dot(vec, vec);
	double hsq = h*h;
	double result = 0;

	//check for neighbourhood
	if ((0 <= rsq) && (rsq <= hsq))
	{
		result = 945 * (hsq - rsq)*(hsq - rsq)*(rsq - 0.75*(hsq - rsq)) / (8 * M_PI*pow(h, 9));
	}
	return result;
}

//gradient of the poly6 kernel
glm::dvec3 poly6GradKernel(glm::dvec3 vec, double h)
{
	double r = glm::length(vec);
	glm::dvec3 result = glm::dvec3(0, 0, 0);

	//check for neighbourhood
	if ((0 <= r) && (r <= h))
	{
		result = vec;
		result *= 945.0 / (32.0*M_PI*pow(h, 9))*(h*h - glm::length(vec)*glm::length(vec))*(h*h - glm::length(vec)*glm::length(vec));
	}
	return result;
}


//poly6 kernel
double poly6Kernel(glm::dvec3 vec, double h)
{
	double result = 0;
	double rsq = glm::dot(vec, vec);
	double hsq = h*h;

	//check for neighbourhood
	if ((0 <= rsq) && (rsq <= hsq))
	{
		result = 315.0 / (64.0*M_PI*pow(h, 9));
		result *= pow((hsq - rsq), 3);
	}
	return result;

}

//laplacien of the viscosity kernel
double laplaceViscKernel(glm::dvec3 vec, double h)
{
	double r = glm::length(vec);
	double result = 0;

	//check for neighbourhood
	if ((0 <= r) && (r <= h))
	{
		result = 45 * (h - r) / (M_PI*pow(h, 6));
	}
	return result;
}

//spiky kernel gradient 
glm::dvec3 spikyGrad(glm::dvec3 vec, double h)
{
	glm::dvec3 result = glm::dvec3(0, 0, 0);
	double r = glm::length(vec);

	//check for neighbourhood
	if ((0 <= r) && (r <= h))
	{
		result = vec;
		result *= -45.0 / (M_PI*pow(h, 6))*(h - r)*(h - r);
	}
	return result;
}



//recalculate density for each particle
void calculateDensity(Particle *particles, int partCount, double h)
{
	for (int i = 0; i < partCount; i++)
	{
		particles[i].density = 0;
		for (int j = 0; j < partCount; j++)
		{
			particles[i].density += particles[j].mass*poly6Kernel(particles[i].x - particles[j].x, h);
		}
	}
}


//recalculate pressure for each particle
void calculatePressure(Particle *particles, int partCount, double h)
{
	for (int i = 0; i < partCount; i++)
	{
		particles[i].pressure = particles[i].k*(particles[i].density - particles[i].densityRest);
	}
}


//compute and apply pressure force on each particle
void applyPressureForce(Particle *particles, int partCount, double h)
{
	for (int i = 0; i < partCount; i++)
	{
		for (int j = 0; j < partCount; j++)
		{
			particles[i].force -= particles[j].mass*(particles[i].pressure + particles[j].pressure) / (2.0*particles[j].density)*spikyGrad(particles[i].x - particles[j].x, h);
		}

	}
}

//compute and apply viscosity force on each particle
void applyViscosityForce(Particle *particles, int partCount, double h, double gamma)
{
	for (int i = 0; i < partCount; i++)
	{
		for (int j = 0; j < partCount; j++)
		{

			particles[i].force += gamma*particles[j].mass*(particles[j].v - particles[i].v) / particles[j].density*laplaceViscKernel(particles[i].x - particles[j].x, h);
		}
	}
}

//calculate acceleration add gravity for each particle
void calculateAcceleration(Particle *particles, int partCount, double dt)
{
	for (int i = 0; i < partCount; i++)
	{
		particles[i].a = particles[i].force / particles[i].density;
		//apply gravity
		particles[i].a += dt*glm::dvec3(0, 9.8, 0)*0.1;
	}
}

//updates velocities for each particle
void updateVelocity(Particle *particles, int partCount, double dt)
{
	for (int i = 0; i < partCount; i++)
	{
		particles[i].vp = particles[i].vm + particles[i].a*dt;
	}
}

//reset forces for each particle
void resetForces(Particle *particles, int partCount)
{
	for (int i = 0; i < partCount; i++)
	{
		particles[i].force = glm::dvec3(0, 0, 0);
	}

}


//calculate color field
void calcColor(Particle *particles, int partCount, int i, double h)
{
	particles[i].c = 0;
	for (int j = 0; j < partCount; j++)
	{
		particles[i].c += particles[j].mass / particles[j].density*poly6Kernel(particles[i].x - particles[j].x, h);
	}
}


//calc colorfield normal for each particle
glm::dvec3 calcNormal(Particle *particles, int partCount, double h, int i)
{
	glm::dvec3 normal = glm::dvec3(0, 0, 0);
	for (int j = 0; j < partCount; j++)
	{
		normal += particles[j].mass / particles[j].density*poly6GradKernel(particles[i].x - particles[j].x, h);
	}
	return normal;

}

//applies surface tension for each particle
void addSurfaceTension(Particle *particles, int partCount, double h, double sigma, double threshold)
{
	for (int i = 0; i < partCount; i++)
	{
		calcColor(particles, partCount, i, h);;
		glm::dvec3 normal = calcNormal(particles, partCount, h, i);

		//compute surface tension forces only if ||normal||>threshold
		if (glm::length(normal)>threshold)
		{
			for (int j = 0; j < partCount; j++)
			{
				particles[i].force -= glm::normalize(normal)*sigma*particles[j].mass / particles[j].density*poly6LaplceKernel(particles[i].x - particles[j].x, h);
			}
		}
	}
}

//calculate next position for each particle
void updateParticles(Particle *particles, int partCount, double h, double eta, double dt, double radius, double threshold, double sigma)
{
	calculateDensity(particles, partCount, h);
	calculatePressure(particles, partCount, h);

	resetForces(particles, partCount);

	applyPressureForce(particles, partCount, h);
	applyViscosityForce(particles, partCount, h, eta);
	addSurfaceTension(particles, partCount, h, sigma, threshold);

	calculateAcceleration(particles, partCount, dt);

	if (first)
	{
		for (int i = 0; i < partCount; i++)
		{
			particles[i].vm = particles[i].v - 0.5*dt*particles[i].a;
		}
		first = false;
	}
	updateVelocity(particles, partCount, dt);

	//update position and velocity
	for (int i = 0; i < partCount; i++)
	{
		particles[i].x += particles[i].vp*dt;
		particles[i].v = (particles[i].vm + particles[i].vp) / 2.0;
		particles[i].vm = particles[i].vp;
	}
}

//check for collision with bounding sphere
double isColliding(glm::dvec3 x, double radius)
{
	return glm::dot(x, x) - radius*radius;
}


//get collision data
void collisionData(glm::dvec3 &collisionPoint, double &depth, glm::dvec3 &normal, glm::dvec3 x, double radius)
{
	collisionPoint = glm::normalize(x)*radius;
	depth = abs(sqrt(glm::dot(x, x)) - radius);

	double sgn = 1.0;
	if (isColliding(x, radius) > 0)
		sgn = 1.0;
	else if (isColliding(x, radius) < 0)
	{
		sgn = -1.0;
	}
	else sgn = 0.0;

	normal = (-1.0)*glm::normalize(x);
}




