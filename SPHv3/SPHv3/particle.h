#ifndef PARTICLE_H
#define PARTICLE_H

#include <list>
#include<glm\glm.hpp>
struct Particle{
	glm::dvec3 x;
	glm::dvec3 vp;
	glm::dvec3 vm;
	glm::dvec3 v;
	glm::dvec3 a;
	glm::dvec3 force;
	std::list<int> neighbours;
	double c;
	double density;
	double mass;
	double pressure;
	double k;
	double densityRest;

};

#endif