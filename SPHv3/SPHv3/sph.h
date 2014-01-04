#ifndef SPH_H
#define SPH_H
#include<glm\glm.hpp>

void updateParticles(Particle *particles, int partCount, double h, double eta, double dt, double radius);
double isColliding(glm::dvec3 x, double radius);
void collisionData(glm::dvec3 &collisionPoint, double &depth, glm::dvec3 &normal, glm::dvec3 x, double radius);
#endif