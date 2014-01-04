#define P1 73856093
#define P2 19349663
#define P3 83492791

#include<glm\glm.hpp>
#include<vector>
#include"particle.h"
#include "hash.h"


int bitXor(unsigned int x, unsigned int y)
{
	return x ^ y;
}

unsigned int hash(glm::ivec3 intvec)
{
	unsigned int hash = 0;

	hash = intvec.x*P1 ^ intvec.y*P2 ^ intvec.z*P3;

	return hash;
}

bool isprime(unsigned int x)
{
	if (x < 2) return false;
	for (int i = 2; i <= sqrt(x); i++)
	{
		if ((x%i) == 0) return false;
	}
	return true;
}

int prime(unsigned int n)
{
	unsigned int k = 2 * n + 1;
	while (k > 0)
	{
		if (isprime(k))
			return k;
		k += 2;
	}

}


glm::ivec3 discretize(glm::dvec3 vec, double l)
{
	glm::ivec3 disc;
	disc.x = (int)(vec.x / l);
	disc.y = (int)(vec.y / l);
	disc.z = (int)(vec.z / l);

	return disc;
}

class HashTable
{
public:
	std::vector<Particle> *hashT;

	HashTable(unsigned int size)
	{
		hashT = new std::vector<Particle>[size];
	}

	unsigned int GetHash(Particle part, double h)
	{
		unsigned int pHash = hash(discretize(part.x, h));
		return pHash;
	}


	std::vector<Particle> Search(Particle part, double h)
	{
		//find the box in which the particle part is
		//find the surrounding boxes -> total of 27 boxes!
		//add every particle with dist<h to neighbour list
		//return neighbour list
	}
};

void fillHashtable(HashTable table, Particle *particles, int partCount, double h)
{

	table = HashTable(prime(partCount));

	for (int i = 0; i < partCount; i++)
	{
		unsigned int pHash = hash(discretize(particles[i].x, h));
		table.hashT[pHash].push_back(particles[i]);
	}
}
