#ifndef VECTOR_H
#define VECTOR_H

#include<math.h>
#include<iostream>
	class Vector3
	{
	public:
		double x, y, z;
		Vector3(double X = 0, double Y = 0, double Z = 0)
		{
			x = X;
			y = Y;
			z = Z;
		}
		~Vector3() {};

		static const Vector3 Zero;

		//returns the vector with switched signs
		Vector3 Inverse()
		{
			return Vector3(-x, -y, -z);
		}

		double Dot(Vector3 vec)
		{
			return x*vec.x + y*vec.y + z*vec.z;
		}

		Vector3 Vector3::operator+(const Vector3 &vec) const
		{
			Vector3 result;
			result.x = vec.x+x;
			result.y = vec.y+y;
			result.z = vec.z+z;

			return result;
		}

		Vector3 Vector3::operator-(const Vector3 &vec) const
		{
			Vector3 result;
			result.x = vec.x - x;
			result.y = vec.y - y;
			result.z = vec.z - z;

			return result;
		}

		Vector3 operator*(double a) const
		{
			Vector3 result;
			result.x = x*a;
			result.y = y*a;
			result.z = z*a;

			return result;
		}
		friend Vector3 operator*(double a, Vector3 const &vec) 
		{
			Vector3 result;
			result.x = vec.x*a;
			result.y = vec.y*a;
			result.z = vec.z*a;

			return result;
		}

		Vector3 operator/(double a) const
		{
			if (a == 0)
				std::cout << "ERROR!"<<std::endl;
			Vector3 result;
			result.x = x/a;
			result.y = y/a;
			result.z = z/a;

			return result;
		}
		friend Vector3 operator/(double a, Vector3 const &vec)
		{
			if (a == 0)
				std::cout << "ERROR!" << std::endl;
			Vector3 result;
			result.x = vec.x/a;
			result.y = vec.y/a;
			result.z = vec.z/a;

			return result;
		}

		double Magnitude()
		{

			return sqrt(x*x + y*y + z*z);
		}

		Vector3 normalized()
		{
			double mag=sqrt(x*x + y*y + z*z);
			if (mag >= 0)
			{
				Vector3 result;
				result.x = x / mag;
				result.y = y / mag;
				result.z = z / mag;
				return result;
			}				
			else			
				return Vector3().Zero;		


		}
	};




#endif 

