#include <cmath>
#include "vector3d.h"

vector3d::vector3d()
{
    x = 0;
    y = 0;
    z = 0;
}

vector3d::vector3d(double le_pos, double wi_pos, double he_pos)
{
    x = le_pos;
    y = wi_pos;
    z = he_pos;
}

double vector3d::abs()
{
    return sqrt(x * x + y * y + z * z);
} 

vector3d& vector3d::operator += (vector3d v)
{ 
    x += v.x ;
    y += v.y ; 
    z += v.z ;

    return *this;
}

vector3d& vector3d::operator *= (double c)
{
    x *= c ; 
    y *= c ;
    z *= c ;

    return *this;
}

vector3d vector3d::operator + (vector3d v)
{
    return vector3d(x + v.x, y + v.y, z + v.z);
}

vector3d vector3d::operator - (vector3d v)
{
    return vector3d(x - v.x, y - v.y, z - v.z);
}

vector3d vector3d::operator - (vector3d v) const
{
    return vector3d(x - v.x, y - v.y, z - v.z);
}

vector3d vector3d::operator * (double c)
{
    return vector3d(x * c, y * c, z * c);
}

vector3d vector3d::operator / (double c)
{
    return vector3d(x / c, y / c, z / c);
}