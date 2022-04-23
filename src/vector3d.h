#ifndef VECTOR3D
#define VECTOR3D

struct vector3d
{   
    vector3d();
    vector3d(double le_pos, double wi_pos, double he_pos);
    
    double x;
    double y;
    double z;

    double abs();

    vector3d operator + (vector3d v);
    vector3d operator - (vector3d v);
    vector3d operator - (vector3d v) const;
    vector3d operator * (double c);
    double operator * (const vector3d &c);
    vector3d operator / (double c);
    vector3d& operator += (vector3d v);
    vector3d& operator *= (double c);

};


#endif