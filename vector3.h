//
// Created by Jeremy Colegrove on 9/17/21.
//

#ifndef LAB02_CLASS_BASICS_VECTOR3_H
#define LAB02_CLASS_BASICS_VECTOR3_H

class Vector3 {

public:
    float x, y, z;

    Vector3();
    Vector3(float _x, float _y, float _z);
    void Print();

    Vector3 operator+(Vector3 vec);
    Vector3 operator+(float a);
    Vector3 operator-(Vector3 vec);
    Vector3 operator-(float a);
    Vector3 operator*(Vector3 vec);
    Vector3 operator*(float a);
    Vector3 operator/(Vector3 vec);
    Vector3 operator/(float a);

    Vector3 operator+=(Vector3 vec);
    Vector3 operator+=(float a);
    Vector3 operator-=(Vector3 vec);
    Vector3 operator-=(float a);
    Vector3 operator*=(Vector3 vec);
    Vector3 operator*=(float a);
    Vector3 operator/=(Vector3 vec);
    Vector3 operator/=(float a);
};


#endif //LAB02_CLASS_BASICS_VECTOR3_H
