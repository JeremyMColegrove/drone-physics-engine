//
// Created by Jeremy Colegrove on 9/17/21.
//
#include <stdio.h>
#include <iostream>
#include "vector3.h"

    Vector3::Vector3() {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector3::Vector3(float _x, float _y, float _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    void Vector3::Print() {
        std::cout << '[' << x << ", " << y << ", " << z << ']' << std::endl;
    }

    Vector3 Vector3::operator+(Vector3 vec) {
        return Vector3(x + vec.x, y + vec.y, z + vec.z);
    }
    Vector3 Vector3::operator+(float a) {
        return Vector3(x + a, y + a, z + a);
    }
    Vector3 Vector3::operator-(Vector3 vec) {
        return Vector3(x - vec.x, y - vec.y, z - vec.z);
    }
    Vector3 Vector3::operator-(float a) {
        return Vector3(x - a, y - a, z - a);
    }
    Vector3 Vector3::operator*(Vector3 vec) {
        return Vector3(x * vec.x, y * vec.y, z * vec.z);
    }
    Vector3 Vector3::operator*(float a) {
        return Vector3(x * a, y * a, z * a);
    }
    Vector3 Vector3::operator/(Vector3 vec) {
        return Vector3(vec.x / x, vec.y / y, vec.z / z);
    }
    Vector3 Vector3::operator/(float a) {
        return Vector3(x / a, y / a, z / a);
    }

    Vector3 Vector3::operator+=(Vector3 vec) {
        return Vector3(x + vec.x, y + vec.y, z + vec.z);
    }
    Vector3 Vector3::operator+=(float a) {
        return Vector3(x + a, y + a, z + a);
    }
    Vector3 Vector3::operator-=(Vector3 vec) {
        return Vector3(x - vec.x, y - vec.y, z - vec.z);
    }
    Vector3 Vector3::operator-=(float a) {
        return Vector3(x - a, y - a, z - a);
    }
    Vector3 Vector3::operator*=(Vector3 vec) {
        return Vector3(x * vec.x, y * vec.y, z * vec.z);
    }
    Vector3 Vector3::operator*=(float a) {
        return Vector3(x * a, y * a, z * a);
    }
    Vector3 Vector3::operator/=(Vector3 vec) {
        return Vector3(vec.x / x, vec.y / y, vec.z / z);
    }
    Vector3 Vector3::operator/=(float a) {
        return Vector3(x / a, y / a, z / a);
    }



