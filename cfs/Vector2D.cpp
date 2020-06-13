//
//  Vector2D.cpp
//  ToolPathTimeEstimator
//
//  Created by Celr on 15/12/10.
//  Copyright © 2015年 Celr. All rights reserved.
//
#include "stdafx.h"
#include "Vector2D.hpp"


Vector2D::Vector2D(double x, double y) {
    this->x = x;
    this->y = y;
	this->direction = false;
    this->value = sqrt(x*x+y*y);
}
Vector2D::Vector2D() {
    
}

double Vector2D::get_value() {
	this->value = sqrt(x*x + y*y);
	return this->value;
}

double Vector2D::maxAlongDirection(Vector2D direction) {
	if (!this->direction) {
		return this->x;
	}
    if (direction.x == 0) {
        return fabs(this->y);
    }
    if (direction.y == 0) {
        return fabs(this->x);
    }
    
	if (fabs(direction.x)>fabs(direction.y)) {
        return sqrt(pow(this->x,2)+pow(this->y*(direction.y/direction.x), 2));
    }else {
        return sqrt(pow(this->y,2)+pow(this->x*(direction.x/direction.y), 2));
    }
}

double Vector2D::tan() {
    return this->x/this->y;
}

double Vector2D::sin() {
    return this->x/this->value;
}

double Vector2D::cos() {
    return this->y/this->value;
}

Vector2D Vector2D::getVectorFromDirection(double value, Vector2D direction) {
    return Vector2D(value*direction.sin(), value*direction.cos());
}

Vector2D Vector2D::operator-(const Vector2D rhs) {
    return Vector2D(this->x-rhs.x, this->y-rhs.y);
}

Vector2D Vector2D::operator+(const Vector2D rhs)
{
	return Vector2D(this->x + rhs.x, this->y + rhs.y);
}

double Vector2D::operator*(const Vector2D rhs) {
    return this->x*rhs.x+this->y*rhs.y;
}

Vector2D Vector2D::operator*(const float rhs)
{
	return Vector2D(this->x*rhs, this->y*rhs);
}

Vector2D Vector2D::normalize() {
    return Vector2D(this->x/this->value, this->y/this->value);
}


