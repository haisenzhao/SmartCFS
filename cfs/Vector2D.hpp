//
//  Vector2D.hpp
//  ToolPathTimeEstimator
//
//  Created by Celr on 15/12/10.
//  Copyright © 2015年 Celr. All rights reserved.
//

#ifndef Vector2D_hpp
#define Vector2D_hpp

#include "stdafx.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>


class Vector2D {
    
public:
    Vector2D();
    Vector2D(double x, double y);
    double x,y;
	bool direction;
    double maxAlongDirection(Vector2D direction);
	double value;
    double get_value();
    double tan();
    double sin();
    double cos();
    Vector2D static getVectorFromDirection(double value, Vector2D direction);
    Vector2D operator-(const Vector2D rhs);
	Vector2D operator+(const Vector2D rhs);
    double operator*(const Vector2D rhs);
	Vector2D operator*(const float rhs);
    Vector2D normalize();
};


#endif /* Vector2D_hpp */
