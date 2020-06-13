#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "Vector2D.hpp"

#pragma once
#define JUMPSPEED 130
#define DEFAULT_NOMINALSPEED 50
#define MIN_NOMINALSPPED 20
#define MAX_NOMINALSPEED 200

class Block
{
public:
	Block(double x, double y)
	{
		this->next = nullptr;
		this->previous = nullptr;
		this->x = x;
		this->y = y;
		this->length = Vector2D(0, 0);
	}
	~Block();
	double getSpeed(double s)
	{
		if (s < s1) return s / s1*(max_speed.get_value() - entry_speed.get_value()) + entry_speed.get_value();

		double d = s1 + s2;
		double d0 = s1 + s2+s3;

		if (s>s1 + s2) return (build_time - s) / s3*(max_speed.get_value() - exit_speed.get_value()) + exit_speed.get_value();
		return max_speed.get_value();
	}
	double x, y;
	float s1, s2, s3;
	bool jump;
	Vector2D max_speed;
	Vector2D nominal_speed;
	Vector2D entry_speed;
	Vector2D exit_speed;
	double build_time;
	Vector2D length;
	Vector2D direction;
	Block* next;
	Block* previous;
};

class ToolPathTimeEstimator
{
public:
	ToolPathTimeEstimator();
	~ToolPathTimeEstimator();
	double calculate();

	void outputSpeed_WAVE(std::string path, Vector3d1 &cnc_3d_path, double time_delta);
	void outputSpeed_OBJ(std::string path, Vector3d1 &cnc_3d_path);

	void detail(FILE * out);
	void timeline(FILE * out, float resolution);
	void prepare();
	void addJump(double x, double y);
	void addBlock(double x, double y);
	void addBlocks(Vector2d1 blocks);

    Vector2D calcJunctionSpeed(Block* A, Block* B);
    double calcBlockBuildTime(Block* blockPtr);
    void backwardIterate(Block* blockPtr);
    void forwardIterate(Block* blockPtr);
	Vector2D calcAllowedSpeed(double acceleration, double target_speed, Vector2D length);

	double length;
	void updateNominalSpeed();
	Vector2D nominal_speed;
    Vector2D allow_speed;
	Vector2D acceleration;
	Vector2D jerk;
	double distance(double acceleration, double entry_speed, double exit_float);
	struct List{
		Block* blockHead;
		Block* blockTail;
	}blocks;
};