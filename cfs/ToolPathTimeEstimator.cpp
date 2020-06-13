#include "stdafx.h"
#include "ToolPathTimeEstimator.hpp"
#define π 3.141592653


ToolPathTimeEstimator::ToolPathTimeEstimator() {

	/*
    this->nominal_speed = Vector2D(0, 0);
    this->blocks.blockHead = new Block(0.0, 0.0);
    this->blocks.blockTail = this->blocks.blockHead;
    this->acceleration = Vector2D(200, 200);
	this->acceleration.direction = false;
	this->nominal_speed = Vector2D(3000.0 / 60.0, 3000.0 / 60.0);
	this->nominal_speed.direction = false;
	this->jerk = Vector2D(1.75,1.75);
	this->jerk.direction = false;
	this->length = 0.0;
	*/

	this->blocks.blockHead = new Block(0.0, 0.0);
	this->blocks.blockTail = this->blocks.blockHead;
	this->acceleration = Vector2D(250, 250);
	this->acceleration.direction = false;
	this->nominal_speed = Vector2D(1000.0 / 60.0, 1000.0 / 60.0);
	this->nominal_speed.direction = false;
	this->jerk = Vector2D(2.0, 2.0);
	this->jerk.direction = false;
	this->length = 0.0;

}

ToolPathTimeEstimator::~ToolPathTimeEstimator(){
}

void ToolPathTimeEstimator::prepare() {
	this->blocks.blockHead = this->blocks.blockHead->next;
	this->addBlock(0, 0);
}

void ToolPathTimeEstimator::addJump(double x, double y)
{
	if ((x == this->blocks.blockTail->x) && (y == this->blocks.blockTail->y)) {
		return;
	}
	Block* T = new Block(x, y);
	T->previous = this->blocks.blockTail;
	this->blocks.blockTail->next = T;
	T->length = Vector2D(x - T->previous->x, y - T->previous->y);
	T->nominal_speed = Vector2D(JUMPSPEED, JUMPSPEED);
	T->jump = true;
	this->blocks.blockTail = T;
}

void ToolPathTimeEstimator::addBlock(double x, double y) {
	if ((x == this->blocks.blockTail->x) && (y == this->blocks.blockTail->y)) {
		return;}
    
	Block* T = new Block(x,y);

    T->previous = this->blocks.blockTail;
    this->blocks.blockTail->next = T;
    T->length = Vector2D(x-T->previous->x, y-T->previous->y);
	T->nominal_speed = this->nominal_speed;
	T->jump = false;
	this->blocks.blockTail = T;
}

void ToolPathTimeEstimator::addBlocks(Vector2d1 blocks)
{
	for (int i = 0; i < blocks.size(); i++){
		if (i == 0){
			this->blocks.blockHead->x = blocks[i].x;
			this->blocks.blockHead->y = blocks[i].y;
		}
		else{
			double x = blocks[i].x;
			double y = blocks[i].y;
			
			Block* T = new Block(x, y);
			T->previous = this->blocks.blockTail;
			this->blocks.blockTail->next = T;
			T->length = Vector2D(x - T->previous->x, y - T->previous->y);
			T->nominal_speed = this->nominal_speed;
			T->jump = false;
			this->blocks.blockTail = T;
		}
	}
}

int nbb = 0;

double ToolPathTimeEstimator::calculate() {

    Block* blockPtr;
    //Processing the 1..n-1 block
	blockPtr = this->blocks.blockHead->next;

    this->blocks.blockHead->exit_speed = Vector2D(0, 0);
    while (blockPtr != this->blocks.blockTail) {
        blockPtr->entry_speed = Vector2D::getVectorFromDirection(blockPtr->previous->exit_speed.value, blockPtr->length);
        blockPtr->exit_speed = calcJunctionSpeed(blockPtr, blockPtr->next);
        if (blockPtr->exit_speed.value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, blockPtr->length).value) {//"-"????
            allow_speed = Vector2D::getVectorFromDirection(calcAllowedSpeed(acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value, blockPtr->length);
            if (allow_speed.value < blockPtr->entry_speed.value) {
                backwardIterate(blockPtr);
			}
			else {
				forwardIterate(blockPtr);
			}
        }
        blockPtr->build_time = calcBlockBuildTime(blockPtr);
        blockPtr = blockPtr->next;
    }
    
    //Processing the n block
    blockPtr->entry_speed = blockPtr->previous->exit_speed;

    allow_speed = Vector2D::getVectorFromDirection(sqrt(2*acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value), blockPtr->length);
    blockPtr->exit_speed = Vector2D(0, 0);
	if (blockPtr->exit_speed.value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, blockPtr->length).value) {//"-"????
		allow_speed = Vector2D::getVectorFromDirection(sqrt(2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value + pow(blockPtr->exit_speed.value, 2)), blockPtr->length);
		if (allow_speed.value < blockPtr->entry_speed.value) {
			backwardIterate(blockPtr);
		}
		else {
			forwardIterate(blockPtr);
		}
	}

    blockPtr->build_time = calcBlockBuildTime(blockPtr);
    


    //Sum up the build time
    double buildTime = 0;
	blockPtr = this->blocks.blockHead->next;


	/***************************************************************/

	int iter = 0;
    while (blockPtr != nullptr) {
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		buildTime += blockPtr->build_time;
		length += blockPtr->length.value;
		iter++;
		blockPtr = blockPtr->next;
	}
	/***************************************************************/

    return buildTime;
}

void ToolPathTimeEstimator::outputSpeed_OBJ(std::string path, Vector3d1 &cnc_3d_path)
{
	std::ofstream out(path);
	
	out << cnc_3d_path.size() << std::endl;

	Block* blockPtr = this->blocks.blockHead->next;
	int iter = 0;

	out << 0.0 << " " << cnc_3d_path[iter][0] << " " << cnc_3d_path[iter][1] << " " << cnc_3d_path[iter][2] << std::endl;

	while (blockPtr != nullptr)
	{
		iter++;
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		double v = blockPtr->length.get_value() / blockPtr->build_time;
		out << v << " " << cnc_3d_path[iter][0] << " " << cnc_3d_path[iter][1] << " " << cnc_3d_path[iter][2] << std::endl;
		blockPtr = blockPtr->next;
	}
	out.clear();
	out.close();

}


void ToolPathTimeEstimator::outputSpeed_WAVE(std::string path, Vector3d1 &cnc_3d_path, double time_delta)
{
	std::ofstream out(path);

	double recordTime = time_delta;
	int nb = 0;
	double buildTime = 0.0;
	Block* blockPtr = this->blocks.blockHead->next;
	int iter = 0;

	out << "i   t   v" << std::endl;
	while (blockPtr != nullptr) {

		blockPtr->build_time = calcBlockBuildTime(blockPtr);

		while (recordTime>buildTime&&recordTime <= buildTime + blockPtr->build_time){
			nb++;
			double v = blockPtr->getSpeed(recordTime - buildTime);
			out <<iter<<" "<< recordTime << " " << v << std::endl;
			recordTime += time_delta;
		}

		buildTime += blockPtr->build_time;
		//length += blockPtr->length.value;
		blockPtr = blockPtr->next;

		iter++;
	}
	out.clear();
	out.close();
}

double ToolPathTimeEstimator::calcBlockBuildTime(Block* blockPtr) {
    if (blockPtr->length.value <= ((2*pow(blockPtr->nominal_speed.maxAlongDirection(blockPtr->length),2)-pow(blockPtr->entry_speed.value,2)-pow(blockPtr->exit_speed.value,2))/2/acceleration.maxAlongDirection(blockPtr->length))) 
	{
        blockPtr->max_speed = Vector2D::getVectorFromDirection(sqrt((2*acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value+pow(blockPtr->exit_speed.value, 2)+pow(blockPtr->entry_speed.value, 2))/2), blockPtr->length);
		blockPtr->s1 = (blockPtr->entry_speed - blockPtr->max_speed).get_value() / acceleration.maxAlongDirection(blockPtr->direction);
		blockPtr->s2 = 0;
		blockPtr->s3 = (blockPtr->exit_speed - blockPtr->max_speed).get_value() / acceleration.maxAlongDirection(blockPtr->direction);

		if (blockPtr->s3 < 0) {
			blockPtr->s3 = 0;
		}

		if (blockPtr->s1 < 0) {
			blockPtr->s1 = 0;
		}

        return (2*blockPtr->max_speed.value-blockPtr->entry_speed.value-blockPtr->exit_speed.value)/acceleration.maxAlongDirection(blockPtr->length);
        
    } else {

		double time = (2 * blockPtr->nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->entry_speed.value - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length) 
			+ (blockPtr->length.value - (2 * pow(blockPtr->nominal_speed.maxAlongDirection(blockPtr->length), 2) - pow(blockPtr->entry_speed.value, 2) - pow(blockPtr->exit_speed.value, 2)) 
			/ 2 / acceleration.maxAlongDirection(blockPtr->length)) / blockPtr->nominal_speed.maxAlongDirection(blockPtr->length);
	
		blockPtr->max_speed = Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length);
		blockPtr->s1 = (blockPtr->entry_speed - blockPtr->max_speed).get_value() / acceleration.maxAlongDirection(blockPtr->direction);
		blockPtr->s3 = (blockPtr->exit_speed - blockPtr->max_speed).get_value() / acceleration.maxAlongDirection(blockPtr->direction);
		blockPtr->s2 = time - blockPtr->s1 - blockPtr->s3;


		if (blockPtr->s3 < 0) {
			blockPtr->s3 = 0;
		}

		if (blockPtr->s2 < 0) {
			blockPtr->s2 = 0;
		}
		
		if (blockPtr->s1 < 0) {
			blockPtr->s1 = 0;
		}
		return time;
    }
}

void ToolPathTimeEstimator::backwardIterate(Block *blockPtr) {
    if (calcJunctionSpeed(blockPtr->previous, blockPtr).value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value) {
        allow_speed = Vector2D::getVectorFromDirection(sqrt(2*acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value+pow(calcJunctionSpeed(blockPtr->previous, blockPtr).value, 2)), blockPtr->length);
        if (allow_speed.value > Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length).value) {
            allow_speed = Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length);
        }
        if (blockPtr->exit_speed.value<=allow_speed.value) {
            return;
        } else {
            forwardIterate(blockPtr);
        }
    } else {
        blockPtr->entry_speed = allow_speed;
        blockPtr = blockPtr->previous;
        blockPtr->exit_speed = Vector2D::getVectorFromDirection(allow_speed.value, blockPtr->length);
        allow_speed = Vector2D::getVectorFromDirection(calcAllowedSpeed(acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value, blockPtr->length);
        if (blockPtr->entry_speed.value<=allow_speed.value) {
			return;
        }
        backwardIterate(blockPtr);
    }
}

void ToolPathTimeEstimator::forwardIterate(Block *blockPtr) {
    
    blockPtr->exit_speed = allow_speed;
    blockPtr = blockPtr->next;
    blockPtr->entry_speed = allow_speed;
    allow_speed = Vector2D::getVectorFromDirection(sqrt(2*acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value+pow(blockPtr->entry_speed.value, 2)), blockPtr->length);
    if (allow_speed.value > Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length).value) {
        allow_speed = Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length);
    }
    if (blockPtr->exit_speed.value<=allow_speed.value) {
        return;
    }
    forwardIterate(blockPtr);
    
}

Vector2D ToolPathTimeEstimator::calcAllowedSpeed(double acceleration, double initialSpeed, Vector2D length)
{
	return Vector2D::getVectorFromDirection(sqrt(2*acceleration*length.value+initialSpeed*initialSpeed), length);
}

void ToolPathTimeEstimator::updateNominalSpeed()
{
	Block * blockPtr;
	blockPtr = this->blocks.blockHead->next;
	while (blockPtr->next != this->blocks.blockTail) {
		if (blockPtr->jump) {
			blockPtr->nominal_speed = Vector2D(JUMPSPEED, JUMPSPEED);
		}
		else {
			blockPtr->nominal_speed = this->nominal_speed;
		}
		blockPtr = blockPtr->next;
	}
}



Vector2D ToolPathTimeEstimator::calcJunctionSpeed(Block *A, Block *B) {
	if (A == nullptr) {
        A = new Block(0,0);
    }
	Vector2D A_speed, junctionSpeed;
	double arc_cos, arc_cos_2;
	A_speed = calcAllowedSpeed(acceleration.maxAlongDirection(A->length), A->entry_speed.value, A->length);
	arc_cos = (A_speed*B->length) / (A_speed.value*B->length.value);
	if (abs(arc_cos - 1) <= 1e-10) {
		if (A_speed.value >= nominal_speed.maxAlongDirection(A->length)) {
			return Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(A->length), A->length);
		}
		return A_speed;
	}
	arc_cos_2 = sqrt((1 + arc_cos) / 2);
	junctionSpeed = Vector2D::getVectorFromDirection(jerk.maxAlongDirection(B->length.normalize()-A->length.normalize()) / (2 * sqrt(1 - pow(arc_cos_2, 2))), A->length);
	if (junctionSpeed.value >= nominal_speed.maxAlongDirection(A->length)) {
		return Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(A->length), A->length);
	}
	return junctionSpeed;
}

void ToolPathTimeEstimator::detail(FILE * out) {
	//Sum up the build time
	double buildTime = 0;
	double length = 0;
	Block * blockPtr;
	blockPtr = this->blocks.blockHead->next;
	while (blockPtr->next != nullptr) {
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		if (abs(blockPtr->entry_speed.value - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length) >= blockPtr->build_time + 1e-4) {
			//can't be small,only equal, this is the minimum time, add a epsilon in case the error
			fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, blockPtr->entry_speed.value, blockPtr->exit_speed.value, blockPtr->length.value, blockPtr->build_time);
		}
		else {
			if (((blockPtr->nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->entry_speed.value) + (blockPtr->nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->exit_speed.value)) /
				acceleration.maxAlongDirection(blockPtr->length) >= blockPtr->build_time + 1e-4) {
				//can not reach the nominal speed or just reach the nominal speed in a very short time
				if (blockPtr->entry_speed.value - blockPtr->exit_speed.value <= 1e-4) {
					double mid_x, mid_y, mid_speed;
					mid_x = (blockPtr->x + blockPtr->previous->x) / 2;
					mid_y = (blockPtr->y + blockPtr->previous->y) / 2;
					mid_speed = (blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, blockPtr->length.value / 2, blockPtr->build_time / 2);
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value / 2, blockPtr->build_time / 2);
				}
				else {
					if (blockPtr->entry_speed.value > blockPtr->exit_speed.value) {
						double mid_x, mid_y, mid_speed, mid_time;
						Vector2D travel;
						mid_speed = (blockPtr->build_time*acceleration.maxAlongDirection(blockPtr->length) + blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
						travel = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
						mid_x = blockPtr->previous->x + travel.x;
						mid_y = blockPtr->previous->y + travel.y;
						mid_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, travel.value, mid_time);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value - travel.value, blockPtr->build_time - mid_time);
					}
					else {
						double mid_x, mid_y, mid_speed, mid_time;
						Vector2D travel;
						mid_speed = (blockPtr->build_time*acceleration.maxAlongDirection(blockPtr->length) + blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
						travel = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
						mid_x = blockPtr->previous->x + travel.x;
						mid_y = blockPtr->previous->y + travel.y;
						mid_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, travel.value, mid_time);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value - travel.value, blockPtr->build_time - mid_time);
					}
				}
			}
			else {
				//Must reach the maximum speed and move with this speed
				double mid1_x, mid1_y, mid_speed, mid1_time;
				double mid2_x, mid2_y, mid2_time_from_exit;
				Vector2D travel_acc, travel_curise;
				mid_speed = nominal_speed.maxAlongDirection(blockPtr->length);
				travel_acc = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
				mid1_x = blockPtr->previous->x + travel_acc.x;
				mid1_y = blockPtr->previous->y + travel_acc.y;
				mid1_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
				mid2_time_from_exit = (mid_speed - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
				travel_curise = Vector2D::getVectorFromDirection(mid_speed*(blockPtr->build_time - mid1_time - mid2_time_from_exit), blockPtr->length);
				mid2_x = mid1_x + travel_curise.x;
				mid2_y = mid1_y + travel_curise.y;
				if (blockPtr->entry_speed.value - mid_speed > 1e-5) {
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid1_x, mid1_y, blockPtr->entry_speed.value, mid_speed, travel_acc.value, mid1_time);
				}
				
				fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid2_x, mid2_y, blockPtr->entry_speed.value, mid_speed, travel_curise.value,
					blockPtr->build_time - mid1_time - mid2_time_from_exit);
				if (mid_speed - blockPtr->exit_speed.value > 1e-5) {
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value,
						(blockPtr->exit_speed.value + mid_speed)*(blockPtr->build_time - mid1_time - mid2_time_from_exit) / 2, mid2_time_from_exit);
				}
			}
		}
			buildTime += blockPtr->build_time;
			length += blockPtr->length.value;
			blockPtr = blockPtr->next;
	}
}

void ToolPathTimeEstimator::timeline(FILE * out, float resolution)
{
	float time_counter = 0;
	float tick = resolution;
	double s1, s2, s3;
	Block * blockPtr;
	Vector2D pre_location;
	Vector2D pre_speed;
	blockPtr = this->blocks.blockHead;
	s1 = blockPtr->s1;
	s2 = blockPtr->s2;
	s3 = blockPtr->s3;
	blockPtr->build_time = this->calcBlockBuildTime(blockPtr);
	while (blockPtr != this->blocks.blockTail) {
		if (blockPtr->build_time < tick) {
			tick -= blockPtr->build_time;
			pre_location = Vector2D(blockPtr->x, blockPtr->y);
			pre_speed = Vector2D(blockPtr->exit_speed);
			blockPtr = blockPtr->next;
			blockPtr->build_time = this->calcBlockBuildTime(blockPtr);
			s1 = blockPtr->s1;
			s2 = blockPtr->s2;
			s3 = blockPtr->s3;

		}
		else {
			if (s1 < tick) {
				if (s1 == 0) {
					if (s2 < tick) {
						if (s2 == 0) {
							if (s3 < tick) {
								if (s3 == 0) {
									//jump to next segment, almost impossible to reach this step
									pre_location = Vector2D(blockPtr->x, blockPtr->y);
									pre_speed = Vector2D(blockPtr->exit_speed);
									blockPtr = blockPtr->next;
									blockPtr->build_time = this->calcBlockBuildTime(blockPtr);
									s1 = blockPtr->s1;
									s2 = blockPtr->s2;
									s3 = blockPtr->s3;
								}
								else {
									//finish the last stage and jump
									tick -= s3;
									pre_location = Vector2D(blockPtr->x, blockPtr->y);
									pre_speed = Vector2D(blockPtr->exit_speed);
									blockPtr = blockPtr->next;
									blockPtr->build_time = this->calcBlockBuildTime(blockPtr);
									s1 = blockPtr->s1;
									s2 = blockPtr->s2;
									s3 = blockPtr->s3;
								}
							}
							else {
								//finish the tick in last stage
								s3 -= tick;
								auto speed = pre_speed - Vector2D::getVectorFromDirection(acceleration.maxAlongDirection(blockPtr->length)*tick, blockPtr->length);
								auto location = pre_location + Vector2D::getVectorFromDirection((pre_speed.get_value() + speed.get_value())*tick / 2, blockPtr->length);
								fprintf(out, "%f %f %f\n", location.x, location.y, speed.value);
								pre_speed = speed;
								pre_location = location;
								tick = resolution;
							}
						}
						else {
							//finish second stage
							auto speed = pre_speed - Vector2D::getVectorFromDirection(blockPtr->max_speed.value*s2, blockPtr->length);
							auto location = pre_location + Vector2D::getVectorFromDirection((pre_speed.get_value() + speed.get_value())*s2 / 2, blockPtr->length);
							pre_speed = speed;
							pre_location = location;
							tick -= s2;
							s2 = 0;
						}
					}
					else {
						//finish the tick in second stage
						s2 -= tick;
						auto speed = pre_speed - Vector2D::getVectorFromDirection(blockPtr->max_speed.value*tick, blockPtr->length);
						auto location = pre_location + Vector2D::getVectorFromDirection((pre_speed.get_value() + speed.get_value())*tick / 2, blockPtr->length);
						fprintf(out, "%f %f %f\n", location.x, location.y, speed.value);
						pre_speed = speed;
						pre_location = location;
						tick = resolution;
					}
				}
				else {
					//finish the first stage
					auto speed = pre_speed + Vector2D::getVectorFromDirection(acceleration.maxAlongDirection(blockPtr->length)*s1, blockPtr->length);
					auto location = pre_location + Vector2D::getVectorFromDirection((pre_speed.get_value() + speed.get_value())*s1 / 2, blockPtr->length);
					pre_speed = speed;
					pre_location = location;
					tick -= s1;
					s1 = 0;
				}
			}
			else {
				//finish tick in first stage
				s1 -= tick;
				auto speed = pre_speed - Vector2D::getVectorFromDirection(acceleration.maxAlongDirection(blockPtr->length)*tick, blockPtr->length);
				auto location = pre_location + Vector2D::getVectorFromDirection((pre_speed.get_value() + speed.get_value())*tick / 2, blockPtr->length);
				fprintf(out, "%f %f %f\n", location.x, location.y, speed.value);
				pre_speed = speed;
				pre_location = location;
				tick = resolution;

			}
		}
	}
}

double ToolPathTimeEstimator::distance(double acceleration, double small_speed, double large_speed) {
	return (pow(large_speed, 2) - pow(small_speed, 2)) / (2 * acceleration);
}




