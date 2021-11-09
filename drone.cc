

#ifndef DRONE_CLASS
#define DRONE_CLASS
#include <tgmath.h>
#include "drone.h"
#include "vector3.h"
#include <math.h>
#include "Eigen/Dense"
#define SCALE 0.01
#define JOYSTICK_SCALE 0.01
#define MASS 3 //KG
#define PI 3.141592653589
#define THRESHOLD 0.0001
#define GRAVITY -9.8
#define DRAG 0.999
#define MAX_ACCELERATION 10.0
#include <algorithm>    // std::min

//Implement the drone.h

using namespace csci3081;

double magnitude(Eigen::Vector3d vector);
double clamp(double value, double min, double max);
//
//
///*
// * TODO: Optomize memory and variables by using pointers and references, pre-compute cos sin for RotationMatrix
// */
//Vector3 position;
//Vector3 velocity;
//Vector3 euler;


double thrust;
double thrust_adjustment;



Eigen::Matrix<double, 3, 3> rotationMatrix;
Eigen::Vector3d acceleration;
Eigen::Vector3d targetAcceleration;
Eigen::Vector3d velocity;
Eigen::Vector3d position;
Eigen::Vector3d euler;
Eigen::Vector3d target;
Eigen::Vector3d repulsion;
Eigen::Vector3d temp;

//Vector3 target;
/**
 * Goal: to use euler angles to calculate angular velocities and translate that into a velocity vector
 * This will allow us to use other more complex algorithms found in most drone mechanics to make a more advanced drone simulation
 */

// Default Constructor
Drone::Drone() {

//    euler.dot(position);
    euler <<    0,
                0,
                0;

    acceleration << 0,
                    0,
                    thrust;

    //Target to try and move to
    target <<   10,
                -100,
                100;



}
Eigen::Vector3d GetAccelerationVector() {
    return acceleration;
}

// Updates the drone position and other dynamic properties
void Drone::Update(double dt) {
    //Time
    double t = 1;
//    if (length (target - position) < 1) {
//        return;
//    }
    //Calculate the desired acceleration, and then the angles to get that acceleration
    //Eigen::Vector3d targetAcceleration;
    //https://stackoverflow.com/questions/42554960/get-xyz-angles-between-vectors


    //this acceleration is based off of the position, but rather we want to calcualte the acceleration
    //
    //    targetAcceleration <<   ((target(0, 0) - position(0, 0)) - velocity(0, 0)*t) / pow(t, 2),
    //                            ((target(1, 0)) - position(1, 0)  - velocity(1, 0)*t) / pow(t, 2),
    //                            ((target(2, 0) - position(2, 0)) - velocity(2, 0)*t) / pow(t, 2);
//    targetAcceleration << 100, 100, 100;
    //First we need to calculate the seperate forces


    //calculate vector that is equal to max acceleration in magnitude going from start point to end point
        //Rules for the acceleration vector
            //1) the magnitude of the acc vector can not be greater than the distance from drone to endpoint
            //2) If the acceleration vector meets the end of the repulsion force, it will point to the endpoint of the repulsion force rather than the desired endpoint (natural deceleration)


            //Use that accerlation vector to change the angle of the drone
            //Watcht the magic happen

    /*
     * SUDO CODE
     * MAX_ACCELERATION = 50
     *
     * //Get direction of vector
     * vector_to_point = target - position;
     * //Only downscale the acceleration vector, but do not upscale

     *
     * //Calculate repulsion force as well
     * //Assume we have a velocity vector velocity that knows the drones current velocity
     * //Repulsion is the velocity vector pointed in reverse
     *
     * //The repulsion force can only grow, not shrink
     * if (repulsion.magnitude < velocity.magnitude):
     *      repulsion = velocity * -1;
     *      repulsion = repulsion * (velocity.magnitude / velocity.magnitude)
     *
     * //Now how does the repulsion affect the desired_acceleration?
     * //Let repulsion force push the desired_acceleration the other direction
     * //Check if repulsion force should be considered
     * //Only update the repulsion vector if the drone is not inside of the repulsion vector already
     *
     * repulsion_vector = <repulsion_x - position_x, repulsion_y - position_y, repulsion_z - position_z>
     *
     *
     *
     * distance_of_repulsion_force = sqrt( pow(repulsion_x - position_x, 2) + pow(repulsion_y - position_y, 2) + pow(repulsion_z - position_z, 2))
     * if (desired_acceleration.magnitude > distance_of_repulsion_force):
     *      desired_acceleration = repulsion_vector;
     *      //Make sure the repulsion_vector is not greater than the max_acceleration
     *      a = sqrt( (pow(repulsion_x_component, 2) + pow(repulsion_y_component, 2) + pow(repulsion_z_component, 2)) / pow(max_acc, 2))
     *      repulsion_vector_scale_by =  Math.Min(1, a)
     *      repulsion_vector = repulsion_vector * repulsion_vector_scale_by
     *      desired_acceleration = repulsion_vector
     *
     * downscale_vector_to_point_factor = Math.min(sqrt((pow(x_component, 2) + pow(y_component, 2) + pow(z_component, 2)) / pow(max_acc, 2)), 1)
     * desired_acceleration = initial_acceleration_vector * downscale_vector_to_point_factor;
     *
     * velocity.x += desired_acceleration.x * dt;
     * velocity.y += desired_acceleration.y * dt;
     * velocity.z += desired_acceleration.z * dt;
     *
     * position.x += velocity.x * dt;
     * position.y += velocity.y * dt;
     * position.z += velocity.z * dt;
     *
     *
     *
     *
     *
     */
//Now that we have the correct target acceleration we just need to translate that into euler angles!
//NEW METHOD

      //Get direction of vector
//      Eigen::Vector3d vector_to_point = target - position;
      //Only downscale the acceleration vector, but do not upscale
      //Calculate repulsion force as well
      //Assume we have a velocity vector velocity that knows the drones current velocity
      //Repulsion is the velocity vector pointed in reverse


//    if (vector_to_point.squaredNorm() < 1) {
//        return;
//    }
//      //The repulsion force can only grow, not shrink
//      if (repulsion.squaredNorm() * (velocity.squaredNorm() / MAX_ACCELERATION) < velocity.squaredNorm()) {
//          repulsion = velocity * -1;
//          repulsion = repulsion * (velocity.squaredNorm() / MAX_ACCELERATION);
//      }
//
//
//      //Now how does the repulsion affect the desired_acceleration?
//      //Let repulsion force push the desired_acceleration the other direction
//      //Check if repulsion force should be considered
//      //Only update the repulsion vector if the drone is not inside of the repulsion vector already
//      Eigen::Vector3d repulsion_vector(repulsion(0, 0)- position(0, 0), repulsion(1, 0) - position(1, 0), repulsion(2, 0) - position(2, 0));
//
////      double distance_of_repulsion_force = sqrt(pow(repulsion_x - position_x, 2) + pow(repulsion_y - position_y, 2) + pow(repulsion_z - position_z, 2));
////       double scale_down;
//        if (vector_to_point.squaredNorm() < repulsion.squaredNorm()) {
////           desired_acceleration = repulsion_vector;
//           vector_to_point = repulsion;
//           //Make sure the repulsion_vector is not greater than the max_acceleration
//           double scale_down = sqrt((pow(repulsion(0, 0), 2) + pow(repulsion(1, 0), 2) + pow(repulsion(2, 0), 2)) / pow(MAX_ACCELERATION, 2));
//           scale_down =  std::min(1.0, scale_down);
//           vector_to_point = vector_to_point * scale_down;
//       }


//        double scale_down;
//      scale_down = min(sqrt((pow(x_component, 2) + pow(y_component, 2) + pow(z_component, 2)) / pow(max_acc, 2)), 1)
//      desired_acceleration = initial_acceleration_vector * downscale_vector_to_point_factor;

//        if (length(target - position) < 2)
//            return;






        temp << (velocity / MAX_ACCELERATION) / dt;

        if (magnitude(temp) > magnitude(repulsion)) {
            repulsion << temp;
        }

        std::cout << magnitude(target - position) << std::endl;

        if (magnitude(target - position) < 2 && magnitude(velocity) < 2) {
            repulsion.setZero();
            acceleration.setZero();
        } else {
            acceleration = target - position - repulsion;
            acceleration *= MAX_ACCELERATION / magnitude(acceleration);
        }

        velocity += acceleration * dt;

        position += velocity * dt;








//    acceleration << targetAcceleration(0, 0),
//                    targetAcceleration(1, 0),
//                    targetAcceleration(2, 0);
//    //First find the angle between the two vectors using dot product
//    //float angle =   targetAcceleration.dot(position);
//    //printf("%f\n", angle);
//    std::cout << targetAcceleration << std::endl;
//    //Here we increment by a small value. Perform matrix multiplication to get x y and z vectors
//
//    //z = psi
//    //x = phi
//    //y = theta
//
//    //https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec29.pdf
//    //https://core.ac.uk/download/pdf/47181924.pdf
////    rotationMatrix <<   0, 0,     cos(euler(2, 0))*sin(euler(1, 0))+cos(euler(1, 0))*sin(euler(0, 0))*sin(euler(2, 0)),
//
//
////    rotationMatrix <<   cos(euler(2, 0))*cos(euler(1, 0)) - sin(euler(0, 0))*sin(euler(1, 0))*sin(euler(2, 0)),    -cos(euler(0, 0))*sin(euler(2, 0)),      cos(euler(2, 0))*sin(euler(1, 0))+cos(euler(1, 0))*sin(euler(0, 0))*sin(euler(2, 0)),
////                        cos(euler(1, 0))*sin(euler(2, 0)) + cos(euler(2, 0))*sin(euler(0, 0))*sin(euler(1, 0)),     cos(euler(0, 0))*cos(euler(2, 0)),      sin(euler(2, 0))*sin(euler(1, 0)) - cos(euler(1, 0))*sin(euler(0, 0))*cos(euler(2, 0)),
////                        -cos(euler(0, 0))*sin(euler(1, 0)),                                             sin(euler(0, 0)),                   cos(euler(0, 0))*cos(euler(1, 0));
////
////    //Our starting vector is a global fixed coordinate vector with size of thrust
////    //If we want to make it relative we pass in the previously changed vector into this original, which will euler.* angular velocity instead of angles
////
////    acceleration << 0,
////                    0,
////                    thrust;
////////
//////    //Rotate our vector by the given euler values using Z X Y Euler angles
////    acceleration = rotationMatrix * acceleration;
//
//    /*
//     *Gravity effects velocity in a squared fashion from the equation  y= at^2 + v*x + p;
//     *Can we automatically increase thrust to try and maintain the current height?
//     *Yes, we know the z component needs to be 9.8, so we simply need to multiply the vector by amount t = original * (target/original)
//     *
//     */
//
////    double thrust_adjustment = thrust/acceleration(2, 0);
////    acceleration = acceleration * thrust_adjustment;
//
//
//    velocity(0, 0) += acceleration(0, 0)*dt;
//    velocity(1, 0) += acceleration(1, 0)*dt;
//    velocity(2, 0) += (MASS*GRAVITY + acceleration(2, 0))*dt;
//
////    std::cout << "Thrust:" << acceleration(2, 0) << ", rx:" << euler(0, 0) << " ry:" << euler(1, 0) << " rz:" << euler(2, 0) << std::endl;
//
//    //Provide drag to each of the acc. vectors bringing it back closer to zero
//    velocity(0, 0) *= DRAG;
//    velocity(1, 0) *= DRAG;
//    velocity(2, 0) *= DRAG;
//
//
//
//    position(0, 0) += velocity(0, 0)*dt;
//    position(1, 0) += velocity(1, 0)*dt;
//    position(2, 0) += velocity(2, 0)*dt;
}

double magnitude(Eigen::Vector3d vector) {
    return sqrt( pow(vector(0, 0), 2) + pow(vector(1, 0), 2) + pow(vector(2, 0), 2));
}
// Gets the drone position.  Index represents the 3D spatial component
double Drone::GetPosition(int index) {
    switch (index) {
        case 0:
            return position(1, 0) * SCALE;
        case 1:
            return position(0, 0) * SCALE;
        case 2:
            return position(2, 0) * SCALE;
        default:
            return 0.0;
    }
}

// Gets the propeller speed.  The index represents the propeller number 0-4.
// Speeds between 0-2 work well, where a speed of 2 is faster than 1.
double Drone::GetPropellerSpeed(int index) {
    double speed = 0;
    switch (index) {
        case 0:
            speed = acceleration(2, 0)/2 + acceleration(1, 0)/4 - acceleration(0, 0)/4;
            speed *= 0.08;
            return speed;
        case 1:
            speed = acceleration(2, 0)/2 - acceleration(1, 0)/4 - acceleration(0, 0)/4;
            speed *= 0.08;
            return speed;
        case 2:
            speed = acceleration(2, 0)/2 + acceleration(1, 0)/4 + acceleration(0, 0)/4;
            speed *= 0.08;
            return speed;
        case 3:
            speed = acceleration(2, 0)/2 - acceleration(1, 0)/4 + acceleration(0, 0)/4;
            speed *= 0.08;
            return speed;
        default:
            return speed;
    }
}

// Sets the direction of the joystick.  For example x = 1, means move in the positive x direction and
// x = -1 means move in the negative x direction.
void Drone::SetJoystick(double x, double y, double z, double rotate) {
    euler(0, 0) -= x*JOYSTICK_SCALE;
    euler(1, 0) += y*JOYSTICK_SCALE;
    euler(2, 0) += rotate*JOYSTICK_SCALE;
    thrust += z*JOYSTICK_SCALE;

    thrust = clamp(thrust, 1, 200);
    euler(0, 0) = clamp(euler(0, 0), -PI/4, PI/4);
    euler(1, 0) = clamp(euler(1, 0), -PI/4, PI/4);
}


double clamp(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}




#endif