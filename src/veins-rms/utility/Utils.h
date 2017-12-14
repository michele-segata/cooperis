//
// Copyright (C) 2017 Agon Memedi <memedi@ccs-labs.org>
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

#pragma once

#include <cmath>

#include <omnetpp.h>

#include <veins/base/utils/Coord.h>

namespace veins {

using namespace omnetpp;

struct Angles {
    // azimuth angle (in radians)
    double phi;
    // elevation angle (in radians)
    double theta;
};

/**
 * Given the position of a node, as well as the position and the orientation of an RMS, returns the phi and theta angles
 * of the node with respect to the RMS
 * @param rms_v1 normalized 3D vector indicating one side of the RMS (used to compute elevation)
 * @param rms_v2 normalized 3D vector indicating the second side of the RMS (together with rms_v1, defines the surface
 * plane of the RMS). v2 is used to compute the azimuth
 * @param rms_vn normalized 3D vector indicating the normal of the surface (cross product between rms_v1 and rms_v2)
 * @param rms_pos position of the RMS
 * @param node position of the node
 * @param rightHanded whether we are considering a classical right handed coordinate system or a left handed one.
 * Set by default to false as omnet is left handed (y values increase moving south of the screen)
 * @return an Angles structure with the theta and phi angles between the node and the RMS
 */

//Angles get_angles(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, const Coord& node, bool rightHanded=false);

Angles spherical_angles(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, const Coord& node);

Coord spherical_point_beam(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, double phi, double theta);


/**
 * Computes the angle between two 3D vectors assuming the lie on the same plane
 * The plane is described using its normal vn
 */
double angle_3d(const Coord& p1, const Coord& p2, const Coord& vn);

/**
 * Computes the projection of a vector on a plane, given its normal
 * TODO: shouldn't there be a plane position as well??? Mabye not because we are reasoning w.r.t. to the origin
 */
Coord projection(const Coord& plane_normal, const Coord& vector);

/**
 * Computes the dot product between two vectors
 */
double dot(const Coord& p1, const Coord& p2);

/**
 * Computes the cross product between two vectors
 */
Coord cross(const Coord& p1, const Coord& p2);

/**
 * Rotates a vector around a generic axis defined by a 3D unit vector
 * @param v vector to rotate
 * @param axis unit vector indicating the rotation axis
 * @param angle rotation angle in radians TODO: direction???
 * @return the rotated vector
 */
Coord rotate(const Coord& v, const Coord& axis, double angle);

/**
 * Finds the scaling factor alpha of a direction vector representing the point of intersection with a given plane
 * i.e., alpha * vectorDir + starting position = intersection point
 * @param planeN normal vector of the plane
 * @param rmsPos starting position of the vector
 * @param planePoint a point on the plane, necessary to identify the plane together with its normal
 * @param vectorDir direction of the vector
 * @return the scaling factor alpha so that alpha * vectorDir + rmsPos is the intersection point of the vector
 * with the plane identified by planeN and planePoint. If such intersection does not exist, the function returns 0
 */
double plane_vector_intersection(const Coord& planeN, const Coord& rmsPos, const Coord& planePoint, const Coord& vectorDir);

/**
 * Returns two points describing the segment begin the 2D projection of a beam on the ground. The 3D beam is assumed
 * to start from the RMS and to go in the direction described by the azimuth phi and elevation theta
 * @param v1 v1 vector of the RMS (first vector describing the surface)
 * @param v2 v2 vector of the RMS (second vector describing the surface)
 * @param vn normal vector of the surface (facing direction)
 * @param rmsPos position of the RMS
 * @param phi azimuth angle in radians (rotation around v1). 0 indicates no rotation, phi > 0 counterclockwise rotation, phi < 0 clockwise rotation
 * @param theta elevation angle in radians (rotation around v2). 0 indicates no rotation, theta > 0 rotation towards the sky, phi < 0 rotation towards the ground
 * @param groundHeight height of the ground
 * @return the two points to be used for projection
 */
std::list<Coord> project_beam(const Coord& v1, const Coord& v2, const Coord& vn, const Coord& rmsPos, double phi, double theta, double groundHeight=0);

/*
 * My angle representation
 * M_PI*1/2 = 90     South
 * M_PI*1 = 180      West
 * M_PI*3/2 = 270    North
 * M_PI*2 = 360      East
 */

/*
 * Transform the heading of the vehicle from traci
 * to my angle representation
 */
double traci2myAngle(double angleRad);

/*
 * Transform the heading of the vehicle from my
 * angle representation to traci
 */

double myAngle2traci(double angleRad);

// Reverse the heading of the vehicle
double reverseTraci(double angleRad);

double rad2deg(double angleRad);

double deg2rad(double angleDeg);

double traci2cartesianAngle(double angleRad);

double utilTrunc(double number);

// Based on close() from FWMath.h
bool close(double one, double two);

// Return BER of OOK at the given SNR.
double getOokBer(double snr);

double getOokPdr(double snr, int packetLength);

} // namespace veins
