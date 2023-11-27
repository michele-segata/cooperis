//
// Copyright (C) 2022 Michele Segata <segata@ccs-labs.org>
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
// Base files derived from Veins VLC by Agon Memedi and contributors

#pragma once

#include <cmath>
#include <list>
#include <string>
#include <sstream>

using namespace std;

#ifndef STANDALONE

#include <veins/base/utils/Coord.h>
namespace veins {

#else

class Coord {
public:
    double x;
    double y;
    double z;
    Coord()
        : x(0)
        , y(0)
        , z(0) {}
    Coord(const Coord& p)
        : x(p.x)
        , y(p.y)
        , z(p.z) {}
    Coord(double x, double y, double z)
        : x(x)
        , y(y)
        , z(z) {}
    friend Coord operator*(const Coord& a, double v)
    {
        Coord n(a);
        n *= v;
        return n;
    }
    Coord& operator+=(const Coord& a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }
    Coord& operator-=(const Coord& a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }
    friend Coord operator+(const Coord& a, const Coord& b)
    {
        Coord tmp(a);
        tmp += b;
        return tmp;
    }
    friend Coord operator-(const Coord& a, const Coord& b)
    {
        Coord tmp(a);
        tmp -= b;
        return tmp;
    }
    Coord inverse() const
    {
        Coord tmp;
        tmp.x = -x;
        tmp.y = -y;
        tmp.z = -z;
        return tmp;
    }
    Coord& operator*=(double v)
    {
        x *= v;
        y *= v;
        z *= v;
        return *this;
    }
    double length()
    {
        return sqrt(x*x + y*y + z*z);
    }
    void normalize()
    {
        double l = length();
        x /= l;
        y /= l;
        z /= l;
    }
    string toString() const
    {
        stringstream ss;
        ss << "x=" << x << " y=" << y << " z=" << z;
        return ss.str();
    }
}; // namespace veins

#endif

struct Angles {
    // azimuth angle (in radians)
    double phi;
    // elevation angle (in radians)
    double theta;
};

/**
 * Given the position of a node, as well as the position and the orientation of an RIS, returns the phi and theta angles
 * of the node with respect to the RIS
 * @param ris_v1 normalized 3D vector indicating one side of the RIS (used to compute elevation)
 * @param ris_v2 normalized 3D vector indicating the second side of the RIS (together with ris_v1, defines the surface
 * plane of the RIS). v2 is used to compute the azimuth
 * @param ris_vn normalized 3D vector indicating the normal of the surface (cross product between ris_v1 and ris_v2)
 * @param ris_pos position of the RIS
 * @param node position of the node
 * @return an Angles structure with the theta and phi angles between the node and the RIS
 */
Angles spherical_angles(const Coord& ris_v1, const Coord& ris_v2, const Coord& ris_vn, const Coord& ris_pos, const Coord& node);

/**
 * Given the position and the orientation of an RIS, phi and theta angles, returns the vector in the direction of the given angles
 * @param ris_v1 normalized 3D vector indicating one side of the RIS (used to compute elevation)
 * @param ris_v2 normalized 3D vector indicating the second side of the RIS (together with ris_v1, defines the surface
 * plane of the RIS). v2 is used to compute the azimuth
 * @param ris_vn normalized 3D vector indicating the normal of the surface (cross product between ris_v1 and ris_v2)
 * @param ris_pos position of the RIS
 * @param phi azimuth
 * @param theta elevation
 * @return a vector indicating the direction of (phi, theta) w.r.t. the RIS
 */
Coord spherical_point_beam(const Coord& ris_v1, const Coord& ris_v2, const Coord& ris_vn, const Coord& ris_pos, double phi, double theta);


/**
 * Computes the angle between two 3D vectors assuming they lie on the same plane
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
 * @param angle rotation angle in radians
 * @param left_handed whether considering a left handed coordinate system or not
 * In a left handed system, rotation is clockwise, in a right handed system, counterclockwise
 * @return the rotated vector
 */
Coord rotate(const Coord& v, const Coord& axis, double angle, bool left_handed= true);

/**
 * Finds the scaling factor alpha of a direction vector representing the point of intersection with a given plane
 * i.e., alpha * vectorDir + starting position = intersection point
 * @param planeN normal vector of the plane
 * @param risPos starting position of the vector
 * @param planePoint a point on the plane, necessary to identify the plane together with its normal
 * @param vectorDir direction of the vector
 * @return the scaling factor alpha so that alpha * vectorDir + risPos is the intersection point of the vector
 * with the plane identified by planeN and planePoint. If such intersection does not exist, the function returns 0
 */
double plane_vector_intersection(const Coord& planeN, const Coord& risPos, const Coord& planePoint, const Coord& vectorDir);

/**
 * Returns two points describing the segment begin the 2D projection of a beam on the ground. The 3D beam is assumed
 * to start from the RIS and to go in the direction described by the azimuth phi and elevation theta
 * @param v1 v1 vector of the RIS (first vector describing the surface)
 * @param v2 v2 vector of the RIS (second vector describing the surface)
 * @param vn normal vector of the surface (facing direction)
 * @param risPos position of the RIS
 * @param phi azimuth angle in radians (rotation around v1). 0 indicates no rotation, phi > 0 counterclockwise rotation, phi < 0 clockwise rotation
 * @param theta elevation angle in radians (rotation around v2). 0 indicates no rotation, theta > 0 rotation towards the sky, phi < 0 rotation towards the ground
 * @param groundHeight height of the ground
 * @return the two points to be used for projection
 */
std::list<Coord> project_beam(const Coord& v1, const Coord& v2, const Coord& vn, const Coord& risPos, double phi, double theta, double groundHeight= 0);

#ifndef STANDALONE
} // namespace veins
#endif
