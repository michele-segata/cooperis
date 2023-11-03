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

#include <cmath>

#ifndef STANDALONE

#include "veins/base/utils/Coord.h"
#include "veins-ris/utility/Utils.h"

using namespace veins;

namespace veins {

#else

#include "Utils.h"

#endif

Angles spherical_angles(const Coord& ris_v1, const Coord& ris_v2, const Coord& ris_vn, const Coord& ris_pos, const Coord& node)
{

    Angles angles;

    Coord node_orig = node - ris_pos;

    // project the node position onto the plane of the metasurface to compute the azimuth phi
    Coord node_surface_proj = projection(ris_vn, node_orig);
    // to get the azimuth angle, measure the angle between the inverse of v2 and the projected node
    Coord nris_v2 = {0, 0, 0};
    nris_v2 -= ris_v2;
    double phi = angle_3d(node_surface_proj, nris_v2, ris_vn);

    // invert rotation of the node by phi, so that the node is in the direction of -v2
    Coord derotated = rotate(node_orig, ris_vn, -phi);

    // measure the angle between vn and the de-rotated node
    double theta = angle_3d(ris_vn, derotated, ris_v1);

    angles.phi = phi;
    angles.theta = theta;

    return angles;

}

Coord spherical_point_beam(const Coord& ris_v1, const Coord& ris_v2, const Coord& ris_vn, const Coord& ris_pos, double phi, double theta)
{
    Coord byElevation = rotate(ris_vn, ris_v1, -theta);
    return rotate(byElevation, ris_vn, phi);
}

double angle_3d(const Coord& p1, const Coord& p2, const Coord& vn)
{
    Coord cp = cross(p1, p2);
    return atan2(dot(cp, vn), dot(p1, p2));
}

Coord projection(const Coord& plane_normal, const Coord& vector)
{
    // projection = v - (v^T * n) * n
    double dot_prod = dot(plane_normal, vector);
    Coord tmp = {plane_normal.x * dot_prod, plane_normal.y * dot_prod, plane_normal.z * dot_prod};
    Coord proj = {vector.x - tmp.x, vector.y - tmp.y, vector.z - tmp.z};
    return proj;
}

double dot(const Coord& p1, const Coord& p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Coord cross(const Coord& p1, const Coord& p2)
{
    Coord c;
    c.x = p1.y * p2.z - p1.z * p2.y;
    c.y = -(p1.x * p2.z - p1.z * p2.x);
    c.z = p1.x * p2.y - p1.y * p2.x;
    return c;
}

Coord rotate(const Coord& v, const Coord& axis, double angle, bool left_handed)
{
    if (left_handed)
        angle = -angle;
    Coord v_cos_angle = v * cos(angle);
    Coord k_cross_v = cross(axis, v);
    k_cross_v *= sin(angle);
    Coord k_k_v = axis * dot(axis, v) * (1 - cos(angle));
    return v_cos_angle + k_cross_v + k_k_v;
}

double plane_vector_intersection(const Coord& planeN, const Coord& risPos, const Coord& planePoint, const Coord& vectorDir)
{
    double Nd = dot(planeN, vectorDir);
    if (std::abs(Nd) < 1e-6)
        // no intersection
        return 0;
    else
        return -dot(planeN, risPos - planePoint) / Nd;
}

std::list<Coord> project_beam(const Coord& v1, const Coord& v2, const Coord& vn, const Coord& risPos, double phi, double theta, double groundHeight)
{
    static Coord groundN(0, 0, 1);
    static Coord groundPoint(1, 1, groundHeight);
    Coord beamDir = spherical_point_beam(v1, v2, vn, risPos, phi, theta);

    double beamLength = plane_vector_intersection(groundN, risPos, groundPoint, beamDir);
    // if the beam is parallel to the ground (no intersection) simply draw a very long beam
    if (beamLength < 1e-9)
        beamLength = 1e6;
    Coord groundIntersection = beamDir * beamLength + risPos;
    Coord risGroundPos(risPos);
    risGroundPos.z = groundHeight;
    std::list<Coord> beam;
    beam.push_back(risGroundPos);
    beam.push_back(groundIntersection);
    return beam;
}

#ifndef STANDALONE
} // namespace veins
#endif
