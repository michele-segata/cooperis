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

#include "veins-rms/utility/Utils.h"
#include "veins/base/utils/Coord.h"

using namespace veins;

namespace veins {

//Angles get_angles(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, const Coord& node, bool rightHanded)
//{
//
//    Angles angles;
//
//    // move the origin to the center of the RMS
//    Coord node_orig = node - rms_pos;
//
//    // project the node onto the plane given by v2 and vn (for azimuth phi)
//    Coord proj_v2_vn = projection(rms_v1, node_orig);
//    // compute the angle between vn and the projection
//    angles.phi = -angle_3d(proj_v2_vn, rms_vn, rms_v1);
//
//    // rotate vn by phi around v1 (moves vn right below or above node)
//    Coord vnr = rotate(rms_vn, rms_v1, angles.phi);
//    // find the normal to the plane made by v1 and vnr
//    Coord v1_vnr_n = cross(rms_v1, vnr);
//
//    // compute the angle between the node and vnr, on the plane characterized by v1 and vnr for theta
//    angles.theta = angle_3d(node_orig, vnr, v1_vnr_n);
//
//    if (!rightHanded) angles.phi = -angles.phi;
//
//    return angles;
//
//}

Angles spherical_angles(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, const Coord& node) {

    Angles angles;

    Coord node_orig = node - rms_pos;

    // project the node position onto the plane of the metasurface to compute the azimuth phi
    Coord node_surface_proj = projection(rms_vn, node_orig);
    // to get the azimuth angle, measure the angle between the inverse of v1 and the projected node
    Coord nrms_v1 = {0, 0, 0};
    nrms_v1 -= rms_v1;
    double phi = angle_3d(node_surface_proj, nrms_v1, rms_vn);

    // invert rotation of the node by phi, so that the node is below vn
    Coord derotated = rotate(node_orig, rms_vn, phi);

    // measure the angle between the de-rotated node and vn
    double theta = angle_3d(derotated, rms_vn, rms_v2);

    angles.phi = phi;
    angles.theta = -theta;

    return angles;

}

Coord spherical_point_beam(const Coord& rms_v1, const Coord& rms_v2, const Coord& rms_vn, const Coord& rms_pos, double phi, double theta)
{
    Coord byElevation = rotate(rms_vn, rms_v2, theta);
    return rotate(byElevation, rms_vn, -phi);
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

Coord rotate(const Coord& v, const Coord& axis, double angle)
{
    Coord v_cos_angle = v * cos(angle);
    Coord k_cross_v = cross(axis, v);
    k_cross_v *= sin(angle);
    Coord k_k_v = axis * dot(axis, v) * (1 - cos(angle));
    return v_cos_angle + k_cross_v + k_k_v;
}

double plane_vector_intersection(const Coord& planeN, const Coord& rmsPos, const Coord& planePoint, const Coord& vectorDir)
{
    double Nd = dot(planeN, vectorDir);
    if (std::abs(Nd) < 1e-6)
        // no intersection
        return 0;
    else
        return -dot(planeN, rmsPos - planePoint) / Nd;
}

std::list<Coord> project_beam(const Coord& v1, const Coord& v2, const Coord& vn, const Coord& rmsPos, double phi, double theta, double groundHeight)
{
    static Coord groundN(0, 0, 1);
    static Coord groundPoint(1, 1, groundHeight);
//    // invert angles to perform rotation according to the convention
//    phi = -phi;
//    theta = -theta;
//    // start with vector normal to the surface
//    Coord beamDir(vn);
//    // rotate by phi around v1
//    beamDir = rotate(beamDir, v1, phi);
//    // rotate v2 by phi around v1 as well
//    Coord rotated_v2(v2);
//    rotated_v2 = rotate(rotated_v2, v1, phi);
//    // rotate the vector by theta around new v2
//    beamDir = rotate(beamDir, rotated_v2, theta);
    Coord beamDir = spherical_point_beam(v1, v2, vn, rmsPos, phi, theta);

    Coord groundIntersection = beamDir * plane_vector_intersection(groundN, rmsPos, groundPoint, beamDir) + rmsPos;
    Coord rmsGroundPos(rmsPos);
    rmsGroundPos.z = groundHeight;
    std::list<Coord> beam;
    beam.push_back(rmsGroundPos);
    beam.push_back(groundIntersection);
    return beam;
}

double traci2myAngle(double angleRad)
{
    if (angleRad >= -M_PI && angleRad < 0)
        angleRad = -angleRad;
    else if (angleRad >= 0 && angleRad <= M_PI)
        angleRad = 2 * M_PI - angleRad;
    else
        throw cRuntimeError("unexpected angle");

    return angleRad;
}

double myAngle2traci(double angleRad)
{
    if (angleRad > 0 && angleRad <= M_PI)
        angleRad = -angleRad;
    else if (angleRad > M_PI && angleRad <= 2 * M_PI)
        angleRad = 2 * M_PI - angleRad;
    else if (angleRad == 0)
        angleRad = 0;
    else
        throw cRuntimeError("unexpected angle");

    return angleRad;
}

double reverseTraci(double angle)
{
    // angle is in radians

    // no need to normalize because traci always returns between -180 and 180
    if (angle < 0)
        angle += M_PI;
    else if (angle >= 0)
        angle -= M_PI;
    else
        throw cRuntimeError("unexpected angle");

    return angle;
}

double rad2deg(double angleRad)
{
    return angleRad * 180.0 / M_PI;
}

double deg2rad(double angleDeg)
{
    return angleDeg * M_PI / 180.0;
}

double traci2cartesianAngle(double angleRad)
{

    double angle = angleRad;

    if (angle >= -M_PI && angle < 0)
        angle += 2 * M_PI;
    else if (angle >= 0 && angle < M_PI)
        angle = angle;
    else
        throw cRuntimeError("unexpected angle");

    return angle;
}

double utilTrunc(double number)
{
    // keep only up to third decimal
    return std::trunc(number * 1000.0) / 1000.0;
}

/*
 * same as FWMath::close(), except EPSILON
 */
bool close(double one, double two)
{
    return fabs(one - two) < 0.0000001;
}

double getOokBer(double snr)
{
    // Modelling OOK based on BPSK of NistErrorRate;
    // Assuming Q-func(sqrt(snr))
    double z = std::sqrt(snr / 2.0);
    double ber = 0.5 * erfc(z);
    return ber;

    // If BER is zero we can receive, else not.
    //    if (ber == 0.0)
    //        return 1.0;
    //    else
    //        return 0.0;
}

double getOokPdr(double snr, int packetLength)
{
    double ber = getOokBer(snr);

    if (ber == 0.0) {
        return 1.0;
    }

    return std::pow(1 - ber, (double) packetLength);
}

} // namespace veins
