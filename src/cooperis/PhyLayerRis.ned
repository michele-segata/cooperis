//
// Copyright (C) 2022-2024 Michele Segata <segata@ccs-labs.org>
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

package org.car2x.cooperis;

import org.car2x.veins.base.phyLayer.BasePhyLayer;

simple PhyLayerRis extends BasePhyLayer
{
    parameters:
        @class(veins::PhyLayerRis);

        // available rates are 120, 180, 240, 360, 480, 720, 960, 1080 Mbps
        int bitrate @unit(bps);

        double centerFrequency @unit(Hz);

        double txPower @unit(mW);

        // determine whether this node is an RIS or a UE
        bool isRIS;
        // just for graphically drawing the beams from the RIS
        double nodesAntennaHeight @unit(m) = default(1.895m);
        
        // vectors used to describe the orientation of the RIS
        // considering the RIS as a rectangle, vector v2 indicates the direction
        // of one side (the one used to calculate azimuth angle phi)
        // v2 is a direction so it should be normalized. the module
        // automatically normalizes it
        double ris_v2_x = default(0);
        double ris_v2_y = default(0);
        double ris_v2_z = default(0);
        // vector v1 represents the other side of the RIS rectangle
        // notice that the orientation of the v1 and v2 vectors determine
        // the direction the RIS is facing (its reflective side).
        // for example, if we define v1 = [1 0 0] and v2 = [0 0 -1]
        // we determine the direction the RIS is facing by computing
        // the cross product between v1 and v2
        // so we obtain v1 x v2 = [0 1 0], an RIS parallel to the x-z
        // plane facing towards positive y values. intuitively, this
        // can be computed using the right-hand rule
        double ris_v1_x = default(0);
        double ris_v1_y = default(0);
        double ris_v1_z = default(0);

        // when to apply the initial angles configuration of the metasurface
        double initialConfigurationTime @unit(s) = default(1us);
        // full path name of node to initially configure incidence/reflection angles
        string focusBeamFrom = default("");
        string pointBeamTo = default("");
        // configure RIS on startup
        // incidence angles are ignored if focusBeamFrom is set
        double initialIncidencePhi @unit(rad) = default(0rad);
        double initialIncidenceTheta @unit(rad) = default(0rad);
        // reflection angles are ignored if pointBeamTo is set
        double initialReflectionPhi @unit(rad) = default(0rad);
        double initialReflectionTheta @unit(rad) = default(0rad);
        // same as the previous parameters, but enables to specify a list of phis and thetas, which will then be combined into a single configuration
        // if this is not set as an empty string (""), this takes precedence over the previous two parameters
        // as for the previous two parameters, this is ignored if pointBeamTo is set
        string initialReflectionsPhi_rad = default("");
        string initialReflectionsTheta_rad = default("");

        // determine whether to use the random strategy when combining multiple configurations or the averaging one
        bool randomCombinationStrategy = default(true);

        // whenever a frame is reflected, reconfigure the RIS to point towards that node
        // this can be set as a space-separated list of nodes to track multiple ones by combining multiple configurations 
        string destinationNodesToTrack = default("");

        // whether to ignore non reflected frame on reception (capture effect, as non reflected signal comes in first)
        bool ignoreNonReflectedSignals = default(false);
        // whether to ignore signals that are shadowed (i.e., crossing through objects)
        bool ignoreShadowedSignals = default(false);
        // whether to receive every packet regardless of noise and interference
        bool ignoreNoiseAndInterference = default(false);

        // whether to use sum or product of distances (near field vs far field) in the free space path loss model
        bool useProductOfDistances = default(false);

        // power below which a RIS will not reflect a signal (for optimization)
        double repropagationThreshold @unit(mW) = default(-100dBm);

        // number of coding states (discretization of phase, n parameter)
        int codingStates = default(4);
        // number of unit cells per lambda
        int cellsPerLambda = default(3);
        // size of a side of the metasurface in units of lambda
        // the total number of elements (assuming a square surfaces) is (lambdaSize * cellsPerLambda)^2
        int lambdaSize = default(5);
        // angular resolution of the model for phi and theta
        double d_phi @unit(rad) = default(1deg);
        double d_theta @unit(rad) = default(1deg);

}
