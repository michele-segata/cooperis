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

package org.car2x.veinsris;

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
        // for example, if we define v2 = [-1 0 0] and v1 = [0 0 1]
        // we determine the direction the RIS is facing by computing
        // the cross product between v2 and v1 (v2 x v1, AND NOT v1 x v2!)
        // so we obtain v2 x v1 = [0 1 0], an RIS parallel to the x-z
        // plane facing towards positive y values. intuitively, this
        // can be computed using the right-hand rule (or better left-hand,
        // as OMNeT is left-handed)
        double ris_v1_x = default(0);
        double ris_v1_y = default(0);
        double ris_v1_z = default(0);
        
        // whether to ignore non reflected frame on reception (capture effect, as non reflected signal comes in first)
        bool ignoreNonReflectedSignals = default(false);
        // whether to receive every packet regardless of noise and interference
        bool ignoreNoiseAndInterference = default(false);

}
