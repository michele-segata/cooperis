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

// Based on AirFrame.msg

import veins.base.messages.AirFrame;
import veins.base.utils.Coord;

namespace veins;

message AirFrameRis extends veins::AirFrame {
    // coordinates describing the orientation of the source RIS (if it is a RIS)
    Coord ris_v1;
    Coord ris_v2;
    Coord ris_vn;
    const cModule* originalSenderModule;
    // list of gains applied by RISs on the path
    double risGains[];
    // total distance travelled, to compute path loss at the end of the path
    double paths[];
    // loss on each path (free space only)
    double loss_dB[];
    double totalDistance;
    // path loss accumulated so far.
    // this can be used to compute the additional pathloss on the new leg of the path,
    // when using the sum of paths
    double currentPathLoss_dB;

    // should the frame be coming reflected by a RIS, the frame will need to
    // include the incidence angle (original transmitter -> RIS)
    double incidencePhi;
    double incidenceTheta;

    double reflectionPhi;
    double reflectionTheta;

    bool reflected;
    int originalId;
}
