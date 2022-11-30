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

package org.car2x.veinsris;

import org.car2x.veinsris.mac.MacLayerRis;
import org.car2x.veinsris.PhyLayerRis;

module NicRis
{
	parameters:
	    
	    bool isRis;
	    // Explicitly specify the Connection Manager for this type of NIC
        // If not specified, MiXiM defaults to the BaseConnectionManager
        string connectionManagerName = default("risConnectionManager");
	    
    gates:
        input radioIn; // radioIn gate for sendDirect
        
        input upperLayerIn;
        output upperLayerOut;

        output upperControlOut;
        input upperControlIn;

    submodules:

        macRis: MacLayerRis if !isRis {
            @display("p=69,82");
        }
		
		// the reflective metasurface only has a phy layer.
        // in fact, there is not even a real phy layer, as it does only reflect a signal and it adds metadata to it
        phyRis: PhyLayerRis {
            parameters:
                isRIS = parent.isRis;
            @display("p=69,218;i=block/process_s");
        }

    connections allowunconnected:
        radioIn --> phyRis.radioIn;
        
        // Bottom up
        phyRis.upperLayerOut --> macRis.lowerLayerIn if !isRis;
        macRis.upperLayerOut --> upperLayerOut if !isRis;
        // Top down
        upperLayerIn --> macRis.upperLayerIn if !isRis;
        macRis.lowerLayerOut --> phyRis.upperLayerIn if !isRis;
        
        // Bottom up Control
        phyRis.upperControlOut --> macRis.lowerControlIn if !isRis;
        macRis.upperControlOut --> upperControlOut if !isRis;
        // Top down Control
        upperControlIn --> macRis.upperControlIn if !isRis;
        macRis.lowerControlOut --> phyRis.upperControlIn if !isRis;
        
}
