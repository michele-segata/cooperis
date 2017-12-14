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


package org.car2x.veinsrms;

import org.car2x.veinsrms.mac.MacLayerRms;
import org.car2x.veinsrms.PhyLayerRms;

module NicRms
{
	parameters:
	    
	    bool isRms;
	    // Explicitly specify the Connection Manager for this type of NIC
        // If not specified, MiXiM defaults to the BaseConnectionManager
        string connectionManagerName = default("rmsConnectionManager");
	    
    gates:
        input radioIn; // radioIn gate for sendDirect
        
        input upperLayerIn;
        output upperLayerOut;

        output upperControlOut;
        input upperControlIn;

    submodules:

        macRms: MacLayerRms if !isRms {
            @display("p=69,82");
        }
		
		// the reflective metasurface only has a phy layer.
        // in fact, there is not even a real phy layer, as it does only reflect a signal and it adds metadata to it
        phyRms: PhyLayerRms {
            parameters:
                isRMS = parent.isRms;
            @display("p=69,218;i=block/process_s");
        }

    connections allowunconnected:
        radioIn --> phyRms.radioIn;
        
        // Bottom up
        phyRms.upperLayerOut --> macRms.lowerLayerIn if !isRms;
        macRms.upperLayerOut --> upperLayerOut if !isRms;
        // Top down
        upperLayerIn --> macRms.upperLayerIn if !isRms;
        macRms.lowerLayerOut --> phyRms.upperLayerIn if !isRms;
        
		// Bottom up Control
		phyRms.upperControlOut --> macRms.lowerControlIn if !isRms;
		macRms.upperControlOut --> upperControlOut if !isRms;
		// Top down Control        
        upperControlIn --> macRms.upperControlIn if !isRms;
        macRms.lowerControlOut --> phyRms.upperControlIn if !isRms;
        
}
