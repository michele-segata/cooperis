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

// Based on Car.ned from Christoph Sommer


package org.car2x.veinsris;

import org.car2x.veins.base.modules.*;
import org.car2x.veinsris.NicRis;
import org.car2x.veinsris.Splitter;

module CarRis
{
    parameters:
        string applType; //type of the application layer
        string nicType = default("Nic80211p"); // type of network interface card
        string veinsmobilityType = default("org.car2x.veins.modules.mobility.traci.TraCIMobility"); //type of the mobility module
        @display("bgb=457,459");
    gates:
        input veinsradioIn; // gate for sendDirect
        input risIn;
    submodules:
        application: <applType> like org.car2x.veins.base.modules.IBaseApplLayer {
            parameters:
                @display("p=163,50");
        }

        nic: <nicType> like org.car2x.veins.modules.nic.INic80211p {
            parameters:
                @display("p=71,243");
        }

        veinsmobility: <veinsmobilityType> like org.car2x.veins.base.modules.IMobility {
            parameters:
                @display("p=368,127;i=block/cogwheel");
        }

        nicRis: NicRis {
            parameters:
                isRis = false;
                @display("p=163,243");
        }

        splitter: Splitter {
            @display("p=163,127");
        }

    connections allowunconnected:
        application.lowerLayerOut --> splitter.applicationIn;
        application.lowerLayerIn <-- splitter.applicationOut;

        splitter.nicOut --> nic.upperLayerIn;
        splitter.nicIn <-- nic.upperLayerOut;

        splitter.nicRisOut --> nicRis.upperLayerIn;
        splitter.nicRisIn <-- nicRis.upperLayerOut;

        veinsradioIn --> nic.radioIn;
        risIn --> nicRis.radioIn;
}
