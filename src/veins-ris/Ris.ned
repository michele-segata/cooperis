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
import org.car2x.veinsris.SUMOBaseMobility;
import org.car2x.veinsris.NicRis;

module Ris
{
    parameters:
        string mobility_type = default("SUMOBaseMobility");
        @display("bgb=457,459");
    gates:
        input risIn;
    submodules:

        mobility: <mobility_type> like IMobility {
            parameters:
                @display("p=368,127;i=block/cogwheel");
        }
        
        nicRis: NicRis {
            parameters:
                isRis = true;
                @display("p=163,243");
        }

    connections allowunconnected:
        risIn --> nicRis.radioIn;
}
