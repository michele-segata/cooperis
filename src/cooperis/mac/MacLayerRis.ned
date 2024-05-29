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

package org.car2x.cooperis.mac;

import org.car2x.veins.base.modules.BaseMacLayer;

simple MacLayerRis extends BaseMacLayer
{
    parameters:
        @class(veins::MacLayerRis);

        // assigning default value to the MAC header
        headerLength @unit(bit) = default(0 bit);

        //the maximum queue size of MAC queue. 0 for unlimited. Queue strategy is "drop if full"
        int queueSize = default(0);

}
