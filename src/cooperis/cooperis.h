//
// Copyright (C) 2024 Michele Segata <segata@ccs-labs.org>
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

#pragma once

#include "veins/veins.h"

// Version number of last release ("major.minor.patch") or an alpha version, if nonzero
#define COOPERIS_VERSION_MAJOR 1
#define COOPERIS_VERSION_MINOR 0
#define COOPERIS_VERSION_PATCH 0
#define COOPERIS_VERSION_ALPHA 0

// Explicitly check Veins version number
#if !(VEINS_VERSION_MAJOR == 5 && VEINS_VERSION_MINOR >= 0)
#error Veins version 5.0 or compatible required
#endif

// COOPERIS_API macro. Allows us to use the same .h files for both building a .dll and linking against it
#if defined(COOPERIS_EXPORT)
#define COOPERIS_API OPP_DLLEXPORT
#elif defined(COOPERIS_IMPORT)
#define COOPERIS_API OPP_DLLIMPORT
#else
#define COOPERIS_API
#endif
