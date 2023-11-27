#!/usr/bin/env python3

#
# Copyright (C) 2013-2019 Christoph Sommer <sommer@ccs-labs.org>
#
# Documentation for these modules is at http://veins.car2x.org/
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"""
Wrapper around cooperis_catch binary
"""

import os
import argparse
import subprocess

# v-- contents of out/config.py go here
# ^-- contents of out/config.py go here

def relpath(s):
    cooperis_root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
    return os.path.relpath(os.path.join(cooperis_root, s), '.')

parser = argparse.ArgumentParser('Run a Cooperis simulation')
parser.add_argument('-d', '--debug', action='store_true', help='Set --mode=debug (deprecated in favor of --mode)')
parser.add_argument('-M', '--mode', metavar='MODE', dest='mode', choices=['', 'release', 'debug', 'sanitize'], help='Instead of opp_run, use opp_run_VARIANT corresponding to MODE (release, debug, sanitize)')
parser.add_argument('-t', '--tool', metavar='TOOL', dest='tool', choices=['lldb', 'gdb', 'memcheck', 'callgrind'], help='Wrap opp_run execution in TOOL (lldb, gdb, memcheck, or callgrind)')
parser.add_argument('-v', '--verbose', action='store_true', help='Print command line before executing')
parser.add_argument('--', dest='arguments', help='Arguments to pass to opp_run')
args, omnet_args = parser.parse_known_args()
if (len(omnet_args) > 0) and omnet_args[0] == '--':
    omnet_args = omnet_args[1:]

opp_run = 'cooperis_catch'
if args.debug:
    args.mode = 'debug'
if args.mode:
    if args.mode == '':
        opp_run = 'cooperis_catch'
    elif args.mode == 'release':
        opp_run = 'cooperis_catch'
    elif args.mode == 'debug':
        opp_run = 'cooperis_catch_dbg'
    elif args.mode == 'sanitize':
        opp_run = 'cooperis_catch_sanitize'
    else:
        assert False, 'unknown --mode option'

prefix = []
if args.tool == 'lldb':
    prefix = ['lldb', '--']
if args.tool == 'gdb':
    prefix = ['gdb', '--args']
if args.tool == 'memcheck':
    prefix = ['valgrind', '--tool=memcheck', '--leak-check=full', '--dsymutil=yes', '--log-file=valgrind.out']
if args.tool == 'callgrind':
    prefix = ['valgrind', '--tool=callgrind', '--log-file=callgrind.out']

cmdline = prefix + [os.path.join(relpath('src'), opp_run)] + omnet_args

if args.verbose:
    print("Running with command line arguments: %s" % ' '.join(['"%s"' % arg for arg in cmdline]))

if os.name == 'nt':
    subprocess.call(['env'] + cmdline)
else:
    os.execvp('env', ['env'] + cmdline)

