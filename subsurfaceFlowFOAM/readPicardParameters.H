/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
/*-----------------------------------------------------------------------*\
Class
    Foam::readPicardParameters[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Reads the input parameters required for Picard Iteration Method.
\*-----------------------------------------------------------------------*/
const label NPl = picardDict.lookupOrDefault<label>("Picard_lower_limit", 1);
const label NPu = picardDict.lookupOrDefault<label>("Picard_upper_limit", 1);
const label NPmax = picardDict.lookupOrDefault<label>("Max_Picard_iterations", 1);
const label Nstab_max = picardDict.lookupOrDefault<label>("Max_Stabilization_Counters", 1);
const scalar delta_f = picardDict.lookupOrDefault<scalar>("Time_factor", 1.0);
const scalar minDeltaT = picardDict.lookupOrDefault<scalar>("minDeltaT", 1.0);
