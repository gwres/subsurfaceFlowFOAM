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
    Foam::readRichardParameters[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Reads the input parameters for "RICHARDS EQUATION"
\*-----------------------------------------------------------------------*/
word RicForm(richardDict.lookup("RicForm"));
word CwForm(richardDict.lookup("CwForm"));
word convergenceCriterion(richardDict.lookup("convergenceCriterion"));
const scalar delta_htol = richardDict.lookupOrDefault<scalar>("h_tolerance", 1.0);
const scalar delta_hretol = richardDict.lookupOrDefault<scalar>("h_relative_tolerance", 1.0);
const scalar delta_thetatol = richardDict.lookupOrDefault<scalar>("theta_tolerance", 1.0);

