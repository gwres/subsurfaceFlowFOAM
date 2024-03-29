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
    Foam::convergenceCriteria[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Calculates the Convergence Criterion for Picard Iteration Method.
\*-----------------------------------------------------------------------*/

if (convergenceCriterion.match("h-based"))
{
	errorSum = gSumMag(Foam::pow((h - hM), 2.0)())/N_ele;
	error = Foam::pow(errorSum, 0.5);
}
else if (convergenceCriterion.match("mixed_h-based"))
{
	scalar sum1 = gSumMag(Foam::pow((h - hM), 2.0)())/N_ele;
	scalar sum2 = gSumMag(Foam::pow(h, 2.0)())/N_ele;
	error = Foam::pow(sum1, 0.5) - delta_hretol*Foam::pow(sum2, 0.5);
}
else if (convergenceCriterion.match("theta-based"))
{
	errorSum = gSumMag(Foam::pow((theta - thetaM), 2.0)())/N_ele;
	error = Foam::pow(errorSum, 0.5);
}
