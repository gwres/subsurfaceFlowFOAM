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
/*---------------------------------------------------------------------------*\
Class
    Foam::stehfestCoefficients

Group
    wellUnconfinedBC/subsurfaceFlowBoundaryConditions

Description
    This calculates the coefficients required in Stehfest Algorithm for 
    the Inverse Laplacian transformation of a variable.
\*-----------------------------------------------------------------------*/
#ifndef stehfestCoefficients_H
#define stehfestCoefficients_H

#include "factorial.H"

namespace Foam
{
	inline scalar stehfestCoef(label n, label Ns)
	{
		const label k1 = (n + 1)/2;
		const label k2 = std::min(n, Ns/2);
		scalar wt = 0;
		label ct = k1;	
		for (label j = k1; j <= k2; j++)
		{
			scalar A1 = fact(2*ct);
			scalar A2 = fact((Ns/2) - ct);
			scalar A3 = fact(ct);
			scalar A4 = fact(ct - 1);
			scalar A5 = fact(n - ct);
			scalar A6 = fact(2*ct - n);
			wt += (A1*std::pow(ct, Ns/2))/(A2*A3*A4*A5*A6);
			ct += 1;
		}
		wt *= std::pow(-1, (n + (Ns/2)));
		return wt;
	}
} //End namespace Foam

#endif
