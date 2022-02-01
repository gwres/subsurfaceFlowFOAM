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
    Foam::solveTan

Group
    wellUnconfinedBC/customBoundaryConditions

Description
    Solves for x: x*tan(x) = a
\*-----------------------------------------------------------------------*/
#ifndef solveTan_H
#define solveTan_H

namespace Foam
{
	inline scalar solTan(scalar x0, scalar a)
	{
		scalar x = x0;
		scalar f = 0.0;
		scalar df = 0.0;
		scalar RMSE = 1.0;
		while (RMSE > 1e-9)
		{
			RMSE = 0.0;
			f = x*sin(x) - a*cos(x);
			df = x*cos(x) + a*sin(x) + sin(x);
			x = x - (f/df);
			RMSE = mag(f/df);
		}
		return x;
	}
}//End namespace Foam

#endif
