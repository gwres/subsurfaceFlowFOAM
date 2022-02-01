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
    Foam::factorial

Group
    wellUnconfinedBC/customBoundaryConditions

Description
    This calculates the factorial of an integer value
\*-----------------------------------------------------------------------*/
#ifndef factorial_H
#define factorial_H

namespace Foam
{
	inline scalar fact(label n)
	{
		scalar F = 1.0;
		if (n == 0 || n < 0)
		{
			F = 1.0;
		}
		else
		{
			label ct = 1;
			for (label i = 1; i <= n; i++)
			{
				F *= ct;
				ct += 1;
			}
		}
		return F;
	}
}//End namespace Foam

#endif
