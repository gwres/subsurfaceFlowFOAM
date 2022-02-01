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
    Foam::besselK0

Group
    wellUnconfinedBC/customBoundaryConditions

Description
    Modified Bessel Function of Second Kind and Order Zero 
\*-----------------------------------------------------------------------*/
#ifndef besselK0_H
#define besselK0_H

#include "besselI0.H"

namespace Foam
{
	inline scalar besselK0(scalar X)
	{
		constexpr scalar A1 = -0.57721566;
		constexpr scalar A2 = 0.42278420;
		constexpr scalar A3 = 0.23069756;
		constexpr scalar A4 = 0.03488590;
		constexpr scalar A5 = 0.00262698;
        constexpr scalar A6 = 0.00010750;
		constexpr scalar A7 = 0.00000740;
		
		constexpr scalar B1 = 1.25331414;
		constexpr scalar B2 = -0.07832358;
		constexpr scalar B3 = 0.02189568;
		constexpr scalar B4 = -0.01062446;
		constexpr scalar B5 = 0.00587872;
		constexpr scalar B6 = -0.00251540;
		constexpr scalar B7 = 0.00053208;
		
		scalar K;
		if (X <= 2.0)
        {
			const scalar T = sqr(X/2.0);
			K = -(log(X/2.0)*besselI0(X)) + A1 + T*(A2 + T*(A3 + T*(A4 + T*(A5 + T*(A6 + T*A7)))));
		}
		else
		{
			const scalar T = 2.0/X;
			K = B1 + T*(B2 + T*(B3 + T*(B4 + T*(B5 + T*(B6 + T*B7)))));
			K /= (exp(X)*sqrt(X));
		}
    return K;
	}
}// End namespace Foam

#endif
