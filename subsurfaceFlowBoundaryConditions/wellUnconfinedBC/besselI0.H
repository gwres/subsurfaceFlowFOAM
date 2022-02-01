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
    Foam::besselI0

Group
    wellUnconfinedBC/customBoundaryConditions

Description
    Modified Bessel Function of First Kind and Order Zero 
\*-----------------------------------------------------------------------*/
#ifndef besselI0_H
#define besselI0_H

namespace Foam
{
	inline scalar besselI0(scalar X)
	{
		constexpr scalar A1 = 1.0;
		constexpr scalar A2 = 3.5156229;
		constexpr scalar A3 = 3.0899424;
		constexpr scalar A4 = 1.2067492;
		constexpr scalar A5 = 0.2659732;
		constexpr scalar A6 = 0.0360768;
		constexpr scalar A7 = 0.0045813;
		
		constexpr scalar B1 = 0.39894228;
		constexpr scalar B2 = 0.01328592;
		constexpr scalar B3 = 0.00225319;
		constexpr scalar B4 = -0.00157565;
		constexpr scalar B5 = 0.00916281;
		constexpr scalar B6 = -0.02057706;
		constexpr scalar B7 = 0.02635537;
		constexpr scalar B8 = -0.01647633;
		constexpr scalar B9 = 0.00392377;
		
		scalar I;
		if (mag(X) < 3.75)
		{
			const scalar T = sqr(X/3.75);
			I = A1 + T*(A2 + T*(A3 + T*(A4 + T*(A5 + T*(A6 + T*A7)))));
		}
		else
		{
			const scalar T = 3.75/mag(X);
			I = B1 + T*(B2 + T*(B3 + T*(B4 + T*(B5 + T*(B6 + T*(B7 + T*(B8 + T*B9)))))));
			I *= (exp(mag(X))/sqrt(mag(X)));
		}
		return I;
	}
}// end namespace Foam

#endif

