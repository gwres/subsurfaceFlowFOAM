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
    Foam::besselI1

Group
    wellUnconfinedBC/subsurfaceFlowBoundaryConditions

Description
    Modified Bessel Function of First Kind and Order One 
\*-----------------------------------------------------------------------*/
#ifndef besselI1_H
#define besselI1_H

namespace Foam
{
	inline scalar besselI1(scalar X)
	{
		constexpr scalar A1 = 0.5;
		constexpr scalar A2 = 0.87890594;
		constexpr scalar A3 = 0.51498869;
		constexpr scalar A4 = 0.15084934;
		constexpr scalar A5 = 0.02658733;
        constexpr scalar A6 = 0.00301532;
		constexpr scalar A7 = 0.00032411;
		
		constexpr scalar B1 = 0.39894228;
		constexpr scalar B2 = -0.03988024;
		constexpr scalar B3 = -0.00362018;
		constexpr scalar B4 = 0.00163801;
		constexpr scalar B5 = -0.01031555;
		constexpr scalar B6 = 0.02282967;
		constexpr scalar B7 = -0.02895312;
		constexpr scalar B8 = 0.01787654;
		constexpr scalar B9 = -0.00420059;
		
		scalar I;
		if (mag(X) < 3.75)
        {
			const scalar T = sqr(X/3.75);
			I = X*(A1 + T*(A2 + T*(A3 + T*(A4 + T*(A5 + T*(A6 + T*A7))))));
		}
		else
		{
			const scalar T = 3.75/mag(X);
			I = B1 + T*(B2 + T*(B3 + T*(B4 + T*(B5 + T*(B6 + T*(B7 + T*(B8 + T*B9)))))));
			I *= (exp(mag(X))/sqrt(mag(X)));
		}
    return I;
	}
}//End namespace Foam

#endif
