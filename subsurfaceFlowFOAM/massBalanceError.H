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
    Foam::massBalanceError[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Calculates the cumulative Mass balance Error after each time-step.
\*-----------------------------------------------------------------------*/
scalar MB1 = 0.0;
scalar MB2_temp = 0.0;
scalar MBE = 0.0;
		
if (nMBE != 1)
{
	//Calculation of Total additional mass in the domain (MB1)
	MB1 = gSum((((theta - theta0) + (Ss*theta/eta)*(h - h0))*mesh.V())());
		
	//Calculation of Total net flux into the domain (MB2)
	forAll (mesh.boundaryMesh(), i)
	{
		forAll(phi.boundaryField()[i], j)
		{
			MB2_temp = MB2_temp + (-phi.boundaryField()[i][j]);
		}
	}
		
	reduce(MB2_temp, sumOp<scalar>());
			
	MB2 = MB2 + MB2_temp*runTime.deltaT().value();
			
	MBE = mag(1.0 - (MB1/MB2))*100;
	
	Info << "\nCummulative Mass Balance Error = " << MBE << " % " << endl;

}
