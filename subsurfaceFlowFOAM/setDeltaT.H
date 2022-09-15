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
    Foam::setDeltaT[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Adjusts the time-step according the number of Picard Iterations.
\*-----------------------------------------------------------------------*/
if (adjustTimeStep)
{	
	//Maximum number of Picard Iteration is exceeded and time step needs to be decreased 
	if ((NP >= NPu) || (NP == NPmax))
	{
		runTime.setDeltaT
		(
			max
			(
				(1.0/delta_f)*runTime.deltaT().value(),
				minDeltaT
			) 
		);
		sc = 0;
	}
	
	//Number of Picard iterations is greater than the Minimum value but less than Maximum value
	if ((NP >= NPl) && (NP < NPu))
	{
		sc = 0;
	}

    //Number of Picard iterations is less than the Minimum value but maximum number of Stabilisation Counters not reached
    if ((NP < NPl) && (sc < Nstab_max)) 
    {
		sc = sc+1;
	}

	//Number of Picard iterations is less than the Minimum value but maximum number of Stabilisation Counters reached
	//hence time-step needs to be increased
	if ((NP < NPl) && (sc == Nstab_max)) 
    {
		runTime.setDeltaT
		(
			min
			(
				delta_f*runTime.deltaT().value(), 
				maxDeltaT
			)
		);
		sc = 0; 
    }
}
