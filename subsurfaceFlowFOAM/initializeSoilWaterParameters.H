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
    Foam::initializeSoilWaterParameters[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Initializes the soil-water characteristic parameters - moisture content,
    unsaturated hydraulic conductivity and specific moisture capacity based 
    on van Genuchten-Mualem model.
\*-----------------------------------------------------------------------*/
//Soil-specific Exponent
volScalarField gammaW = 1.0 - (1.0/betaW);
    
//--------------------------------------------------------------------------------------------------------------------------//
//Degree of Saturation
volScalarField Se = 0.5*((1.0 + sign(h)) + (1.0 - sign(h))*(Foam::pow((1.0 + Foam::pow(mag(alphaW*h), betaW)), -gammaW))); 

//--------------------------------------------------------------------------------------------------------------------------//
//Moisture Content
theta = thetaR + (thetaS - thetaR)*Se;
theta.correctBoundaryConditions();

volScalarField theta0 = theta;    										//Initial moisture content
volScalarField thetaM = theta; 											//Previous Picard Iteration level moisture content
volScalarField thetaN = theta; 											//Previous time level moisture content
volScalarField thetaNm1 = theta; 										//Previous to previous time level moisture content

//--------------------------------------------------------------------------------------------------------------------------//
//Check whether any Ss-field value is zero 
forAll(mesh.C(), i)
{
	if (Ss.internalField()[i] < 1e-10)
	{
		nSs = 1;
	}
	break;
}

//--------------------------------------------------------------------------------------------------------------------------//
//Specific Moisture Capacity
if (RicForm.match("h-based"))
{
	if (CwForm.match("derivative")) 									//Analytical Derivative form
	{
		Cw = (mag(alphaW)*gammaW*(thetaS - thetaR)/(1.0 - gammaW))*(Foam::pow(Se, (1.0/gammaW)))
		      *(Foam::pow((1.0 - Foam::pow(Se, (1.0/gammaW))), gammaW)); 
	}
	else if (CwForm.match("scs"))										//Standard Chord Slope Approximation form
	{
		if (nSs == 1)
		{
			Cw = (mag(alphaW)*gammaW*(thetaS - thetaR)/(1.0 - gammaW))*(Foam::pow(Se, (1.0/gammaW)))
			      *(Foam::pow((1.0 - Foam::pow(Se, (1.0/gammaW))), gammaW));
		}
		else
		{
			Cw = 0.5*((1.0 - sign(h))*((theta - thetaN)/((hU*pos0(h - hN)*pos0(hN - h)) + (h - hN))));
		}
	}
}
else if (RicForm.match("mixed"))										//Analytical Derivative form
{
	Cw = (mag(alphaW)*gammaW*(thetaS - thetaR)/(1.0 - gammaW))*(Foam::pow(Se, (1.0/gammaW)))
		  *(Foam::pow((1.0 - Foam::pow(Se, (1.0/gammaW))), gammaW));
}
Cw.correctBoundaryConditions();

//--------------------------------------------------------------------------------------------------------------------------//
//Unsaturated Hydraulic Conductivity
K = Ks*(Foam::pow(Se, 0.5)*Foam::pow((1.0 - Foam::pow((1.0 - Foam::pow(Se, (1.0/gammaW))), gammaW)), 2.0));
K.correctBoundaryConditions();

//---------------------------------------------------------------------------------------------------------------------------//
volScalarField F = Ss*theta/eta + Cw;									//Initializing F

volVectorField gTerm = K & grad_z;										//Initializing the Gravity term
	
Kfc = fvc::interpolate(K, "interpolate_K");								//Initializing face-centred Unsaturated Hydraulic Conductivity
	
volVectorField U = -K & (fvc::grad(h) + grad_z);						//Initializing the cell-centred and face-centred flux
phi = fvc::interpolate(U) & mesh.Sf();
