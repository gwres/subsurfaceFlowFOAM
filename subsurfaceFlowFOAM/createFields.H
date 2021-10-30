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
    Foam::createFields[subsurfaceFlowFOAM]

Group
    subsurfaceFlowFOAM

Description
    Modeling of Subsurface Flow [Richards Equation]
/*------------------------------------------------------------------------------*/
Info<< "Reading PRESSURE HEAD field h\n" << endl;

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SOURCE/SINK TERM field qS\n" << endl;

volScalarField qS
(
    IOobject
    (
        "qS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Creating FACE-CENTRED FLUX field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading UNIT PRESSURE HEAD field hU\n" << endl;

volScalarField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading ELEVATION field z\n" << endl;

volScalarField z
(
	IOobject
    (
		"z",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.C().component(2)
);
/*------------------------------------------------------------------------------*/
Info<< "Creating GRADIENT OF ELEVATION field grad_z\n" << endl;

volVectorField grad_z
(
	IOobject
    (
		"grad_z",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(z)
);
/*------------------------------------------------------------------------------*/
// FIELDS FOR SOIL PROPERTIES
/*------------------------------------------------------------------------------*/
Info<< "Reading SATURATED HYDRAULIC CONDUCTIVITY field Ks\n" << endl;

volTensorField Ks
(
    IOobject
    (
        "Ks",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Creating UNSATURATED HYDRAULIC CONDUCTIVITY field K\n" << endl;

volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Creating FACE-CENTRED UNSATURATED HYDRAULIC CONDUCTIVITY field Kfc\n" << endl;

surfaceTensorField Kfc
(
    IOobject
    (
        "Kfc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SATURATED MOISTURE CONTENT field thetaS\n" << endl;

volScalarField thetaS
(
    IOobject
    (
        "thetaS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading RESIDUAL MOISTURE CONTENT field thetaR\n" << endl;

volScalarField thetaR
(
    IOobject
    (
        "thetaR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Creating MOISTURE CONTENT field theta\n" << endl;

volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SOIL POROSITY field eta\n" << endl;

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SPECIFIC STORAGE field Ss\n" << endl;

volScalarField Ss
(
    IOobject
    (
        "Ss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SOIL-SPECIFIC COEFFICIENT field alphaW\n" << endl;

volScalarField alphaW
(
    IOobject
    (
        "alphaW",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading SOIL-SPECIFIC EXPONENT field betaW\n" << endl;

volScalarField betaW
(
    IOobject
    (
        "betaW",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Creating SPECIFIC MOISTURE CAPACITY field Cw\n" << endl;

volScalarField Cw
(
    IOobject
    (
        "Cw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
// Reads user-defined data from picardDict and richardDict in system directory

IOdictionary picardDict
(
	IOobject
    (
		"picardDict",
		runTime.system(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    )
);

IOdictionary richardDict
(
	IOobject
    (
		"richardDict",
		runTime.system(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    )
);
/*------------------------------------------------------------------------------*/

