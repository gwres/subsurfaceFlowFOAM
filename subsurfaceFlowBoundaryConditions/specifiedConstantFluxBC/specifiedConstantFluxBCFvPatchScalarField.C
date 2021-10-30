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
    Foam::specifiedConstantFluxBCFvPatchScalarField

Group
    specifiedConstantFluxBC/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the surface normal gradient for 
    pressure head due to time-constant flux specified at that boundary.
\*-----------------------------------------------------------------------*/
#include "specifiedConstantFluxBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	specifiedConstantFluxBCFvPatchScalarField::specifiedConstantFluxBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_(p.size(), 0.0)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantFluxBCFvPatchScalarField::specifiedConstantFluxBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_("qC", dict, p.size())
		{
			if (dict.found("gradient"))
			{
				gradient() = scalarField("gradient", dict, p.size());
				fixedGradientFvPatchScalarField::updateCoeffs();
				fixedGradientFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
				gradient() = 0.0;
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	specifiedConstantFluxBCFvPatchScalarField::specifiedConstantFluxBCFvPatchScalarField
	(
		const specifiedConstantFluxBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		qC_(ptf.qC_, mapper)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantFluxBCFvPatchScalarField::specifiedConstantFluxBCFvPatchScalarField
	(
		const specifiedConstantFluxBCFvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		qC_(ptf.qC_)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantFluxBCFvPatchScalarField::specifiedConstantFluxBCFvPatchScalarField
	(
		const specifiedConstantFluxBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		qC_(ptf.qC_)
		{}
	
	//Member Functions
	
	void specifiedConstantFluxBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const fvPatchField<tensor>& K_C = patch().lookupPatchField<volTensorField, tensor>("K");
		const scalarField K_nf = (K_C & cmptMag(patch().nf())) & cmptMag(patch().nf());
        
        const fvPatchField<vector>& grad_z_C = patch().lookupPatchField<volVectorField, vector>("grad_z");
                
        gradient() = qC_/-K_nf - (grad_z_C & patch().nf());
        
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void specifiedConstantFluxBCFvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		qC_.writeEntry("qC", os);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, specifiedConstantFluxBCFvPatchScalarField);
			
} //End namespace Foam
