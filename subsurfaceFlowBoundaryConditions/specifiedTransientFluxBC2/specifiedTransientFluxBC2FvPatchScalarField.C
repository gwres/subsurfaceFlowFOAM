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
    Foam::specifiedTransientFluxBC2FvPatchScalarField

Group
    specifiedTransientFluxBC2/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the normal gradient of pressure head 
    due to time varying specified flux value at that boundary. The flux is 
    a Function1 type entry ie. it can be read either as a table entry or 
    polynomial function or from a CSV file. 
\*-----------------------------------------------------------------------*/
#include "specifiedTransientFluxBC2FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "one.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	specifiedTransientFluxBC2FvPatchScalarField::specifiedTransientFluxBC2FvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_()
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC2FvPatchScalarField::specifiedTransientFluxBC2FvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF)
		{
			if (dict.found("qC"))
			{
                qC_ = Function1<vector>::New("qC", dict);
			}
			else
			{
				FatalIOErrorInFunction(dict)
				<< "Please supply qC" << nl << exit(FatalIOError);
			}
			
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
	specifiedTransientFluxBC2FvPatchScalarField::specifiedTransientFluxBC2FvPatchScalarField
	(
		const specifiedTransientFluxBC2FvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		qC_(ptf.qC_.clone())
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC2FvPatchScalarField::specifiedTransientFluxBC2FvPatchScalarField
	(
		const specifiedTransientFluxBC2FvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		qC_(ptf.qC_.clone())
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC2FvPatchScalarField::specifiedTransientFluxBC2FvPatchScalarField
	(
		const specifiedTransientFluxBC2FvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		qC_(ptf.qC_.clone())
		{}
	
	//Member Functions
	
	void specifiedTransientFluxBC2FvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const scalar t = db().time().timeOutputValue();
		
		const fvPatchField<tensor>& K_C = patch().lookupPatchField<volTensorField, tensor>("K");
		
        const fvPatchField<vector>& grad_z_C = patch().lookupPatchField<volVectorField, vector>("grad_z");
        
        gradient() = -((qC_->value(t) & inv(K_C)) & patch().nf()) - (grad_z_C & patch().nf());
                
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void specifiedTransientFluxBC2FvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		if (qC_.valid())
		{
			qC_->writeData(os);
		}
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, specifiedTransientFluxBC2FvPatchScalarField);
			
} //End namespace Foam
