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
    Foam::specifiedTransientFluxBC1FvPatchScalarField

Group
    specifiedTransientFluxBC1/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the normal gradient of pressure head 
    due to transient specified flux value at that boundary as in Paniconi et al.(1991).
\*-----------------------------------------------------------------------*/
#include "specifiedTransientFluxBC1FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	specifiedTransientFluxBC1FvPatchScalarField::specifiedTransientFluxBC1FvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_(0.0)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC1FvPatchScalarField::specifiedTransientFluxBC1FvPatchScalarField
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
	specifiedTransientFluxBC1FvPatchScalarField::specifiedTransientFluxBC1FvPatchScalarField
	(
		const specifiedTransientFluxBC1FvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		qC_(ptf.qC_, mapper)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC1FvPatchScalarField::specifiedTransientFluxBC1FvPatchScalarField
	(
		const specifiedTransientFluxBC1FvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		qC_(ptf.qC_)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedTransientFluxBC1FvPatchScalarField::specifiedTransientFluxBC1FvPatchScalarField
	(
		const specifiedTransientFluxBC1FvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		qC_(ptf.qC_)
		{}
	
	//Member Functions
	
	void specifiedTransientFluxBC1FvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		scalar t = this->db().time().value();
		
		const fvPatchField<tensor>& K_C = patch().lookupPatchField<volTensorField, tensor>("K");
		
        const fvPatchField<vector>& grad_z_C = patch().lookupPatchField<volVectorField, vector>("grad_z");
        
        gradient() = -(((qC_*(t/8.2944e8)) & inv(K_C)) & patch().nf()) - (grad_z_C & patch().nf());
                
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void specifiedTransientFluxBC1FvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		qC_.writeEntry("qC", os);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, specifiedTransientFluxBC1FvPatchScalarField);
			
} //End namespace Foam
