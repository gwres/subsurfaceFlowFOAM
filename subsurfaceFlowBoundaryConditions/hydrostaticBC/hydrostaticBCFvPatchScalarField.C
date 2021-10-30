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
    Foam::hydrostaticBCFvPatchScalarField
    
Group
	hydrostaticBC/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the hydrostatic pressure head to 
    be specified on a boundary.
\*---------------------------------------------------------------------------*/
#include "hydrostaticBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
	//Constructors
	
	hydrostaticBCFvPatchScalarField::hydrostaticBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		zWT_(p.size(), 0.0)
		{}
	/*--------------------------------------------------------------------------*/
	hydrostaticBCFvPatchScalarField::hydrostaticBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		zWT_("zWT", dict, p.size())
		{
			if (dict.found("value"))
			{
				fvPatchField<scalar>::operator = (scalarField("value", dict, p.size()));
				fixedValueFvPatchScalarField::updateCoeffs();
				fixedValueFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
			}
		}
	/*--------------------------------------------------------------------------*/
	hydrostaticBCFvPatchScalarField::hydrostaticBCFvPatchScalarField
	(
		const hydrostaticBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedValueFvPatchScalarField(ptf, p, iF, mapper),
		zWT_(ptf.zWT_, mapper)
		{}
	/*--------------------------------------------------------------------------*/
	hydrostaticBCFvPatchScalarField::hydrostaticBCFvPatchScalarField
	(
		const hydrostaticBCFvPatchScalarField& hpsf
	)
	:
		fixedValueFvPatchScalarField(hpsf),
		zWT_(hpsf.zWT_)
		{}
	/*--------------------------------------------------------------------------*/
	hydrostaticBCFvPatchScalarField::hydrostaticBCFvPatchScalarField
	(
		const hydrostaticBCFvPatchScalarField& hpsf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(hpsf, iF),
		zWT_(hpsf.zWT_)
		{}

	//Member Functions

	void hydrostaticBCFvPatchScalarField::autoMap
	(
		const fvPatchFieldMapper& m
	)
	{
		fixedValueFvPatchScalarField::autoMap(m);
		zWT_.autoMap(m);
	}

	void hydrostaticBCFvPatchScalarField::rmap
	(
		const fvPatchScalarField& ptf,
		const labelList& addr
	)
	{
		fixedValueFvPatchScalarField::rmap(ptf, addr);

		const hydrostaticBCFvPatchScalarField& hptf = refCast<const hydrostaticBCFvPatchScalarField>(ptf);
		zWT_.rmap(hptf.zWT_, addr);
	}

	void hydrostaticBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}

		const fvPatchField<scalar>& z = patch().lookupPatchField<volScalarField, scalar>("z");
		
		operator == (zWT_ - z);

		fixedValueFvPatchScalarField::updateCoeffs();
	}

	void hydrostaticBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		zWT_.writeEntry("zWT", os);
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchScalarField, hydrostaticBCFvPatchScalarField);

} //End namespace Foam

