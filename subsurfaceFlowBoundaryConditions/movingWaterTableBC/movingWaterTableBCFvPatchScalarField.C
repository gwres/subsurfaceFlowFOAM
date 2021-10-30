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
    Foam::movingWaterTableBCFvPatchScalarField
    
Group
	movingWaterTableBC/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the pressure head to be specified
    at the bottom boundary for a moving water table.
\*---------------------------------------------------------------------------*/
#include "movingWaterTableBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
	//Constructors
	
	movingWaterTableBCFvPatchScalarField::movingWaterTableBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		zWT_(p.size(), 0.0),
		vWT_(p.size(), 0.0)
		{}
	/*--------------------------------------------------------------------------*/
	movingWaterTableBCFvPatchScalarField::movingWaterTableBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		zWT_("zWT", dict, p.size()),
		vWT_("vWT", dict, p.size())
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
	movingWaterTableBCFvPatchScalarField::movingWaterTableBCFvPatchScalarField
	(
		const movingWaterTableBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedValueFvPatchScalarField(ptf, p, iF, mapper),
		zWT_(ptf.zWT_, mapper),
		vWT_(ptf.vWT_, mapper)
		{}
	/*--------------------------------------------------------------------------*/
	movingWaterTableBCFvPatchScalarField::movingWaterTableBCFvPatchScalarField
	(
		const movingWaterTableBCFvPatchScalarField& mwtpsf
	)
	:
		fixedValueFvPatchScalarField(mwtpsf),
		zWT_(mwtpsf.zWT_),
		vWT_(mwtpsf.vWT_)
		{}
	/*--------------------------------------------------------------------------*/
	movingWaterTableBCFvPatchScalarField::movingWaterTableBCFvPatchScalarField
	(
		const movingWaterTableBCFvPatchScalarField& mwtpsf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(mwtpsf, iF),
		zWT_(mwtpsf.zWT_),
		vWT_(mwtpsf.vWT_)
		{}

	//Member Functions

	void movingWaterTableBCFvPatchScalarField::autoMap
	(
		const fvPatchFieldMapper& m
	)
	{
		fixedValueFvPatchScalarField::autoMap(m);
		zWT_.autoMap(m);
		vWT_.autoMap(m);
	}

	void movingWaterTableBCFvPatchScalarField::rmap
	(
		const fvPatchScalarField& ptf,
		const labelList& addr
	)
	{
		fixedValueFvPatchScalarField::rmap(ptf, addr);

		const movingWaterTableBCFvPatchScalarField& mwtptf = refCast<const movingWaterTableBCFvPatchScalarField>(ptf);
		zWT_.rmap(mwtptf.zWT_, addr);
		vWT_.rmap(mwtptf.vWT_, addr);
	}

	void movingWaterTableBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}

		scalar t = this->db().time().value();
		
		operator== (zWT_ + vWT_*t);

		fixedValueFvPatchScalarField::updateCoeffs();
	}

	void movingWaterTableBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		zWT_.writeEntry("zWT", os);
		vWT_.writeEntry("vWT", os);
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchScalarField, movingWaterTableBCFvPatchScalarField);

} //End namespace Foam

