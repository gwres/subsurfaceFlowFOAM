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
#ifndef movingWaterTableBCFvPatchScalarField_H
#define movingWaterTableBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{
class movingWaterTableBCFvPatchScalarField
:
    public fixedValueFvPatchScalarField
	{
    
		//Private data
		scalarField zWT_;												//Initial height of the water table
		scalarField vWT_;												//Velocity of the moving water table

		public:

			//Runtime type information
			TypeName("movingWaterTableBC");

			//Constructors

				//Construct from patch and internal field
				movingWaterTableBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct from patch, internal field and dictionary
				movingWaterTableBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);

				//Construct by mapping given movingWaterTableBCFvPatchScalarField onto a new patch
				movingWaterTableBCFvPatchScalarField
				(
					const movingWaterTableBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);

				//Construct as copy
				movingWaterTableBCFvPatchScalarField
				(
					const movingWaterTableBCFvPatchScalarField&
				);

				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new movingWaterTableBCFvPatchScalarField(*this)
					);
				}

				//Construct as copy setting internal field reference
				movingWaterTableBCFvPatchScalarField
				(
					const movingWaterTableBCFvPatchScalarField&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct and return a clone setting internal field reference
				virtual tmp<fvPatchScalarField> clone
				(
					const DimensionedField<scalar, volMesh>& iF
				) const
				{
					return tmp<fvPatchScalarField>
					(
						new movingWaterTableBCFvPatchScalarField(*this, iF)
					);
				}

			//Member functions

				//Update Coefficients
				virtual void updateCoeffs();

				//Write
				virtual void write(Ostream&) const;
				
				//Return the Initial height of the water table
				const scalarField& zWT() const
				{
					return zWT_;
				}
				
				scalarField& zWT()
				{
					return zWT_;
				}
				
				//Return the Velocity of the moving water table
				const scalarField& vWT() const
				{
					return vWT_;
				}
				
				scalarField& vWT()
				{
					return vWT_;
				}

				//Mapping functions

					//Map (and resize as needed) from self given a mapping object
					virtual void autoMap
					(
						const fvPatchFieldMapper&
					);

					//Reverse map the given fvPatchField onto this fvPatchField
					virtual void rmap
					(
						const fvPatchScalarField&,
						const labelList&
					);

	};

} // End namespace Foam

#endif
