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
    This boundary condition calculates the normal gradient of pressure head 
    due to time constant specified flux value at that boundary.
\*-----------------------------------------------------------------------*/
#ifndef specifiedConstantFluxBCFvPatchScalarField_H
#define specifiedConstantFluxBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{
class specifiedConstantFluxBCFvPatchScalarField
:
	public fixedGradientFvPatchScalarField
	{
		//Private Data
		vectorField qC_; //Time Constant Flux
		
		public:

			//Runtime type information
			TypeName("specifiedConstantFluxBC");
			
			//Constructors
				
				//Construct from patch and internal field
				specifiedConstantFluxBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);
				
				//Construct from patch,internal field and dictionary
				specifiedConstantFluxBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);
				
				//Construct by mapping given specifiedConstantFluxBCFvPatchScalarField onto a new patch
				specifiedConstantFluxBCFvPatchScalarField
				(
					const specifiedConstantFluxBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);
				
				//- Construct as copy
				specifiedConstantFluxBCFvPatchScalarField
				(
					const specifiedConstantFluxBCFvPatchScalarField&
				);
				
				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new specifiedConstantFluxBCFvPatchScalarField(*this)
					);
				}
				
				//Construct as copy setting internal field reference
				specifiedConstantFluxBCFvPatchScalarField
				(
					const specifiedConstantFluxBCFvPatchScalarField&,
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
						new specifiedConstantFluxBCFvPatchScalarField(*this, iF)
					);
				}
				
			//Member Functions
			
				//Update Coefficients
				virtual void updateCoeffs();
				
				//Write
				virtual void write(Ostream&) const;
				
				//Return Time Constant Flux
				virtual const vectorField& qC() const
				{
					return qC_;
				}

				virtual vectorField& qC()
				{
					return qC_;
				}
				
	};
	
} // End namespace Foam

#endif
