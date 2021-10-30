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
    Foam::wellUnconfinedBCFvPatchScalarField

Group
    wellUnconfinedBC/subsurfaceFlowBoundaryConditions

Description
    This boundary condition calculates the drawdown at the well-boundary
    and at the pumping well for unsteady groundwater flow through 
    an unconfined aquifer.
\*-----------------------------------------------------------------------*/
#ifndef wellUnconfinedBCFvPatchScalarField_H
#define wellUnconfinedBCFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{
class wellUnconfinedBCFvPatchScalarField
:
	public fixedValueFvPatchScalarField
	{
		//Private Data
		scalar Kr_; 			//Radial hydraulic conductivity
		scalar Kz_; 			//Vertical hydraulic conductivity
		scalar Ss_; 			//Specific Storage
		scalar Sy_; 			//Specific Yield
		scalar Qp_; 			//Volumetric pumping rate
		scalar rW_; 			//Radius of the pumping well
		scalar r_; 			    //Radius of the bounding circle around the well
		scalar zWT_;  			//Initial saturated thickness
		scalar d_;  			//Top of screened level of pumping well
		scalar l_;  			//Bottom of screened level of pumping well
		scalar dS_; 			//Well-bore skin thickness
		scalar Ksw_; 			//Well-bore skin hydraulic conductivity
		scalar rWC_; 			//Inside radius of the pumping well
		
		public:

			//Runtime type information
			TypeName("wellUnconfinedBC");
			
			//Constructors
				
				//Construct from patch and internal field
				wellUnconfinedBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);
				
				//Construct from patch,internal field and dictionary
				wellUnconfinedBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);
				
				//Construct by mapping given wellUnconfinedBCFvPatchScalarField onto a new patch
				wellUnconfinedBCFvPatchScalarField
				(
					const wellUnconfinedBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);
				
				//Construct as copy
				wellUnconfinedBCFvPatchScalarField
				(
					const wellUnconfinedBCFvPatchScalarField&
				);
				
				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new wellUnconfinedBCFvPatchScalarField(*this)
					);
				}
				
				//Construct as copy setting internal field reference
				wellUnconfinedBCFvPatchScalarField
				(
					const wellUnconfinedBCFvPatchScalarField&,
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
						new wellUnconfinedBCFvPatchScalarField(*this, iF)
					);
				}
				
			//Member Functions
			
		        //Return Radial hydraulic conductivity
		        const scalar& Kr() const
				{
					return Kr_;
				}
		        scalar& Kr()
				{
					return Kr_;
				}
				
				//Return Vertical hydraulic conductivity
				const scalar& Kz() const
				{
					return Kz_;
				}
		        scalar& Kz()
				{
					return Kz_;
				}
				
				//Return Specific Storage
				const scalar& Ss() const
				{
					return Ss_;
				}
				scalar& Ss()
				{
					return Ss_;
				}
				
				//Return Specific Yield
				const scalar& Sy() const
				{
					return Sy_;
				}
				scalar& Sy()
				{
					return Sy_;
				}
		        
		        //Return Volumetric pumping rate
		        const scalar& Qp() const
				{
					return Qp_;
				}
				scalar& Qp()
				{
					return Qp_;
				}
			
				//Return Radius of the pumping well
				const scalar& rW() const
				{
					return rW_;
				}
				scalar& rW()
				{
					return rW_;
				}
				
				//Return Radius of the bounding circle around the well
				const scalar& r() const
				{
					return r_;
				}
				scalar& r()
				{
					return r_;
				}
		
				//Return Initial saturated thickness
				const scalar& zWT() const
				{
					return zWT_;
				}
				scalar& zWT()
				{
					return zWT_;
				}
		
		        //Return Top of screened level of pumping well
		        const scalar& d() const
				{
					return d_;
				}
				scalar& d()
				{
					return d_;
				}
		
		        //Return Bottom of screened level of pumping well
		        const scalar& l() const
				{
					return l_;
				}
		        scalar& l()
				{
					return l_;
				}
				
		        //Return Well-bore skin thickness
		        const scalar& dS() const
				{
					return dS_;
				}
		        scalar& dS()
				{
					return dS_;
				}
				
				//Return Well-bore skin hydraulic conductivity
		        const scalar& Ksw() const
				{
					return Ksw_;
				}
		        scalar& Ksw()
				{
					return Ksw_;
				}
		
				//Return Inside radius of the pumping well
				const scalar& rWC() const
				{
					return rWC_;
				}
				scalar& rWC()
				{
					return rWC_;
				}
				
				//Update Coefficients
				virtual void updateCoeffs();
				
				//Write
				virtual void write(Ostream&) const;
	};
	
} // End namespace Foam

#endif
