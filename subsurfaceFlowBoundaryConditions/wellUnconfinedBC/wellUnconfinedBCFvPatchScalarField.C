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
#include "wellUnconfinedBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "solveTan.H"
#include "besselK0.H"
#include "besselK1.H"
#include "stehfestCoefficients.H"

namespace Foam
{
	//Constructors
	
	wellUnconfinedBCFvPatchScalarField::wellUnconfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		Kr_(0.0),
		Kz_(0.0),
		Ss_(0.0),
		Sy_(0.0),
		Qp_(0.0),
		rW_(0.0),
		r_(0.0),
		zWT_(0.0),
		d_(0.0),
		l_(0.0),
		dS_(0.0),
		Ksw_(0.0),
		rWC_(0.0)
		{}
	/*--------------------------------------------------------------------------*/
	wellUnconfinedBCFvPatchScalarField::wellUnconfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		Kr_(readScalar(dict.lookup("Kr"))),
		Kz_(readScalar(dict.lookup("Kz"))),
		Ss_(readScalar(dict.lookup("Ss"))),
		Sy_(readScalar(dict.lookup("Sy"))),
		Qp_(readScalar(dict.lookup("Qp"))),
		rW_(readScalar(dict.lookup("rW"))),
		r_(readScalar(dict.lookup("r"))),
		zWT_(readScalar(dict.lookup("zWT"))),
		d_(readScalar(dict.lookup("d"))),
		l_(readScalar(dict.lookup("l"))),
		dS_(readScalar(dict.lookup("dS"))),
		Ksw_(readScalar(dict.lookup("Ksw"))),
		rWC_(readScalar(dict.lookup("rWC")))
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
	wellUnconfinedBCFvPatchScalarField::wellUnconfinedBCFvPatchScalarField
	(
		const wellUnconfinedBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedValueFvPatchScalarField(ptf, p, iF, mapper),
		Kr_(ptf.Kr_),
		Kz_(ptf.Kz_),
		Ss_(ptf.Ss_),
		Sy_(ptf.Sy_),
		Qp_(ptf.Qp_),
		rW_(ptf.rW_),
		r_(ptf.r_),
		zWT_(ptf.zWT_),
		d_(ptf.d_),
		l_(ptf.l_),
		dS_(ptf.dS_),
		Ksw_(ptf.Ksw_),
		rWC_(ptf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	wellUnconfinedBCFvPatchScalarField::wellUnconfinedBCFvPatchScalarField
	(
		const wellUnconfinedBCFvPatchScalarField& wucpsf
	)
	:
		fixedValueFvPatchScalarField(wucpsf),
		Kr_(wucpsf.Kr_),
		Kz_(wucpsf.Kz_),
		Ss_(wucpsf.Ss_),
		Sy_(wucpsf.Sy_),
		Qp_(wucpsf.Qp_),
		rW_(wucpsf.rW_),
		r_(wucpsf.r_),
		zWT_(wucpsf.zWT_),
		d_(wucpsf.d_),
		l_(wucpsf.l_),
		dS_(wucpsf.dS_),
		Ksw_(wucpsf.Ksw_),
		rWC_(wucpsf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	wellUnconfinedBCFvPatchScalarField::wellUnconfinedBCFvPatchScalarField
	(
		const wellUnconfinedBCFvPatchScalarField& wucpsf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(wucpsf, iF),
		Kr_(wucpsf.Kr_),
		Kz_(wucpsf.Kz_),
		Ss_(wucpsf.Ss_),
		Sy_(wucpsf.Sy_),
		Qp_(wucpsf.Qp_),
		rW_(wucpsf.rW_),
		r_(wucpsf.r_),
		zWT_(wucpsf.zWT_),
		d_(wucpsf.d_),
		l_(wucpsf.l_),
		dS_(wucpsf.dS_),
		Ksw_(wucpsf.Ksw_),
		rWC_(wucpsf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	//Member Functions
	/*--------------------------------------------------------------------------*/
	void wellUnconfinedBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		constexpr scalar Pi = 3.141592654;
		constexpr scalar maxQ = 708.0;
		constexpr scalar alphaUC = 1e9;
		constexpr label Ns = 8;
		constexpr label Nc_min = 4;
		constexpr label Nc_max = 30;
		
		const fvPatchField<scalar>& z = patch().lookupPatchField<volScalarField, scalar>("z");
		
		scalar t = this->db().time().value();
		
		/*if (t < 1e-6)
		{
			t = 1e-6;
		}*/
		
		//Dimensionless variables
		scalar rD = r_/rW_;									
		scalar rWD = rW_/zWT_;					  				
		scalar dD = d_/zWT_;					  				
		scalar lD = l_/zWT_;									
		scalarField zD = z/zWT_;
							  			 
		scalar gammaUC = alphaUC*zWT_*Sy_/Kz_;  				//Dimensionless empirical drainage constant
		scalar sigmaUC = Ss_*zWT_/Sy_;
		scalar betaUCD = (Kz_/Kr_)*sqr(rWD);
		scalar betaUC = betaUCD*sqr(rD);
		scalar SW = (Kr_/Ksw_)*(dS_/rW_);						//Dimensionless well-bore skin parameter
		scalar WD = sqr(rWC_)/(2.0*Ss_*(l_ - d_)*sqr(rW_)); 	//Dimensionless well-bore storage parameter 
		scalar tD = Kr_*t/(Ss_*sqr(rW_));						//Dimensionless time
		
		//Number of terms for the finite summation of drawdown
		label Nc = Nc_max*(pow(2.0,(-log10(betaUC)-2.0)));		
		if (Nc < Nc_min)
		{
			Nc = Nc_min;
		}
		else if (Nc > Nc_max)
		{
			Nc = Nc_max;
		}
		
		scalar saWD = 0.0;
		scalarField saOD = 0.0*z;
		label ct = 1;
		    
		for (label i = 1; i <= Ns; i++)
		{
			scalar p = ct*log(2.0)/tD;
				
			//Calculation of the Laplacian of the drwadowns at well and observation locations
			scalar ar = p/(sigmaUC*betaUCD + (p/gammaUC));
			scalar eps0 = Pi/2.0;
			scalar A = 0.0;
			scalarField E = 0.0*z;
				
			for(label j = 1; j <= Nc; j++)
			{
				scalar eps = solTan(eps0, ar);
				
				scalar qN1 = sqrt(p + betaUCD*sqr(eps));
				if (qN1 > maxQ)
				{
					qN1 = maxQ;
				}
				
				scalar qN2 = qN1*rD;
				if (qN2 > maxQ)
				{
					qN2 = maxQ;
				}
				
				scalar K0a = besselK0(qN1);
				scalar K0e = besselK0(qN2);
				scalar K1 = besselK1(qN1); 
				
				scalar coef1 = sin(eps*(1.0 - dD)) - sin(eps*(1.0 - lD));
				scalar coef2 = qN1*K1*(eps + 0.5*sin(2.0*eps));
				
				A = A + (2.0/(lD - dD))*(K0a*sqr(coef1)/(eps*coef2));
				E = E + 2.0*(K0e*cos(eps*zD)*coef1/coef2);
					
				eps0 = eps + Pi;  
			} 
			
			scalar lap_saWD = (A + SW)/(p*(1.0 + WD*p*(A + SW)));
			scalarField lap_saOD = E/(p*(1.0 + WD*p*(A + SW)));
			
			scalar wt = stehfestCoef(ct, Ns);

			saWD = saWD + wt*lap_saWD;
			saOD = saOD + wt*lap_saOD;
			
			ct += 1;
		} 
			
		saWD = (log(2.0)/tD)*(2.0/(lD - dD))*saWD;
		saOD = (log(2.0)/tD)*(2.0/(lD - dD))*saOD;
	
		scalar zaW = zWT_ - saWD*Qp_/(4.0*Pi*Kr_*zWT_);
		scalarField zaO = zWT_ - saOD*Qp_/(4.0*Pi*Kr_*zWT_);
		
		Info << "\nHydraulic head at the well =" << zaW << endl;
		
		operator == (zaO - z);
		
		fixedValueFvPatchScalarField::updateCoeffs();		
	}
	/*--------------------------------------------------------------------------*/
	//Write
	void wellUnconfinedBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		os.writeEntry("Kr", Kr_);
		os.writeEntry("Kz", Kz_);
		os.writeEntry("Ss", Ss_);
		os.writeEntry("Sy", Sy_);
		os.writeEntry("Qp", Qp_);
		os.writeEntry("rW", rW_);
		os.writeEntry("r", r_);
		os.writeEntry("zWT", zWT_);
		os.writeEntry("d", d_);
		os.writeEntry("l", l_);
		os.writeEntry("dS", dS_);
		os.writeEntry("Ksw", Ksw_);
		os.writeEntry("rWC", rWC_);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, wellUnconfinedBCFvPatchScalarField);
			
} //End namespace Foam
