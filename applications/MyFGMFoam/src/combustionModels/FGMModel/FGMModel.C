/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FGMModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGMModel<ReactionThermo>::FGMModel
(
    const word& modelType, ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    solver_(tableSolver(this->mesh(), tables())),
    PVsolver_(PVtableSolver(this->mesh(), PVtables())),
    Sct_(this->coeffs().lookupOrDefault("Sct", 0.8)),  // Liu	
     
    T_(this->thermo().T()),                 //--- senbin 
    psi_(const_cast<volScalarField&>(this->thermo().psi())),
    alpha_(const_cast<volScalarField&>(this->thermo().alpha())),           //--- senbin 

    mu_(const_cast<volScalarField&>(this->thermo().mu()())),                  
    sourceYc_(const_cast<volScalarField&>(this->turbulence().sourceYc())),                     
    YWI_(const_cast<volScalarField&>(this->turbulence().YWI())),
    YuWI_(const_cast<volScalarField&>(this->turbulence().YuWI())),
    YbWI_(const_cast<volScalarField&>(this->turbulence().YbWI())),
    Cmean_(const_cast<volScalarField&>(this->turbulence().Cmean())),  
    Cfluc_(const_cast<volScalarField&>(this->turbulence().Cfluc())),                
    Zvar_(const_cast<volScalarField&>(this->turbulence().Zvar())),
    Cvar_(const_cast<volScalarField&>(this->turbulence().Cvar())),   
    Zmean_(const_cast<volScalarField&>(this->turbulence().Zmean())),
    Zfluc_(const_cast<volScalarField&>(this->turbulence().Zfluc())),
    Ycmean_(const_cast<volScalarField&>(this->turbulence().Ycmean())),                
    Ycfluc_(const_cast<volScalarField&>(this->turbulence().Ycfluc())),          
    
    useProgressVariableVariance_(this->coeffs().lookupOrDefault("useProgressVariableVariance",true))      
{
    // Sct_ = dimensioned<scalar>::lookupOrDefault  //Liu
    // (
    //     "Sct",
    //     this->coeffDict(),
    //     1.0
    // );

    findUscaledPV();
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGMModel<ReactionThermo>::~FGMModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
hashedWordList FGMModel<ReactionThermo>::tables()
{
    // - Read useProgressVariableVariance_  
    read();

	// - FGM tables to be read
	hashedWordList tableNames;
	tableNames.clear();
	tableNames.append("T");   //--- Read 'T' table
	tableNames.append("psi");
	tableNames.append("mu");
	tableNames.append("alpha");
	tableNames.append("sourceYc"); //--- Read source term of PV
        
    if (useProgressVariableVariance_)
	{
           tableNames.append("YWI"); //--- Read source term of PV
           tableNames.append("YuWI"); //--- To calculate scalar dissipation rate for Yc
           tableNames.append("YbWI"); //--- To calculate scalar dissipation rate for Yc   
	}
        return tableNames;
}

template<class ReactionThermo>
hashedWordList FGMModel<ReactionThermo>::PVtables()
{
    // - Read useProgressVariableVariance_  
    read();

    // - PV tables to be read
	hashedWordList PVtableNames;
	PVtableNames.clear();
	PVtableNames.append("minYc");
	PVtableNames.append("maxYc");
    
	if (useProgressVariableVariance_)
	{
	   PVtableNames.append("Yu2I");
	   PVtableNames.append("YuYbI");
	   PVtableNames.append("Yb2I");
	}
       return PVtableNames;
}


template<class ReactionThermo>
void FGMModel<ReactionThermo>::correct()
{
    Info << "Entering FGMModel correct()" << endl;    //--- Added L.Ma, 14-10-2014
    const scalarField& ZmeanCells = Zmean_.internalField();  
    const scalarField& ZflucCells = Zfluc_.internalField(); 
    const scalarField& YcmeanCells = Ycmean_.internalField();
    const scalarField& YcflucCells = Ycfluc_.internalField();
    
    scalarField& TCells = T_.ref();
    scalarField& sourceYcCells = sourceYc_.ref();       
    scalarField& CmeanCells = Cmean_.ref();       
    scalarField& YWICells = YWI_.ref();               
    scalarField& YuWICells = YuWI_.ref();              
    scalarField& YbWICells = YbWI_.ref();             
    scalarField& ZvarCells = Zvar_.ref();           
    scalarField& CvarCells = Cvar_.ref();     
    scalarField& CflucCells = Cfluc_.ref();  
    scalarField& muCells = mu_.ref();                   
    scalarField& alphaCells = alpha_.ref();             
    scalarField& psiCells = psi_.ref();             

    //- Update the species and temperature field
    if(this->active())
    {
       //- 2D table for minimum and maximum progress variable, (varPV, PV)
       scalarList params2D(2, 0.0); 
           
       //- 4D independent variables (varPV, PV, varZ, Z)
       scalarList params4D(4, 0.0);
       
       //- Upper bounds and position for table interpolation
       List<int>  ubParams4dIF_, ubParams2dIF_, ubParams4dBF_, ubParams2dBF_;     //senbin 06-11-18         //--- Added L.Ma, 06-10-2014
       scalarList posParams4dIF_, posParams2dIF_, posParams4dBF_, posParams2dBF_;  //senbin 06-11-18          //--- Added L.Ma, 06-10-2014

       //- Minimum and Maximum unscaled progress variable
       scalar YcMinCells, YcMaxCells, YcMinPatchs, YcMaxPatchs;
       scalar Cmean_i, Cfluc_i;
       //- Scale progress variable and its variance
       scalar fc, gc, hc;
       scalar Yu2I, YuYbI, Yb2I;
       
       //- A small number to prevent divide by zero
       scalar smallVar(1e-5);

       forAll(ZmeanCells, cellI)    // 1st: for internal field
       {

            //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as a independent variable
            ZvarCells[cellI] = max(min(ZflucCells[cellI]/max(ZmeanCells[cellI]*(1.0 - ZmeanCells[cellI]), smallVar),0.99), 0.0);
                   
            //- Calculate scaled progress variable   
	        params2D[0] = ZvarCells[cellI];
            params2D[1] = max(min(ZmeanCells[cellI], 1.0), 0.0);
      
            ubParams2dIF_ = PVsolver_.upperBounds(params2D);
            posParams2dIF_ = PVsolver_.position(ubParams2dIF_, params2D);
                 
            //- Update minimum and maximum unscaled progress variable
            YcMinCells = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 0);
            YcMaxCells = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 1);
                
            //- Update scaled progress variable
            Cmean_i = (YcmeanCells[cellI] - YcMinCells) / max((YcMaxCells - YcMinCells), smallVar);
            CmeanCells[cellI] =  max(min(Cmean_i, 1.0), 0.0);    // = C''2

             //- Update scaled progress variable variance?
            if (useProgressVariableVariance_)
            {
                //- calculate scaled progress variable variance from unscaled propress variance
                //- Reference 1: A progress variable approach based on premixed flamelets for turbulent combustion modeling, B.A. Albrecht, W.J.S.Ramaekers et al
                //- Eq.(10) and Eq.(17)  (in numerator should be -2gc)
                //- Reference 2: A premixed Flamelet-PDF Model for Biomass Combustion in a Grate Furnace, Energy & Fuels, 2008, Albrecht, Oijen et al
                //- Eq.(10) and Eq.(11) (in denomenator should be (Yb-Yu)^2)
                //- Reference 3: Private communication, Scaling of the reaction progress variable
                
                Yu2I = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 2);
                YuYbI = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 3);
                Yb2I = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 4);
                
                fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2.0*YuYbI + Yu2I;
                Cfluc_i =  (YcflucCells[cellI] + sqr(YcmeanCells[cellI]) - fc - 2.0*gc*CmeanCells[cellI]) / max(hc, smallVar) - sqr(CmeanCells[cellI]);
                CflucCells[cellI] = max(Cfluc_i, 0.0);
                if (hc < smallVar)        
                {
                    CflucCells[cellI] = 0;
                }     
            }
            else
            {
                CflucCells[cellI] = 0.0;
            }
                   
            CvarCells[cellI] = max(min(CflucCells[cellI]/max(CmeanCells[cellI]*(1.0 -CmeanCells[cellI]), smallVar), 0.99), 0.0);
                                    
            //- Prepare independent variables for table look up          
             params4D[0] = CvarCells[cellI];
             params4D[1] = CmeanCells[cellI];
	         params4D[2] = max(min(ZvarCells[cellI], 0.99), 0.0);
             params4D[3] = max(min(ZmeanCells[cellI], 1.0), 0.0);
         
	        //- find up-bound and pos for FGM table interpolation
             ubParams4dIF_ = solver_.upperBounds(params4D);
             posParams4dIF_ = solver_.position(ubParams4dIF_, params4D);
 
             //- Update temperature, sourcePV, flame sensor internal field
             TCells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 0);               
	         psiCells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 1);     
             muCells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 2);
	         alphaCells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 3);
	         sourceYcCells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 4);
	     
            if (useProgressVariableVariance_)
            {
                YWICells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 5);
                YuWICells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 6);
                YbWICells[cellI] = solver_.interpolate(ubParams4dIF_, posParams4dIF_, 7);
            }
	        else
            {
                YWICells[cellI] = 0.0;
                YuWICells[cellI] = 0.0;
                YbWICells[cellI] = 0.0;
	        }  
        }

        // Interpolate for patches
        forAll(T_.boundaryField(), patchi)    // Changed L.Ma, 24-06-2014, 2nd, for patches.
        {
            const fvPatchScalarField& pYcfluc = Ycfluc_.boundaryField()[patchi];
	        const fvPatchScalarField& pYcmean = Ycmean_.boundaryField()[patchi];
            const fvPatchScalarField& pZfluc = Zfluc_.boundaryField()[patchi];
            const fvPatchScalarField& pZmean = Zmean_.boundaryField()[patchi];
	      
            fvPatchScalarField& pT = T_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014 
            fvPatchScalarField& psourceYc = sourceYc_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
            fvPatchScalarField& pCmean = Cmean_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
            fvPatchScalarField& pYWI = YWI_.boundaryFieldRef()[patchi];
	        fvPatchScalarField& pYuWI = YuWI_.boundaryFieldRef()[patchi];
    	    fvPatchScalarField& pYbWI = YbWI_.boundaryFieldRef()[patchi];
	        fvPatchScalarField& pZvar = Zvar_.boundaryFieldRef()[patchi];	  
	        fvPatchScalarField& pCfluc = Cfluc_.boundaryFieldRef()[patchi];	  
	  
            fvPatchScalarField& pCvar = Cvar_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
            fvPatchScalarField& pMu = mu_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
            fvPatchScalarField& pAlpha = alpha_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
            fvPatchScalarField& ppsi = psi_.boundaryFieldRef()[patchi];
          
            forAll(pZmean, facei)
            {
      
               //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as a independent variable
               pZvar[facei] = max(min(pZfluc[facei]/max(pZmean[facei]*(1.0 - pZmean[facei]), smallVar), 0.99), 0.0);

               //- Calculate scaled progress variable 
               params2D[0] = pZvar[facei];
               params2D[1] = max(min(pZmean[facei], 1.0), 0.0);
      
               ubParams2dBF_ = PVsolver_.upperBounds(params2D);
               posParams2dBF_ = PVsolver_.position(ubParams2dBF_, params2D);
      
               //- Update minimum and maximum unscaled progress variable
               YcMinPatchs = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 0);
               YcMaxPatchs = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 1);
      
               //- Update scaled progress variable
	           Cmean_i =  (pYcmean[facei] - YcMinPatchs)/max((YcMaxPatchs - YcMinPatchs), smallVar);
               pCmean[facei] = max(min(Cmean_i, 1.0), 0.0);

               //- How about the progress variable variance?
                if (useProgressVariableVariance_)
                {    
                    Yu2I = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 2);
                    YuYbI = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 3);
                    Yb2I = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 4);
                    
                    fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2.0*YuYbI + Yu2I;
                    Cfluc_i =  (pYcfluc[facei] + sqr(pYcmean[facei]) - fc - 2.0*gc*pCmean[facei])/max(hc, smallVar) - sqr(pCmean[facei]);
                    
                    pCfluc[facei] =  max(Cfluc_i, 0.0);
                            
                    //- For the situation where Pvmin = PVmax, this normally happens at the fuel or oxidizer inlet
                    //- The scaled progress variable is set to zero in this situation
                    if (hc < smallVar)        
                        {
                            pCfluc[facei] = 0;
                        }		       
                }
                else
                {
                    pCfluc[facei] = 0.0;
                }
		     
                pCvar[facei] = max(min(pCfluc[facei]/max(pCmean[facei]*(1.0 - pCmean[facei]), smallVar), 0.99), 0.0);
                    
                //- Prepare independent variables for table look up
                params4D[0] = pCvar[facei];
                params4D[1] = pCmean[facei];
                params4D[2] = max(min(pZvar[facei], 0.99), 0.0);
                params4D[3] = max(min(pZmean[facei], 1.0), 0.0);
                ubParams4dBF_ = solver_.upperBounds(params4D);
                posParams4dBF_ = solver_.position(ubParams4dBF_, params4D);

                //- Updated temperature patch field
                pT[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 0);
                ppsi[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 1);  
                pMu[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 2);
                pAlpha[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 3);
                psourceYc[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 4);
                
                if (useProgressVariableVariance_)
                {     
                        pYWI[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 5);
                        pYuWI[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 6);
                        pYbWI[facei] = solver_.interpolate(ubParams4dBF_, posParams4dBF_, 7);
                }
                else
                {
                        pYWI[facei] = 0.0;
                        pYuWI[facei] = 0.0;
                        pYbWI[facei] = 0.0;
                }
            }
        }
  
    }
}

template<class ReactionThermo>
Switch FGMModel<ReactionThermo>::correctDensity()
{
	return true;
}
//---------------------------------------------------------------
template<class ReactionThermo>
void FGMModel<ReactionThermo>::findUscaledPV()
{
    Info << endl;
    Info << "Initializing unscaled progress variable" << endl;
    Info << endl;

    List<int> ubParams2dIF_, ubParams2dBF_;  // IF = internal Field; BF = Boundary Field
    scalarList posParams2dIF_, posParams2dBF_;
    scalar YcMinIF, YcMaxIF, YcMinBF, YcMaxBF;
    scalar varZIF, varZBF;
    scalarList params2D(2, 0.0);
    
    const scalarField& meanZIF = Zmean_.internalField();
    const scalarField& flucZIF = Zfluc_.internalField();  
    const scalarField& meanCIF = Cmean_.internalField(); 
    scalarField& meanYcIF = Ycmean_.ref();  // L
    
    forAll(meanZIF, cellI)
    {
        varZIF = flucZIF[cellI]/max(meanZIF[cellI]*(1.0 - meanZIF[cellI]), SMALL);

        params2D[0] = max(min(varZIF, 0.99), 0.0);
	    params2D[1] = max(min(meanZIF[cellI], 1.0), 0.0);
	
        ubParams2dIF_ = PVsolver_.upperBounds(params2D);
        posParams2dIF_ = PVsolver_.position(ubParams2dIF_, params2D);

        //- Find minimum and maximum unscaled progress variable
        YcMinIF = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 0);
        YcMaxIF = PVsolver_.interpolate(ubParams2dIF_, posParams2dIF_, 1);
	
        meanYcIF[cellI]  =  YcMinIF + meanCIF[cellI]*(YcMaxIF-YcMinIF);
    }
    
    forAll(Zmean_.boundaryField(), patchi)    
    {
        const fvPatchScalarField& meanZBF = Zmean_.boundaryField()[patchi];
        const fvPatchScalarField& flucZBF = Zfluc_.boundaryField()[patchi];      
        const fvPatchScalarField& meanCBF = Cmean_.boundaryField()[patchi];  
        fvPatchScalarField& meanYcBF = Ycmean_.boundaryFieldRef()[patchi];        
        
        forAll(meanZBF , facei)
        {
            varZBF = flucZBF[facei]/max(meanZBF[facei]*(1.0 - meanZBF[facei]), SMALL);
            //- Mixture fraction
            params2D[0] = max(min(varZBF, 0.99), 0.0);
            params2D[1] = max(min(meanZBF[facei], 1.0), 0.0);

            ubParams2dBF_ = PVsolver_.upperBounds(params2D);
            posParams2dBF_ = PVsolver_.position(ubParams2dBF_, params2D);

            //- Find minimum and maximum unscaled progress variable
            YcMinBF = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 0);
            YcMaxBF = PVsolver_.interpolate(ubParams2dBF_, posParams2dBF_, 1);

            meanYcBF[facei] = YcMinBF + meanCBF[facei]*(YcMaxBF - YcMinBF);
        
        }
    }
}
//\-------------------------------------------------------------------

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
FGMModel<ReactionThermo>::R
(
    volScalarField& Y                 //---Changed L.Ma, 03-07-2014
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

// Progress variable source 
template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGMModel<ReactionThermo>::SourceYc() const          //Added L.Ma, 08-10-2014
{
  return sourceYc_;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGMModel<ReactionThermo>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tSh;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGMModel<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimVolume/dimTime, 0.0), //here the unit is right for HRR
            zeroGradientFvPatchScalarField::typeName
        )
    );

}

template<class ReactionThermo>
bool FGMModel<ReactionThermo>::read()
{
    if (ThermoCombustion<ReactionThermo>::read())
    {
        this->coeffs().lookup("useProgressVariableVariance") >> useProgressVariableVariance_; 
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
