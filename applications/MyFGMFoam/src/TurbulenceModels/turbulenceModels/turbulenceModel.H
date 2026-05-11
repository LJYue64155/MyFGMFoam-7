/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::turbulenceModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    turbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceModel_H
#define turbulenceModel_H

#include "IOdictionary.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "nearWallDist.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class turbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModel
:
    public IOdictionary
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& alphaRhoPhi_;
        const surfaceScalarField& phi_;

        //- Near wall distance boundary field
        nearWallDist y_;

        //- Mixture Fraction []
        volScalarField Zmean_;

        //- Mixture Fraction variance []
        volScalarField Zfluc_;

        //- Progress variable []
        volScalarField Ycmean_;                  //--- Added L.Ma, 07-10-2014

        //- Progress variable variance []
        volScalarField Ycfluc_;               //--- Added L.Ma, 07-10-2014

        //- Progress variable source []
        volScalarField sourceYc_;            //--- Added L.Ma, 07-10-2014
    
        //- For PV variance source and scalar dissipation rate
        volScalarField YWI_;            //--- Added L.Ma, 06-03-2015

        //- For PV variance source and scalar dissipation rate
        volScalarField YuWI_;            //--- Added L.Ma, 06-03-2015

        //- For PV variance source and scalar dissipation rate
        volScalarField YbWI_;            //--- Added L.Ma, 06-03-2015

        //- Scaled Progress variable[]
        volScalarField Cmean_;            //--- Added L.Ma, 15-10-2014

        //- Scaled Progress variable[]
        volScalarField Cfluc_;            //--- Added L.Ma, 15-10-2014

        //- Normalized Mixture fraction variance[]
        volScalarField Zvar_;            //--- Added L.Ma, 06-03-2015

        //- Normalized Scaled Progress variable variance[]
        volScalarField Cvar_;            //--- Added L.Ma, 06-03-2015

        //- Scalar dissipation rate of mixture fraction[1/s]
        volScalarField ChiZ_;

        //- Scalar dissipation rate of Progress variable[1/s]
        volScalarField ChiC_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulenceModel(const turbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModel&);


public:

    //- Runtime type information
    TypeName("turbulenceModel");

    //- Default name of the turbulence properties dictionary
    static const word propertiesName;


    // Constructors

        //- Construct from components
        turbulenceModel
        (
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const word& propertiesName
        );


    //- Destructor
    virtual ~turbulenceModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        const Time& time() const
        {
            return runTime_;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;

        //- Helper function to return the name of the turbulence G field
        inline word GName() const
        {
            return word(type() + ":G");
        }

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to phase flux field
        inline const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Return the near wall distances
        const nearWallDist& y() const
        {
            return y_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity on patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const = 0;

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const = 0;

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const = 0;

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const = 0;

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const = 0;

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const = 0;

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

            //- Mixture Fraction []
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& Zmean()
	        {
                    return Zmean_; 
                } 
	        //- Mixture Fraction []
	        virtual const volScalarField& Zmean() const
	        {
                    return Zmean_; 
                } 
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& Zfluc()
	        {
                    return Zfluc_; 
                } 

	        //- Mixture Fraction []
	        virtual const volScalarField& Zfluc() const		
		{
			return Zfluc_;
		}		
	        
		//- Progress variable []
	        virtual volScalarField& Ycmean()               //--- Added L.Ma, 07-07-2014
		{
			return Ycmean_;
		}
	        //- Progress variable []
	        virtual const volScalarField& Ycmean() const   //--- Added L.Ma, 07-07-2014
		{
			return Ycmean_;
		}
	
	        //- Progress variable variance []
	        virtual volScalarField& Ycfluc()             //--- Added L.Ma, 07-07-2014
 		{
			return Ycfluc_;
		}	
	        //- Progress variable variance []
	        virtual const volScalarField& Ycfluc() const //--- Added L.Ma, 07-07-2014
 		{
			return Ycfluc_;
		}	
	        //- Progress variable source []
	        virtual volScalarField& sourceYc()           //--- Added L.Ma, 07-07-2014
 		{
			return sourceYc_;
		}	

	        //- Progress variable source []
	        virtual const volScalarField& sourceYc() const   //--- Added L.Ma, 07-07-2014
 		{
			return sourceYc_;
		}	

	        //- For PV source and scalar dissipation rate []
	        virtual volScalarField& YWI()           //--- Added L.Ma, 06-03-2015
 		{
			return YWI_;
		}	
	        //- For PV source and scalar dissipation rate []
	        virtual const volScalarField& YWI() const           //--- Added L.Ma, 06-03-2015
 		{
			return YWI_;
		}	

	        //- For PV source and scalar dissipation rate []
	        virtual volScalarField& YuWI()           //--- Added L.Ma, 06-03-2015
 		{
			return YuWI_;
		}	

	        //- For PV source and scalar dissipation rate []
	        virtual const volScalarField& YuWI() const           //--- Added L.Ma, 06-03-2015
 		{
			return YuWI_;
		}	

	        //- For PV source and scalar dissipation rate []
	        virtual volScalarField& YbWI()           //--- Added L.Ma, 06-03-2015
 		{
			return YbWI_;
		}	

	        //- For PV source and scalar dissipation rate []
	        virtual const volScalarField& YbWI() const           //--- Added L.Ma, 06-03-2015
 		{
			return YbWI_;
		}	

        	//- Scaled Progress variable []
	        virtual volScalarField& Cmean()           //--- Added L.Ma, 15-10-2014
 		{
			return Cmean_;
		}	

	        //- Scaled Progress variable [] 
	        virtual const volScalarField& Cmean() const   //--- Added L.Ma, 15-10-2014
 		{
			return Cmean_;
		}	

	        //- Scaled Progress variable []
	        virtual volScalarField& Cfluc()           //--- Added L.Ma, 15-10-2014
 		{
			return Cfluc_;
		}	

	        //- Scaled Progress variable [] 
	        virtual const volScalarField& Cfluc() const   //--- Added L.Ma, 15-10-2014
 		{
			return Cfluc_;
		}	

	        //- Normalized mixture fraction variance []
	        virtual volScalarField& Zvar()           //--- Added L.Ma, 06-03-2015
 		{
			return Zvar_;
		}	

	        //- Normalized mixture fraction variance [] 
	        virtual const volScalarField& Zvar() const   //--- Added L.Ma, 06-03-2015
 		{
			return Zvar_;
		}	

	        //- Normalized scaled Progress variable variance []
	        virtual volScalarField& Cvar()           //--- Added L.Ma, 06-03-2015
 		{
			return Cvar_;
		}	

	        //- Normalized scaled Progress variable variance [] 
	        virtual const volScalarField& Cvar() const   //--- Added L.Ma, 06-03-2015
 		{
			return Cvar_;
		}	


	        //  Non-const access allowed for transport equations
	        virtual volScalarField& ChiZ() 
 		{
			return ChiZ_;
		}	

	        //- Scalar dissipation rate [1/s]
	        virtual const volScalarField& ChiZ() const 
 		{
			return ChiZ_;
		}	

	        //  Non-const access allowed for transport equations
	        virtual volScalarField& ChiC() 
 		{
			return ChiC_;
		}	

	        //- Scalar dissipation rate [1/s]
	        virtual const volScalarField& ChiC() const 
 		{
			return ChiC_;
		}	

        //- Validate the turbulence fields after construction
        //  Update derived fields as required
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Calculate the mixture fraction variance
        virtual void correctVarZ(); 

        //- Calculate the scalar dissipation rate 
        virtual void correctChiZ(); 

        //- Calculate the scalar dissipation rate 
        virtual void correctChiYc(); 
	
        //- Calculate the progress variable variance
        virtual void correctVarYc();   

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
