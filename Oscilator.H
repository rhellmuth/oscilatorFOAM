/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    oscilator
Description
    Ordinary differential equation for three degrees of freedom
    solid body motion. Adapted from translationODE of OF-extend 3.0.
Author
    Rudolf Hellmuth
Source Files
    OscilatorI.H
    Oscilator.C
Credits
    Hrvoje Jasak
    Dubravko Matijasevic

\*---------------------------------------------------------------------------*/

#ifndef OSCILATOR_H
#define OSCILATOR_H

#include "ODESystem.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Oscilator Declaration
\*---------------------------------------------------------------------------*/


class Oscilator
:
    public IOdictionary,
    public ODESystem
{
    // Private data

        // Body data

            //- Mass
            dimensionedScalar mass_;


        // Platform variables

            //- Spring equilibrium position for translation
            const dimensionedVector xEquilibrium_;

            //- Linear spring coeffs
            const dimensioned<diagTensor> linSpringCoeffs_;

            //- Linear damping coeffs
            const dimensioned<diagTensor> linDampingCoeffs_;


        // Body position and rotation variables

            //- Displacement relative to spring equilibrium
            dimensionedVector Xrel_;

            //- Velocity of mass centroid
            dimensionedVector U_;

            //- Velocity of mass centroid at previous time-step
            dimensionedVector Uold_;


        // External forces

            //- Force driving the motion
            dimensionedVector force_;


        //- ODE coefficients
        scalarField coeffs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Oscilator(const Oscilator&);

        //- Disallow default bitwise assignment
        void operator=(const Oscilator&);


        //- Set ODE coefficients from position and rotation
        inline void setCoeffs();


        // Variables in relative coordinate system (solved for)

            //- Return acceleration in relative coordinate system
            //  given current values of relative displacement and velocity
            dimensionedVector A
            (
                const dimensionedVector& xR,
                const dimensionedVector& uR
            ) const;


public:

//     //- Runtime type information
//     TypeName("Oscilator");


    // Constructors

        //- Construct from dictionary
        Oscilator(const IOobject& io);


    // Destructor

        virtual ~Oscilator();


    // Member Functions

        //- Return mass
        inline const dimensionedScalar& mass() const;


        // Variables in relative coordinate system (solved for)

            //- Return displacement in relative coordinate system
            inline const dimensionedVector& Xrel() const;


        // Displacement and rotation in the absolute coordinate system

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector X() const;

            //- Return velocity of origin
            inline const dimensionedVector& U() const;

            //- Return velocity of origin for the previous time-step
            inline const dimensionedVector& Uold() const;

            //- Return acceleration of origin
            inline dimensionedVector A() const;


        // Average motion per time-step

            //- Return average velocity of origin
            inline dimensionedVector Uaverage() const;


        // Force

            //- Return force
            inline const dimensionedVector& force() const;

            //- Return access to force
            inline dimensionedVector& force();


        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const
            {
                return 6;
            }

            //- Return reference to interpolation coefficients
            virtual scalarField& coeffs()
            {
                return coeffs_;
            }

            //- Return reference to interpolation coefficients
            virtual const scalarField& coeffs() const
            {
                return coeffs_;
            }

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);


        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const Oscilator&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OscilatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif  // OSCILATOR_H

// ************************************************************************* //
