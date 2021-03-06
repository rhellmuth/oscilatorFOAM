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
    Oscilator
Description
    Ordinary differential equation for three degrees of freedom
    solid body motion. Adapted from translationODE of OF-extend 3.0.
Author
    Rudolf Hellmuth
Credits
    Hrvoje Jasak
    Dubravko Matijasevic

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

#include "Oscilator.H"
#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(Oscilator, 0);
// }

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Oscilator::setCoeffs()
{
    // Set ODE coefficients from position and rotation

    // Linear displacement in relative coordinate system
    {
        const vector& xVal = Xrel_.value();
        coeffs_[0] = xVal.x();
        coeffs_[1] = xVal.y();
        coeffs_[2] = xVal.z();
    }

    // Linear velocity in relative coordinate system
    {
        const vector& uVal = U_.value();
        coeffs_[3] = uVal.x();
        coeffs_[4] = uVal.y();
        coeffs_[5] = uVal.z();
    }
}


Foam::dimensionedVector Foam::Oscilator::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR
) const
{
    return
    (
       - (linSpringCoeffs_ & xR)    // spring
       - (linDampingCoeffs_ & uR)   // damping
       + force()
    )/mass_;                        // gravity
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::Oscilator::Oscilator
(
    const IOobject& io
)
:
    IOdictionary(io),
    mass_(lookup("mass")),
    xEquilibrium_(lookup("equilibriumPosition")),
    linSpringCoeffs_(lookup("linearSpring")),
    linDampingCoeffs_(lookup("linearDamping")),
    Xrel_(lookup("Xrel")),
    U_(lookup("U")),
    Uold_(lookup("Uold")),
    force_(lookup("force")),
    coeffs_(6, 0.0)
{
    setCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Oscilator::~Oscilator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Oscilator::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // Set the derivatives for displacement
    dydx[0] = y[3];
    dydx[1] = y[4];
    dydx[2] = y[5];

    dimensionedVector curX("curX", dimLength, vector(y[0], y[1], y[2]));
    dimensionedVector curU("curU", dimLength/dimTime, vector(y[3], y[4], y[5]));

    const vector& accel = A(curX, curU).value();

    dydx[3] = accel.x();
    dydx[4] = accel.y();
    dydx[5] = accel.z();
}


void Foam::Oscilator::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("Oscilator::jacobian(...) const");
}


void Foam::Oscilator::update(const scalar delta)
{
    // Update position
    vector& Xval = Xrel_.value();

    Xval.x() = coeffs_[0];
    Xval.y() = coeffs_[1];
    Xval.z() = coeffs_[2];

    // Update velocity
    Uold_ = U_;

    vector& Uval = U_.value();

    Uval.x() = coeffs_[3];
    Uval.y() = coeffs_[4];
    Uval.z() = coeffs_[5];
}


bool Foam::Oscilator::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const Oscilator& sds)
{
    os.writeKeyword("mass") << sds.mass_ << token::END_STATEMENT << endl;
    os.writeKeyword("equilibriumPosition") << sds.xEquilibrium_
        << token::END_STATEMENT << endl;
    os.writeKeyword("linearSpring") << sds.linSpringCoeffs_
        << token::END_STATEMENT << endl;
    os.writeKeyword("linearDamping") << sds.linDampingCoeffs_
        << token::END_STATEMENT << endl;

    os.writeKeyword("Xrel") << sds.Xrel() << token::END_STATEMENT << endl;
    os.writeKeyword("U") << sds.U() << token::END_STATEMENT << endl;
    os.writeKeyword("Uold") << tab << sds.Uold() << token::END_STATEMENT << nl;

    os.writeKeyword("force") << sds.force() << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //
