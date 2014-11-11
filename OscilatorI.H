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
Credits
    Hrvoje Jasak
    Dubravko Matijasevic

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::Oscilator::mass() const
{
    return mass_;
}


const Foam::dimensionedVector& Foam::Oscilator::Xrel() const
{
    return Xrel_;
}


Foam::dimensionedVector Foam::Oscilator::X() const
{
    return xEquilibrium_ + Xrel();
}


const Foam::dimensionedVector& Foam::Oscilator::U() const
{
    return U_;
}


const Foam::dimensionedVector& Foam::Oscilator::Uold() const
{
    return Uold_;
}


Foam::dimensionedVector Foam::Oscilator::A() const
{
    return A(Xrel(), U());
}


Foam::dimensionedVector Foam::Oscilator::Uaverage() const
{
    return 0.5*(U_ + Uold_);
}


const Foam::dimensionedVector& Foam::Oscilator::force() const
{
    return force_;
}


Foam::dimensionedVector& Foam::Oscilator::force()
{
    return force_;
}


// ************************************************************************* //
