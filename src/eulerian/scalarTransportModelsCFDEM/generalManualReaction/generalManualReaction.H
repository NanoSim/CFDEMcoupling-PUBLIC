/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2014- Stefan Radl, TU Graz, Austria

    Description
    Evolves an arbitrary number of scalar quantities, exchanges sources
    with the granular phase

\*---------------------------------------------------------------------------*/

#ifndef generalManualReaction_H
#define generalManualReaction_H

#include "forceModel.H"
#include "generalManual.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalManual Declaration
\*---------------------------------------------------------------------------*/

class generalManualReaction
:
    public generalManual
{

protected:

    bool                            haveReaction_; 

    dictionary                      reactionDict_;

    const wordList                  reactantsList_;     //the names of the reactants for easy selection

    mutable std::vector<int>        reactantsIndex_;    //saves the index of the reactants

    const scalarList                reactantsStoichiometricFactors_;    

    const scalarList                reactantsExponents_;
    
    const dimensionedScalar         reactionRateConstant_;

    mutable volScalarField          mReactionRate_;       //the reaction rate


public:

    //- Runtime type information
    TypeName("generalManualReaction");


    // Constructors

    //- Construct from components
    generalManualReaction
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

    ~generalManualReaction();


    // Member Functions

    void                 update();

    void                 calculateReactionRate();

    void                 modifySourceTerms();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
