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

\*---------------------------------------------------------------------------*/

#include "generalManualReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalManualReaction, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	generalManualReaction,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
generalManualReaction::generalManualReaction
(
    const dictionary& dict,
    cfdemCloud&       sm
)
:
    generalManual(dict,sm),
    haveReaction_(false),
    reactionDict_(propsDict_.subDict("ReactionParameters")),
    reactantsList_(reactionDict_.lookup("reactants")),
    reactantsStoichiometricFactors_(reactionDict_.lookup("stoichiometricFactors")),
    reactantsExponents_(reactionDict_.lookup("exponents")),
    reactionRateConstant_(reactionDict_.lookupOrDefault<dimensionedScalar>("reactionRateConstant", 1.0)),
    mReactionRate_
    (   IOobject
        (
            reactionDict_.lookupOrDefault<word>("reactionRateFieldName","generalManualReactionRate"),
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("dummy", dimensionSet(0,0,-1,0,0), 0)
    )
{
    if(reactantsList_.size()) 
        haveReaction_ =true;
    else
        return;
    

    Info << "generalManualReaction is now setting the reaction parameters for " 
         <<  reactantsList_.size()
         << " reactants." << endl;

    if( reactantsList_.size()!= reactantsStoichiometricFactors_.size() ||
        reactantsList_.size()!= reactantsExponents_.size())
        FatalError <<"The list of reactants, exponents, and stoichiometric factors is of different length. " 
                   << abort(FatalError);  
    for (int i=0;i<reactantsList_.size();i++)
    {
        bool foundIndex = false;
        for (int j=0;j<eulerianFieldList_.size();j++)
        {
            if(eulerianFieldList_[j] == reactantsList_[i])    
            {
                foundIndex = true;
                reactantsIndex_.push_back(j);
            }
        }
        if(!foundIndex)   FatalError <<"generalManualReaction cannot find your reactant '" 
                                     << reactantsList_[i] << ". Please check name of reactant" << endl
                                     << abort(FatalError);  
        Info << "...for reactant " 
             << reactantsList_[i] 
             << " the index: " 
             << reactantsIndex_.back()
             << ", the exponent " << reactantsExponents_[i]
             << ", and the stoichiometric factors " << reactantsStoichiometricFactors_[i]
             << " was identified." << endl;
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
generalManualReaction::~generalManualReaction()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// ************************************************************
void generalManualReaction::update()
{
    //Re-set the sources due to particle-fluid interactions
    generalManual::setSources();

    //Calculate the reaction rates and modify the sources    
    calculateReactionRate();
    modifySourceTerms();

    //to stuff for standard scalar transport
    generalManual::evolveFields();
}

// ************************************************************
void generalManualReaction::calculateReactionRate()
{
    if(!haveReaction_)  return;

    //Calculate the reaction rate (this is per CONTINUOUS unit volume, NOT per TOTAL unit volume since reaction is only in CONTINUOUS phase)
    //This uses a simple Pi-Product formulation
    mReactionRate_ = reactionRateConstant_ 
                   * Foam::pow(eulerianScalarF(reactantsIndex_[0]).m(), reactantsExponents_[0]);
    for (int i=1;i<reactantsList_.size();i++) //start from 1 since first elements resets th rate
    {
        if(reactantsExponents_[i]!=0)
            mReactionRate_ *= Foam::pow(eulerianScalarF(reactantsIndex_[i]).m(), reactantsExponents_[i]);
    }
}

// ************************************************************
void generalManualReaction::modifySourceTerms()
{
    if(!haveReaction_)  return;

    //Set the sources due to particle-fluid interactions
    //This MUST be per TOTAL unit volume, so we need to multiply with the voidage)
    const volScalarField&     voidfraction(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));
    for (int i=0;i<reactantsList_.size();i++) 
    {
        eulerianScalarF(reactantsIndex_[i]).mSource() += voidfraction 
                                                       * reactantsStoichiometricFactors_[i] * mReactionRate_;
    }
}

} // End namespace Foam

// ************************************************************************* //
