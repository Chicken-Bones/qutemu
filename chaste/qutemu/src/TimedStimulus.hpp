#pragma once

#include "AbstractStimulusFunction.hpp"
#include <vector>

class TimedStimulus : public AbstractStimulusFunction {
public:
    /** The 'height' of the square wave applied */
    double mMagnitudeOfStimulus;
    /** The length of the square wave */
    double mDuration;
    /** The activation times */
    std::vector<double> mTimes;

private:
    unsigned mStimIndex; ///Optimisation
    double mLastTime; ///Cache variable
    double mLastStim; ///Cache variable

public:
    TimedStimulus(double magnitudeOfStimulus, double duration, std::vector<double> times) :
            mMagnitudeOfStimulus(magnitudeOfStimulus),
            mDuration(duration),
            mTimes(times),
            mStimIndex(0),
            mLastTime(0),
            mLastStim(0)
    {}

    double GetStimulus(double time) override;
};