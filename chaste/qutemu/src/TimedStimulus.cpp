#include "TimedStimulus.hpp"

double TimedStimulus::GetStimulus(double time) {
    if (time == mLastTime)
        return mLastStim;

    if (time < mLastTime)
        mStimIndex = 0;

    mLastTime = time;

    for (unsigned i = mStimIndex; i < mTimes.size(); i++) {
        double startTime = mTimes[i];
        if (startTime >= time)
            break;

        mStimIndex = i;
        if (time >= startTime && time <= startTime + mDuration)
            return mLastStim = mMagnitudeOfStimulus;
    }

    return mLastStim = 0;
}
