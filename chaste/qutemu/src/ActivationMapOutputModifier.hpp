#pragma once

#include "AbstractOutputModifier.hpp"

class ActivationMapOutputModifier : public AbstractOutputModifier
{
private:
    class Variable
    {
    public:
        std::string mName;
        std::vector<float> mArr;
        hid_t mVarId;

        Variable(std::string name) : mName(name)
        {}

        float& operator[](int i) {
            return mArr[i];
        }
    };

    double mThresholdVoltage; ///< The user-defined threshold at which activation is to be measured
    double mRestingVoltage; ///< The user-defined resting potiential for APD90 calculation
    std::vector<double> mSnapshotTimes;

    unsigned mNumNodes;    ///< Global problem size
    unsigned mLo;          ///< Local ownership of PETSc node vector
    unsigned mHi;          ///< Local ownership of PETSc node vector
    unsigned mNumberOwned; ///< mNumberOwned=#mHi-#mLo

    hid_t mFileId;

    bool mAnyActivated = false;
    unsigned mActivationIndex = 0; ///< The index of this snapshot, increases every time a cell is reactivated
    double mCurStartTime = 0; ///< The time of the activation that triggered this snapshot. Used for detecting reactivations
    double mLastProcessedTime = 0;
    std::vector<bool> mActivationState; ///< Local per-node vector. True if cell was active last timestep
    std::vector<float> mCurrentPeak; ///< Local per-node vector. Peak value for current activation
    Variable mActivationTime; ///< Local per-node vector. Time of most recent activation
    Variable mPeakVoltage; ///< Local per-node vector. Peak voltage of last activation (reset on threshold cross)
    Variable mActionPotentialDuration; ///< Local per-node vector. APD90 of last repolarisation

    std::vector<Variable*> mVariables;
public:
    ActivationMapOutputModifier(const std::string &rFilename, double thresholdVoltage, double restingVoltage) :
            AbstractOutputModifier(rFilename),
            mThresholdVoltage(thresholdVoltage),
            mRestingVoltage(restingVoltage),
            mActivationTime("Activation"),
            mPeakVoltage("Peak"),
            mActionPotentialDuration("APD"),
            mVariables{&mActivationTime, &mPeakVoltage, &mActionPotentialDuration}
    {};

    ActivationMapOutputModifier(const std::string &rFilename, double thresholdVoltage, double restingVoltage, const std::vector<double> &rSnapshotTimes) :
            ActivationMapOutputModifier(rFilename, thresholdVoltage, restingVoltage)
    {
        mSnapshotTimes = rSnapshotTimes;
    };

    ~ActivationMapOutputModifier() override;

    void InitialiseAtStart(DistributedVectorFactory *pVectorFactory) override;
    void FinaliseAtEnd() override;
    void ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim) override;
    void ProcessPdeSolutionAtTimeStep(double time, Vec solution, unsigned problemDim) override;

private:
    void Close();
    bool IsSnapshotTime(float time, double* pSolution, unsigned problemDim);
    void SaveSnapshot();

    void CreateDataset(Variable* var);
    void SaveDataset(Variable* var);
};
