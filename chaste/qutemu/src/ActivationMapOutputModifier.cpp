#include <boost/foreach.hpp>
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"

#include "ActivationMapOutputModifier.hpp"

ActivationMapOutputModifier::~ActivationMapOutputModifier() {
    Close();
}

void ActivationMapOutputModifier::InitialiseAtStart(DistributedVectorFactory *pVectorFactory) {
    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
    std::string file_name = output_file_handler.FindFile(mFilename).GetAbsolutePath();

    // Set up a property list saying how we'll open the file
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, PETSC_COMM_WORLD, MPI_INFO_NULL);

    //create file
    hid_t fcpl = H5Pcreate(H5P_FILE_CREATE);
    mFileId = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, fcpl, fapl);
    H5Pclose(fcpl);

    H5Pclose(fapl);

    if (mFileId < 0)
        EXCEPTION("Failed to Create H5F " << file_name << " error code = " << mFileId);

    mNumNodes = pVectorFactory->GetProblemSize();
    mLo = pVectorFactory->GetLow();
    mHi = pVectorFactory->GetHigh();
    mNumberOwned = pVectorFactory->GetLocalOwnership();

    mActivationState.assign(mNumberOwned, false);
    BOOST_FOREACH(Variable* var, mVariables)
        CreateDataset(var);
}

void ActivationMapOutputModifier::CreateDataset(Variable* var) {
    var->mArr.assign(mNumberOwned, std::numeric_limits<float>::quiet_NaN());

    hsize_t data_dims[2] = {1, mNumNodes};
    hsize_t max_dims[2] = {H5S_UNLIMITED, mNumNodes};
    hsize_t chunking[2] = {1, mNumNodes};//one window per chunk, (~1MB for 256k nodes)

    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunking);
    hid_t filespace = H5Screate_simple(2, data_dims, max_dims);
    var->mVarId = H5Dcreate(mFileId, var->mName.c_str(), H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Sclose(filespace);
    H5Pclose(dcpl);
}

void ActivationMapOutputModifier::Close() {
    if (!mFileId)
        return;

    BOOST_FOREACH(Variable* var, mVariables)
        H5Dclose(var->mVarId);

    H5Fclose(mFileId);
    mFileId = 0;
}

void ActivationMapOutputModifier::SaveWindow() {
    BOOST_FOREACH(Variable* var, mVariables)
        SaveDataset(var);

    if (PetscTools::AmMaster())
        std::cout << "Saving Window (" << mCurStartTime << "-" << mLastProcessedTime << ")" << std::endl;
}

void ActivationMapOutputModifier::SaveDataset(Variable* var) {
    hsize_t dims[2] = {mActivationIndex+1, mNumNodes};
    H5Dset_extent(var->mVarId, dims );

    hid_t memspace, hyperslab_space;
    if (mNumberOwned != 0)
    {
        hsize_t v_size[1] = {mNumberOwned};
        memspace = H5Screate_simple(1, v_size, nullptr);

        hsize_t start[2] = {mActivationIndex, mLo};
        hsize_t count[2] = {1, mNumberOwned};

        hyperslab_space = H5Dget_space(var->mVarId);
        H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, nullptr, count, nullptr);
    }
    else
    {
        memspace = H5Screate(H5S_NULL);
        hyperslab_space = H5Screate(H5S_NULL);
    }

    // Create property list for collective dataset write
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    // Write!
    H5Dwrite(var->mVarId, H5T_NATIVE_FLOAT, memspace, hyperslab_space, property_list_id, &(*var)[0]);

    // Tidy up
    H5Sclose(memspace);
    H5Sclose(hyperslab_space);
    H5Pclose(property_list_id);
}

void ActivationMapOutputModifier::FinaliseAtEnd() {
    SaveWindow();
    Close();
}

void ActivationMapOutputModifier::ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim) {
    if (time <= mLastProcessedTime)
        return;
    mLastProcessedTime = time;

    double* p_solution;
    VecGetArray(solution, &p_solution);

    unsigned new_window = false;
    for (unsigned local_index=0; local_index < mNumberOwned; local_index++)
    {
        double v = p_solution[local_index*problemDim];
        if (!mActivationState[local_index] && v > mThresholdVoltage && //activation
                mActivationTime[local_index] >= mCurStartTime)//reactivation (NaN comparison returns false)
        {
            new_window = true;
            std::cout << "Window trigger." <<
                      " node: " << mLo + local_index <<
                      " on " << PetscTools::GetMyRank() <<
                      ", period: " << time - mActivationTime[local_index] << std::endl;
            break;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &new_window, 1, MPI_UNSIGNED, MPI_LOR, PETSC_COMM_WORLD);

    if (new_window) {
        SaveWindow();
        mActivationIndex++;
        mCurStartTime = time;
    }

    for (unsigned local_index=0; local_index < mNumberOwned; local_index++)
    {
        double v = p_solution[local_index*problemDim];
        float& activation_time = mActivationTime[local_index];
        float& peak = mPeakVoltage[local_index];
        float& apd = mActionPotentialDuration[local_index];

        if (!mActivationState[local_index] && v > mThresholdVoltage) {//activation
            mActivationState[local_index] = true;
            activation_time = (float)time;
            peak = (float)v; //reset peak voltage
        }
        if (mActivationState[local_index]) {
            //update peak
            if (v > peak)
                peak = (float)v;
            // APD90, deactivation
            if (v < peak - (peak - mRestingVoltage) * 0.9) {
                mActivationState[local_index] = false;
                apd = (float)(time - activation_time);
            }
        }
    }
    VecRestoreArray(solution, &p_solution);
}

void ActivationMapOutputModifier::ProcessPdeSolutionAtTimeStep(double time, Vec solution, unsigned problemDim) {
    ProcessSolutionAtTimeStep(time, solution, problemDim);
}
