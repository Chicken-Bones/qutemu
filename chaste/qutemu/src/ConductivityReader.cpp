#include <H5Fpublic.h>
#include <string>

#include "ConductivityReader.hpp"
#include "Exception.hpp"

std::vector<float> ConductivityReader::ReadConductivities(const FileFinder &h5_file) {
    std::string file_name = h5_file.GetAbsolutePath();
    if (!h5_file.Exists())
        EXCEPTION("Could not open " << file_name << " , as it does not exist.");

    hid_t mFileId = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (mFileId <= 0)
        EXCEPTION("Could not open " << file_name << " , H5Fopen error code = " << mFileId);

    hid_t datasetId = H5Dopen(mFileId, "Conductivity", H5P_DEFAULT);
    if (datasetId <= 0)
    {
        H5Fclose(mFileId);
        EXCEPTION("Opened " << file_name << " but could not find the dataset 'Conductivity', H5Dopen error code = " << datasetId);
    }

    hid_t dspace = H5Dget_space(datasetId);
    std::vector<float> data(H5Sget_simple_extent_npoints(dspace));
    H5Sclose(dspace);

    H5Dread(datasetId, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    H5Dclose(datasetId);
    H5Fclose(mFileId);
    return data;
}
