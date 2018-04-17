#pragma once

#include <vector>
#include <hdf5.h>

#include "FileFinder.hpp"

class ConductivityReader
{
private:
    hid_t mFileId;

public:
    static std::vector<float> ReadConductivities(const FileFinder& rFileFinder);
};
