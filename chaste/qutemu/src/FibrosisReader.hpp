#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "UblasIncludes.hpp"
#include "FileFinder.hpp"

class FibrosisReader
{
private:
    /** File stream to use for GetTokensAtNextLine */
    std::ifstream mDataFile;

    /** Absolute path of the file being read */
    std::string mFilePath;

    /** Number of lines of data in the file, read from the first line of the file */
    unsigned mNumLinesOfData;
    unsigned mCurLine;

    void ReadHeader();
    std::string NextLine();
    unsigned ReadNextItemLine(double& item);

public:
    /**
     * Create a new FibreReader.
     *
     * @param rFileFinder  the path to the fibre direction file
     */
    FibrosisReader(const FileFinder& rFileFinder);

    /**
     *  Destructor closes file.
     */
    ~FibrosisReader();

    unsigned GetNumLinesOfData()
    {
        return mNumLinesOfData;
    }

    void GetAll(std::vector<double>& fibrosis);
};
