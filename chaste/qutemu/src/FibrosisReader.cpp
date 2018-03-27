#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include "FibrosisReader.hpp"

//#include <sstream>
#include "Exception.hpp"

template<typename T> bool TryParse(std::string s, T& out)
{
    try {
        out = boost::lexical_cast<T>(s);
        return true;
    } catch (boost::bad_lexical_cast e) {
        return false;
    }
}

FibrosisReader::FibrosisReader(const FileFinder& rFileFinder)
{
    mFilePath = rFileFinder.GetAbsolutePath();
    mDataFile.open(mFilePath.c_str());
    if (!mDataFile.is_open())
        EXCEPTION("Failed to open fibre file " + rFileFinder.GetAbsolutePath());

    mCurLine = 0;
    // Note: this method will close the file on error
    ReadHeader();
}

FibrosisReader::~FibrosisReader()
{
    mDataFile.close();
}

void FibrosisReader::GetAll(std::vector<double>& fibrosis)
{
    assert(fibrosis.empty());
    fibrosis.reserve(mNumLinesOfData);
    for (unsigned i=0; i<mNumLinesOfData; i++)
    {
        std::string line = NextLine();
        double item;
        if (!TryParse(line, item))
            EXCEPTION("Invalid Line(" << mCurLine << "): " << line);

        fibrosis.push_back(item);
    }
}

std::string FibrosisReader::NextLine()
{
    std::string line;
    while (true)
    {
        getline(mDataFile, line);
        mCurLine++;
        if (line.empty() && mDataFile.eof())
        {
            mDataFile.close();
            EXCEPTION("End of file " + mFilePath);
        }

        line = line.substr(0, line.find('#', 0));
        boost::algorithm::trim(line);
        if (line.length() > 0)
            return line;
    }
}

void FibrosisReader::ReadHeader()
{
    try
    {
        std::string line = NextLine();
        if (!TryParse(line, mNumLinesOfData))
            EXCEPTION("Invalid Header Line: " << line);
    }
    catch (Exception e){
        mDataFile.close();
        EXCEPTION(e.GetMessage());
    }
}