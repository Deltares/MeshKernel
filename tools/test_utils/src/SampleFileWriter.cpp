#include "SampleFileWriter.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

void WriteSampleFileToAsc(const std::vector<std::vector<double>>& sampleDataMatrix,
                          double sampleDelta,
                          std::string const& filePath)
{
    std::ofstream outputFile;

    outputFile.open(filePath);
    outputFile << "ncols " << sampleDataMatrix[0].size() << std::endl;
    outputFile << "nrows " << sampleDataMatrix.size() << std::endl;
    outputFile << "xllcentre " << 0.0 << std::endl;
    outputFile << "xllcentre " << 0.0 << std::endl;
    outputFile << "cellsize " << sampleDelta << std::endl;
    outputFile << "NODATA_value " << 0.0 << std::endl;

    for (int r = 0; r < sampleDataMatrix.size(); r++) // columns
    {
        for (int c = 0; c < sampleDataMatrix[0].size(); c++) // rows
        {
            outputFile << std::fixed << std::setprecision(12) << sampleDataMatrix[r][c] << " ";
        }
        outputFile << std::endl;
    }

    outputFile.close();
}
