#include "CsvReader.h"

#include <fstream>
#include <sstream>

void CsvReader::read(const string& file)
{

    rows = 0;
    columns = 0;
    data.clear();

    ifstream infile;
    infile.open(file);

    string line;

    while (getline(infile, line)) {
        vector<double> row;
        stringstream tokenizer(line);
        string element;
        while (getline(tokenizer, element, ',')) {
            row.push_back(stod(element));
        }
        data.push_back(row);
        if (columns == 0)
            columns = row.size();
        rows++;
    }
    infile.close();

}

double CsvReader::get(int r, int c) const
{
    return data[r][c];
}

unsigned long CsvReader::getRows() const
{
    return rows;
}

unsigned long CsvReader::getColumns() const
{
    return columns;
}

CsvReader::CsvReader() = default;

CsvReader::CsvReader(const string& file)
{
    read(file);
}

void CsvReader::print() const
{
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
            std::printf("%+.4f ", data[r][c]);
        }
        std::printf("\n");
    }
}
