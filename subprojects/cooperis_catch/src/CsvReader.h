#pragma once

#include <string>
#include <vector>

using namespace std;

class CsvReader {

public:
    CsvReader();
    CsvReader(const string& file);

    void read(const string& file);

    double get(int r, int c) const;
    unsigned long getRows() const;
    unsigned long getColumns() const;

    void print() const;

protected:
    vector<vector<double>> data;
    unsigned long rows = 0;
    unsigned long columns = 0;
};
