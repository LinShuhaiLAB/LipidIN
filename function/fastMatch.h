#ifndef DATA_H
#define DATA_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <unordered_map>
#include <utility>
#include <vector>
class ms2 {
public:
    ms2();
    ms2(std::vector<double> mz, std::vector<double> intensity, std::vector<double> min_mz, std::vector<double> max_mz, std::string info, bool count_total_intensity = 1);

public:
    std::pair<double, double> const compare(ms2& library_ms2);

public:
    std::vector<double>& GetMz();
    std::vector<double>& GetIntneisty();
    std::vector<double>& GetMinMz();
    std::vector<double>& GetMaxMz();
    std::string GetInfo();

private:
    std::vector<double> m_mz;
    std::vector<double> m_intensity;
    std::vector<double> m_min_mz;
    std::vector<double> m_max_mz;
    std::string m_info;
    double m_total_intensity;
};

class ms1 {
public:
    ms1();
    ms1(double precursor_mz, double precursor_intensity, double rt, ms2 ms2, std::string index);

public:
    void compare(std::vector<double>& min_mz, std::vector<double>& max_mz, std::vector<ms2>& ms2_vector, std::ofstream& data_file);

public:
    double GetPrecursorMz();
    double GetPrecursorIntensity();
    double GetRt();
    ms2 GetMs2();

private:
    std::string m_index;
    double m_precursor_mz;
    double m_precursor_intensity;
    double m_rt;
    ms2 m_ms2;
};

class fastMatch {
public:
    fastMatch();
    fastMatch(std::vector<ms1>& ms1_vecotr, std::vector<ms2>& ms2_vector, std::vector<double>& ms1_database_mz, double ppm, std::string output_path);

public:
    void compare();
    std::vector<ms1> GetMs1Vector();
    std::vector<ms2> GetMs2Vector();
    std::vector<double> GetLibraryMinMz();
    std::vector<double> GetLibraryMaxMz();
    std::string GetOutputPath();

private:
    double m_ppm = 5;
    std::vector<ms1> m_ms1_vecotr;
    std::vector<ms2> m_ms2_vector;
    std::vector<double> m_library_min_mz;
    std::vector<double> m_library_max_mz;
    std::string m_output_path = "output.csv";
};

#endif // DATA_H
