<<<<<<< HEAD
﻿#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <Rcpp.h>

class MS2 {
public:
	MS2();
	MS2(std::vector<double>&& mz, std::vector<double>&& intensity, std::string&& info, int ppm);
public:
	// 友元函数
	friend std::pair<double, double> CompareMS2(MS2& sample_ms2, MS2& library_ms2);
public:
	std::vector<double>& GetMz();
	std::vector<double>& GetIntensity();
	std::vector<double>& GetMinMz();
	std::vector<double>& GetMaxMz();
	std::string& GetInfo();
	double GetTotalIntensity();
private:
	std::vector<double> m_mz;
	std::vector<double> m_intensity;
	std::vector<double> m_min_mz;
	std::vector<double> m_max_mz;
	std::string m_info;
	double m_total_intensity;
};

class CombineSpectrum {
public:
	CombineSpectrum();
	CombineSpectrum(double mz, double intensity, double rt, MS2&& ms2, int ppm, std::string&& index);
public:
	double GetMz();
	double GetMinMz();
	double GetMaxMz();
	double GetIntensity();
	double GetRt();
	MS2& GetMS2();
	std::string& GetIndex();
protected:
	double m_mz;
	double m_min_mz;
	double m_max_mz;
	double m_intensity;
	double m_rt;
	MS2 m_ms2;
	std::string m_index;
};

class CompareResult {
public:
	CompareResult();
public:
	void AddResult(CombineSpectrum& sample, CombineSpectrum& library, double score_1, double score_2);
	void Print();
	void OutputCsv(std::string output_path);
	int GetSize();
private:
	std::vector<CombineSpectrum*> m_sample;
	std::vector<CombineSpectrum*> m_library;
	std::vector<double> m_score_1;
	std::vector<double> m_score_2;
};

class Comparator {
public:
	Comparator();
public:
	void Compare();
public:
	void Test();
	int GetSampleSize();
	int GetLibrarySize();
public:
	// 需要按照mz排序好
	std::vector<CombineSpectrum> m_sample;
	std::vector<CombineSpectrum> m_library;
	CompareResult m_compare_res;
};

// 比较两个二级
std::pair<double, double> CompareMS2(MS2& sample_ms2, MS2& library_ms2);

// 加载数据
void Load(Comparator* m_comparator, Rcpp::List& sample_List, std::vector<double>& library_mz, double ms1_ppm, double ms2_ppm, Rcpp::List& library_list);

// 输出csv文件
=======
﻿#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <Rcpp.h>

class MS2 {
public:
	MS2();
	MS2(std::vector<double>&& mz, std::vector<double>&& intensity, std::string&& info, int ppm);
public:
	// 友元函数
	friend std::pair<double, double> CompareMS2(MS2& sample_ms2, MS2& library_ms2);
public:
	std::vector<double>& GetMz();
	std::vector<double>& GetIntensity();
	std::vector<double>& GetMinMz();
	std::vector<double>& GetMaxMz();
	std::string& GetInfo();
	double GetTotalIntensity();
private:
	std::vector<double> m_mz;
	std::vector<double> m_intensity;
	std::vector<double> m_min_mz;
	std::vector<double> m_max_mz;
	std::string m_info;
	double m_total_intensity;
};

class CombineSpectrum {
public:
	CombineSpectrum();
	CombineSpectrum(double mz, double intensity, double rt, MS2&& ms2, int ppm, std::string&& index);
public:
	double GetMz();
	double GetMinMz();
	double GetMaxMz();
	double GetIntensity();
	double GetRt();
	MS2& GetMS2();
	std::string& GetIndex();
protected:
	double m_mz;
	double m_min_mz;
	double m_max_mz;
	double m_intensity;
	double m_rt;
	MS2 m_ms2;
	std::string m_index;
};

class CompareResult {
public:
	CompareResult();
public:
	void AddResult(CombineSpectrum& sample, CombineSpectrum& library, double score_1, double score_2);
	void Print();
	void OutputCsv(std::string output_path);
	int GetSize();
private:
	std::vector<CombineSpectrum*> m_sample;
	std::vector<CombineSpectrum*> m_library;
	std::vector<double> m_score_1;
	std::vector<double> m_score_2;
};

class Comparator {
public:
	Comparator();
public:
	void Compare();
public:
	void Test();
	int GetSampleSize();
	int GetLibrarySize();
public:
	// 需要按照mz排序好
	std::vector<CombineSpectrum> m_sample;
	std::vector<CombineSpectrum> m_library;
	CompareResult m_compare_res;
};

// 比较两个二级
std::pair<double, double> CompareMS2(MS2& sample_ms2, MS2& library_ms2);

// 加载数据
void Load(Comparator* m_comparator, Rcpp::List& sample_List, std::vector<double>& library_mz, double ms1_ppm, double ms2_ppm, Rcpp::List& library_list);

// 输出csv文件
>>>>>>> master
void OutputCsv(Comparator* m_comparator, std::string output_path);