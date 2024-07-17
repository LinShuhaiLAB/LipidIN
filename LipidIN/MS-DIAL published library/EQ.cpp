#include "SpeedKing.h"

using namespace std;

MS2::MS2()
{
	this->m_mz = {};
	this->m_intensity = {};
	this->m_info = "";
	this->m_min_mz = {};
	this->m_max_mz = {};
	this->m_total_intensity = 0;
}

MS2::MS2(std::vector<double>&& mz, std::vector<double>&& intensity, std::string&& info, int ppm)
{
	this->m_mz = move(mz);
	this->m_intensity = move(intensity);
	this->m_info = info;


	for (auto i : this->m_mz) {
		this->m_min_mz.push_back(i - i * ppm / 1000000);
		this->m_max_mz.push_back(i + i * ppm / 1000000);
	}


	double total_intensity = 0;
	for (double num : m_intensity) {
		total_intensity += num;
	}
	this->m_total_intensity = total_intensity;
}

std::vector<double>& MS2::GetMz()
{
	return this->m_mz;
}

std::vector<double>& MS2::GetIntensity()
{
	return this->m_intensity;
}

std::vector<double>& MS2::GetMinMz()
{
	return this->m_min_mz;
}

std::vector<double>& MS2::GetMaxMz()
{
	return this->m_max_mz;
}

std::string& MS2::GetInfo()
{
	return this->m_info;
}

double MS2::GetTotalIntensity()
{
	return this->m_total_intensity;
}

CombineSpectrum::CombineSpectrum()
{
	this->m_mz = 0;
	this->m_intensity = 0;
	this->m_rt = 0;
	this->m_ms2 = MS2();
	this->m_min_mz = 0;
	this->m_max_mz = 0;
}

CombineSpectrum::CombineSpectrum(double mz, double intensity, double rt, MS2&& ms2, int ppm, std::string&& index)
{
	this->m_mz = mz;
	this->m_intensity = intensity;
	this->m_rt = rt;
	this->m_ms2 = move(ms2);
	this->m_index = move(index);
	// 计算min_mz和max_mz
	this->m_min_mz = this->m_mz - this->m_mz * ppm / 1000000;
	this->m_max_mz = this->m_mz + this->m_mz * ppm / 1000000;
}

double CombineSpectrum::GetMz()
{
	return this->m_mz;
}

double CombineSpectrum::GetMinMz()
{
	return this->m_min_mz;
}

double CombineSpectrum::GetMaxMz()
{
	return this->m_max_mz;
}

double CombineSpectrum::GetIntensity()
{
	return this->m_intensity;
}

double CombineSpectrum::GetRt()
{
	return this->m_rt;
}

MS2& CombineSpectrum::GetMS2()
{
	return this->m_ms2;
}

std::string& CombineSpectrum::GetIndex()
{
	return this->m_index;
}

Comparator::Comparator()
{
	Rcpp::Rcout << "####################" << endl;
	Rcpp::Rcout << "#   ⚡⚡⚡⚡⚡   #" << endl;
	Rcpp::Rcout << "#⚡LipidIN  EQ #" << endl;
	Rcpp::Rcout << "#   ⚡⚡⚡⚡⚡   #" << endl;
	Rcpp::Rcout << "####################" << endl;
}



void Comparator::Compare()
{
	Rcpp::Rcout << "Sample size:" << this->GetSampleSize() << endl;
	Rcpp::Rcout << "Library size:" << this->GetLibrarySize() << endl;

	clock_t start, end;
	start = clock();

	int library_left = 0;
	int library_right = m_library.size() - 1;
	int sample_left = 0;
	int sample_right = m_sample.size() - 1;


	while (sample_left <= sample_right)
	{
		int left = library_left;
		int right = library_right;
		bool find_sample_left = 0;
		bool find_sample_right = 0;

		while (left <= right) {
			int mid = (left + right) / 2;
			double sample_mz = m_sample[sample_left].GetMz();

			if (sample_mz >= m_library[mid].GetMinMz() && sample_mz <= m_library[mid].GetMaxMz()) {

				find_sample_left = 1;

				auto compare_ms2_res = CompareMS2(m_sample[sample_left].GetMS2(), m_library[mid].GetMS2());

				if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
					m_compare_res.AddResult(m_sample[sample_left], m_library[mid], compare_ms2_res.first, compare_ms2_res.second);
				}

				int left_t = mid - 1;
				int right_t = mid + 1;
				while (left_t >= 0) {
					if (m_library[left_t].GetMinMz() <= sample_mz && sample_mz <= m_library[left_t].GetMaxMz()) {

						compare_ms2_res = CompareMS2(m_sample[sample_left].GetMS2(), m_library[left_t].GetMS2());

						if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
							m_compare_res.AddResult(m_sample[sample_left], m_library[left_t], compare_ms2_res.first, compare_ms2_res.second);
						}

						left_t--;
					}
					else {
						break;
					}
				}
				while (right_t <= int(m_library.size() - 1)) {
					if (m_library[right_t].GetMinMz() <= sample_mz && sample_mz <= m_library[right_t].GetMaxMz()) {

						compare_ms2_res = CompareMS2(m_sample[sample_left].GetMS2(), m_library[right_t].GetMS2());

						if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
							m_compare_res.AddResult(m_sample[sample_left], m_library[right_t], compare_ms2_res.first, compare_ms2_res.second);
						}

						right_t++;
					}
					else {
						break;
					}
				}
	
				left = left_t + 1;
				library_left = left_t + 1;
				break; 
			}
			else if (sample_mz <= m_library[mid].GetMinMz()) {
				right = mid - 1;
			}
			else if (sample_mz >= m_library[mid].GetMaxMz()) {
				left = mid + 1;
			}
		}


		if (sample_left == sample_right) {
			break;
		}

		else if (left > int(m_library.size() - 1)) {
			break;
		}


		if (!find_sample_left) {
			library_left = left;
		}


		right = library_right;



		while (left <= right) {
			int mid = (left + right) / 2;
			double sample_mz = m_sample[sample_right].GetMz();

			if (sample_mz >= m_library[mid].GetMinMz() && sample_mz <= m_library[mid].GetMaxMz()) {

				find_sample_right = 1;

				auto compare_ms2_res = CompareMS2(m_sample[sample_right].GetMS2(), m_library[mid].GetMS2());

				if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
					m_compare_res.AddResult(m_sample[sample_right], m_library[mid], compare_ms2_res.first, compare_ms2_res.second);
				}

				int left_t = mid - 1;
				int right_t = mid + 1;
				while (left_t >= 0) {
					if (m_library[left_t].GetMinMz() <= sample_mz && sample_mz <= m_library[left_t].GetMaxMz()) {

						compare_ms2_res = CompareMS2(m_sample[sample_right].GetMS2(), m_library[left_t].GetMS2());

						if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
							m_compare_res.AddResult(m_sample[sample_right], m_library[left_t], compare_ms2_res.first, compare_ms2_res.second);
						}

						left_t--;
					}
					else {
						break;
					}
				}

				while (right_t <= int(m_library.size() - 1)) {
					if (m_library[right_t].GetMinMz() <= sample_mz && sample_mz <= m_library[right_t].GetMaxMz()) {

						compare_ms2_res = CompareMS2(m_sample[sample_right].GetMS2(), m_library[right_t].GetMS2());

						if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
							m_compare_res.AddResult(m_sample[sample_right], m_library[right_t], compare_ms2_res.first, compare_ms2_res.second);
						}

						right_t++;
					}
					else {
						break;
					}
				}

				library_right = right_t - 1;
				break; 
			}
			else if (sample_mz <= m_library[mid].GetMinMz()) {
				right = mid - 1;
			}
			else if (sample_mz >= m_library[mid].GetMaxMz()) {
				left = mid + 1;
			}
		}


		if (!find_sample_right) {
			library_right = right;
		}


		sample_left++;
		sample_right--;
	}


	end = clock();
	Rcpp::Rcout << "Comparing cost:" << (double)(end - start) / CLOCKS_PER_SEC << "s ⚡" << endl;

	Rcpp::Rcout << "Comparing result get:" << m_compare_res.GetSize() << endl;
}



void Comparator::Test()
{
	cout << "Runing Compartor testing..." << endl;

	MS2 ms2_1({ 100,200 }, { 999,50 }, "test1", 30);
	MS2 ms2_2({ 300,400 }, { 999,50 }, "test2", 30);
	MS2 ms2_3({ 500,600 }, { 999,50 }, "test3", 30);

	MS2 ms2_4({ 100,200 }, { 999,50 }, "test1", 30);
	MS2 ms2_5({ 300,400 }, { 999,50 }, "test2", 30);
	MS2 ms2_6({ 500,600 }, { 999,50 }, "test3", 30);


	CombineSpectrum sample_1(100, 100, 100, move(ms2_1), 5, "100");
	CombineSpectrum sample_2(200, 200, 200, move(ms2_2), 5, "200");
	CombineSpectrum sample_3(300, 300, 300, move(ms2_3), 5, "300");

	m_sample.push_back(sample_1);
	m_sample.push_back(sample_2);
	m_sample.push_back(sample_3);


	CombineSpectrum library_1(100, 100, 100, move(ms2_4), 5, "100");
	CombineSpectrum library_2(200, 200, 200, move(ms2_5), 5, "200");
	CombineSpectrum library_3(300, 300, 300, move(ms2_6), 5, "300");

	m_library.push_back(library_1);
	m_library.push_back(library_2);
	m_library.push_back(library_3);


	this->Compare();

	this->m_compare_res.Print();

	cout << "Testing finish" << endl;
}

int Comparator::GetSampleSize()
{
	return this->m_sample.size();
}

int Comparator::GetLibrarySize()
{
	return this->m_library.size();
}

std::pair<double, double> CompareMS2(MS2& sample_ms2, MS2& library_ms2)
{

	vector<double>& library_min_mz = library_ms2.GetMinMz();
	vector<double>& library_max_mz = library_ms2.GetMaxMz();
	vector<double>& library_mz = library_ms2.GetMz();


	unordered_map<double, int> m_hash_map;


	pair<double, double> res;


	for (auto itr = sample_ms2.m_mz.begin(); itr != sample_ms2.m_mz.end(); itr++) {
		int left = 0;
		int right = library_min_mz.size() - 1;

		while (left <= right) {
			int mid = (left + right) / 2;
			if (library_min_mz.at(mid) <= *(itr) && *(itr) <= library_max_mz.at(mid)) {

				auto find_itr = m_hash_map.find(library_mz.at(mid));

				if (find_itr == m_hash_map.end()) {
					m_hash_map[library_mz.at(mid)] = itr - sample_ms2.m_mz.begin();
				}
				else {

					double origin_dif = abs(find_itr->first - sample_ms2.m_mz.at(find_itr->second));
					double new_diff = abs(find_itr->first - (*itr));

					if (new_diff <= origin_dif) {
						m_hash_map[find_itr->first] = itr - sample_ms2.m_mz.begin();
					}
				}
				break;
			}
			else if (*(itr) <= library_min_mz.at(mid)) {
				right = mid - 1;
			}
			else if (*(itr) >= library_max_mz.at(mid)) {
				left = mid + 1;
			}
		}
	}

	double score1 = 0;
	double score2 = 0;


	for (auto itr = m_hash_map.begin(); itr != m_hash_map.end(); itr++) {
		score2 += sample_ms2.m_intensity.at(itr->second);
	}
	score2 = score2 / sample_ms2.m_total_intensity;


	score1 = double(m_hash_map.size()) / library_mz.size();

	res.first = score1;
	res.second = score2;
	return res;
}

void Load(Comparator* m_comparator, Rcpp::List& sample_List, std::vector<double>& library_mz, double ms1_ppm, double ms2_ppm, Rcpp::List& library_list)
{

	m_comparator->m_sample.reserve(sample_List.size());
	m_comparator->m_library.reserve(library_list.size());

	clock_t start, end;
	start = clock();

	for (int i = 0; i < sample_List.size(); i++) {
		Rcpp::List single_List = sample_List.at(i);
		double precursor_mz = single_List["PrecursorMZ"];
		double precursor_intensity = single_List["PrecursorIntensity"];
		double rt = single_List["rt"];
		string index = single_List["num"];

		Rcpp::List ms2_List = single_List["MS2mz"];
		vector<double> ms2_mz = ms2_List["da.temp.mz"];
		vector<double> ms2_intensity = ms2_List["da.temp.intensity"];


		MS2 this_ms2(move(ms2_mz), move(ms2_intensity), "", ms2_ppm);


		m_comparator->m_sample.emplace_back(CombineSpectrum(precursor_mz, precursor_intensity, rt, move(this_ms2), ms1_ppm, move(index)));
	}
	end = clock();
	Rcpp::Rcout << "ms1 sample copying cost:" << ((double)(end - start) / CLOCKS_PER_SEC) << "s" << endl;


	start = clock();

	for (int i = 0; i < library_mz.size(); i++) {
		Rcpp::List single_List = library_list.at(i);
		string info = single_List["inf"];


		Rcpp::List ms2_List = single_List["ms2.spe"];

		vector<double> ms2_mz = ms2_List["mz"];
		vector<double> ms2_intensity = ms2_List["intensity"];


		MS2 this_ms2(move(ms2_mz), move(ms2_intensity), move(info), ms2_ppm);


		m_comparator->m_library.emplace_back(CombineSpectrum(library_mz[i], 0, 0, move(this_ms2), ms1_ppm, ""));
	}
	end = clock();
	Rcpp::Rcout << "ms2 library copying cost:" << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void Load2(Comparator* m_comparator, Rcpp::List& sample_List, double ms1_ppm, double ms2_ppm)
{

  m_comparator->m_sample.reserve(sample_List.size());
  
  clock_t start, end;
  start = clock();

  for (int i = 0; i < sample_List.size(); i++) {
    Rcpp::List single_List = sample_List.at(i);
    double precursor_mz = single_List["PrecursorMZ"];
    double precursor_intensity = single_List["PrecursorIntensity"];
    double rt = single_List["rt"];
    string index = single_List["num"];

    Rcpp::List ms2_List = single_List["MS2mz"];
    vector<double> ms2_mz = ms2_List["da.temp.mz"];
    vector<double> ms2_intensity = ms2_List["da.temp.intensity"];
    

    MS2 this_ms2(move(ms2_mz), move(ms2_intensity), "", ms2_ppm);
    

    m_comparator->m_sample.emplace_back(CombineSpectrum(precursor_mz, precursor_intensity, rt, move(this_ms2), ms1_ppm, move(index)));
  }
  end = clock();
  Rcpp::Rcout << "ms1 sample copying cost:" << ((double)(end - start) / CLOCKS_PER_SEC) << "s" << endl;
  
}

void OutputCsv(Comparator* m_comparator, std::string output_path)
{
	m_comparator->m_compare_res.OutputCsv(output_path);
}


RCPP_MODULE(unif_module) {
	using namespace Rcpp;
	class_<Comparator>("Comparator")
		.constructor()
		// .field("min", &Uniform::min)
		// .field("max", &Uniform::max)
		//.method("draw", &Uniform::draw)
		.method("Load", &Load)
    .method("Load2", &Load2)
		.method("OutputCsv", &OutputCsv)
		.method("Compare", &Comparator::Compare)
		;
}

CompareResult::CompareResult()
{
}

void CompareResult::AddResult(CombineSpectrum& sample, CombineSpectrum& library, double score_1, double score_2)
{
	this->m_sample.push_back(&sample);
	this->m_library.push_back(&library);
	this->m_score_1.push_back(score_1);
	this->m_score_2.push_back(score_2);
}

void CompareResult::Print()
{
	string res;
	for (size_t i = 0; i < m_sample.size(); i++) {
		res += "SampleMz:" + to_string(m_sample[i]->GetMz()) + " " + "SampleRt:" + to_string(m_sample[i]->GetRt());
		res += "   ";
		res += "LibraryMz:" + to_string(m_library[i]->GetMz()) + " " + "LibraryRt:" + to_string(m_library[i]->GetRt());
		res += "   ";
		res += "score_1:" + to_string(m_score_1[i]) + " ";
		res += "score_2:" + to_string(m_score_2[i]);
		res += "\n";
	}
	cout << res << endl;
}

void CompareResult::OutputCsv(std::string output_path)
{

	ofstream dataFile(output_path);

	dataFile << "index"
		<< ","
		<< "mz"
		<< ","
		<< "rt"
		<< ","
		<< "intensity"
		<< ","
		<< "title"
		<< ","
		<< "score1"
		<< ","
		<< "score2" << endl;


	for (int i = 0; i < m_sample.size(); i++) {
		dataFile << m_sample[i]->GetIndex() << "," << m_sample[i]->GetMz() << "," << m_sample[i]->GetRt() << "," << m_sample[i]->GetIntensity() << ",";
		dataFile << m_library[i]->GetMS2().GetInfo() << ",";
		dataFile << m_score_1[i] << ",";
		dataFile << m_score_2[i] << endl;
	}
	dataFile.close();
	m_sample.clear();m_library.clear();m_score_1.clear();m_score_2.clear();
}

int CompareResult::GetSize()
{
	return this->m_sample.size();
}

