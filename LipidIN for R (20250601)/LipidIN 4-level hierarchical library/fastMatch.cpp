#include "fastMatch.h"
#include <Rcpp.h>

using namespace std;

fastMatch::fastMatch()
{
}

fastMatch::fastMatch(std::vector<ms1>& ms1_vecotr, std::vector<ms2>& ms2_vector, std::vector<double>& ms1_database_mz, double ppm, string output_path)
{
    this->m_ms1_vecotr = ms1_vecotr;
    this->m_ms2_vector = ms2_vector;
    for (double num : ms1_database_mz) {
        this->m_library_min_mz.push_back(num - num * ppm / 1000000);
        this->m_library_max_mz.push_back(num + num * ppm / 1000000);
    }
    this->m_output_path = output_path;
}

void fastMatch::compare()
{
    // 输出的csv文件
    ofstream dataFile(this->m_output_path);
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
    // 进行查找
    for (auto itr = this->m_ms1_vecotr.begin(); itr != this->m_ms1_vecotr.end(); itr++) {
        itr->compare(this->m_library_min_mz, this->m_library_max_mz, this->m_ms2_vector, dataFile);
    }
    dataFile.close();
}

std::vector<ms1> fastMatch::GetMs1Vector()
{
    return this->m_ms1_vecotr;
}

std::vector<ms2> fastMatch::GetMs2Vector()
{
    return this->m_ms2_vector;
}

std::vector<double> fastMatch::GetLibraryMinMz()
{
    return this->m_library_min_mz;
}

std::vector<double> fastMatch::GetLibraryMaxMz()
{
    return this->m_library_max_mz;
}

string fastMatch::GetOutputPath()
{
    return this->m_output_path;
}

ms1::ms1()
{
}

ms1::ms1(double precursor_mz, double precursor_intensity, double rt, ms2 ms2, string index)
{
    this->m_precursor_mz = precursor_mz;
    this->m_precursor_intensity = precursor_intensity;
    this->m_rt = rt;
    this->m_ms2 = move(ms2);
    this->m_index = move(index);
}

void ms1::compare(std::vector<double>& min_mz, std::vector<double>& max_mz, std::vector<ms2>& ms2_vector, std::ofstream& data_file)
{
    // 二分查找,min_mz需要已经排序好
    int left = 0;
    int right = min_mz.size() - 1;

    while (left <= right) {
        int mid = (left + right) / 2;
        if (min_mz.at(mid) <= this->m_precursor_mz && this->m_precursor_mz <= max_mz.at(mid)) {
            // 比较两个二级
            pair<double, double> compare_ms2_res = this->m_ms2.compare(ms2_vector.at(mid));
            if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
                // 输出信息
                data_file << this->m_index << "," << this->m_precursor_mz << "," << this->m_rt << "," << this->m_precursor_intensity << "," << ms2_vector.at(mid).GetInfo() << "," << compare_ms2_res.first << "," << compare_ms2_res.second << endl;
            }
            // 向左和右继续寻找
            int left_t = mid - 1;
            int right_t = mid + 1;

            while (left_t >= 0) {
                if (min_mz.at(left_t) <= this->m_precursor_mz && this->m_precursor_mz <= max_mz.at(left_t)) {
                    // 比较两个二级
                    compare_ms2_res = this->m_ms2.compare(ms2_vector.at(left_t));
                    if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
                        // 输出信息
                        data_file << this->m_index << "," << this->m_precursor_mz << "," << this->m_rt << "," << this->m_precursor_intensity << "," << ms2_vector.at(left_t).GetInfo() << "," << compare_ms2_res.first << "," << compare_ms2_res.second << endl;
                    }
                    // 向左继续找
                    left_t--;
                } else {
                    break;
                }
            }
            while (right_t <= int(min_mz.size() - 1)) {
                if (min_mz.at(right_t) <= this->m_precursor_mz && this->m_precursor_mz <= max_mz.at(right_t)) {
                    // 比较两个二级
                    compare_ms2_res = this->m_ms2.compare(ms2_vector.at(right_t));
                    if (compare_ms2_res.first != 0 && compare_ms2_res.second != 0) {
                        // 输出信息
                        data_file << this->m_index << "," << this->m_precursor_mz << "," << this->m_rt << "," << this->m_precursor_intensity << "," << ms2_vector.at(right_t).GetInfo() << "," << compare_ms2_res.first << "," << compare_ms2_res.second << endl;
                    }
                    // 向左继续找
                    right_t++;
                } else {
                    break;
                }
            }
            break; // 找完左边和右边的，就说明ms1_itr的元素已经找完了，跳出这个while
        } else if (this->m_precursor_mz <= min_mz.at(mid)) {
            right = mid - 1;
        } else if (this->m_precursor_mz >= max_mz.at(mid)) {
            left = mid + 1;
        }
    }
}

double ms1::GetPrecursorMz()
{
    return this->m_precursor_mz;
}

double ms1::GetPrecursorIntensity()
{
    return this->m_precursor_intensity;
}

double ms1::GetRt()
{
    return this->m_rt;
}

ms2 ms1::GetMs2()
{
    return this->m_ms2;
}

ms2::ms2()
{
}

ms2::ms2(std::vector<double> mz, std::vector<double> intensity, std::vector<double> min_mz, std::vector<double> max_mz, std::string info, bool count_total_intensity)
{
    this->m_mz = move(mz);
    this->m_intensity = move(intensity);
    this->m_min_mz = move(min_mz);
    this->m_max_mz = move(max_mz);
    this->m_info = move(info);

    if (count_total_intensity) {
        // 计算总强度
        double total_intensity = 0;
        for (double num : m_intensity) {
            total_intensity += num;
        }
        this->m_total_intensity = total_intensity;
    } else {
        this->m_total_intensity = 0;
    }
}

std::pair<double, double> const ms2::compare(ms2& library_ms2)
{
    // 二分查找,library_ms2的min_mz和max_mz需要已经排序好
    vector<double>& library_min_mz = library_ms2.GetMinMz();
    vector<double>& library_max_mz = library_ms2.GetMaxMz();
    vector<double>& library_mz = library_ms2.GetMz();

    // 用于去重
    unordered_map<double, int> m_hash_map;

    // first存储score1，second存储score2
    pair<double, double> res;

    for (auto itr = this->m_mz.begin(); itr != this->m_mz.end(); itr++) {
        int left = 0;
        int right = library_min_mz.size() - 1;

        while (left <= right) {
            int mid = (left + right) / 2;
            if (library_min_mz.at(mid) <= *(itr) && *(itr) <= library_max_mz.at(mid)) {
                auto find_itr = m_hash_map.find(library_mz.at(mid)); // 寻找这个理论库的这个碎片是否被匹配过
                // 如果未被匹配过,将理论库的mz作为键，样本匹配上的碎片的位置存入哈希表
                if (find_itr == m_hash_map.end()) {
                    m_hash_map[library_mz.at(mid)] = itr - this->m_mz.begin();
                } else {
                    // 计算上一个匹配的diff和新匹配上的diff
                    double origin_dif = abs(find_itr->first - this->m_mz.at(find_itr->second));
                    double new_diff = abs(find_itr->first - (*itr));
                    // 如果新的diff比上一次匹配的diff小，则更新哈希表对应的值
                    if (new_diff <= origin_dif) {
                        m_hash_map[find_itr->first] = itr - this->m_mz.begin();
                    }
                }
                break;
            } else if (*(itr) <= library_min_mz.at(mid)) {
                right = mid - 1;
            } else if (*(itr) >= library_max_mz.at(mid)) {
                left = mid + 1;
            }
        }
    }

    double score1 = 0;
    double score2 = 0;

    // 计算score1
    for (auto itr = m_hash_map.begin(); itr != m_hash_map.end(); itr++) {
        score2 += this->m_intensity.at(itr->second);
    }
    score2 = score2 / this->m_total_intensity;

    // 计算score2
    score1 = double(m_hash_map.size()) / library_mz.size();

    res.first = score1;
    res.second = score2;
    return res;
}

std::vector<double>& ms2::GetMz()
{
    return this->m_mz;
}

std::vector<double>& ms2::GetIntneisty()
{
    return this->m_intensity;
}

std::vector<double>& ms2::GetMinMz()
{
    return this->m_min_mz;
}

std::vector<double>& ms2::GetMaxMz()
{
    return this->m_max_mz;
}

string ms2::GetInfo()
{
    return this->m_info;
}

// [[Rcpp::export]]
void SuperFastCompare(Rcpp::List& sample_List, std::vector<double>& library_mz, double ppm, Rcpp::List& library_list, std::string output_path)
{
    vector<ms1> ms1_vector;
    vector<ms2> ms2_vector;
    clock_t start, end;

    start = clock();
    // 载入一级向量
    for (int i = 0; i < sample_List.size(); i++) {
        Rcpp::List single_List = sample_List.at(i);
        double precursor_mz = single_List["PrecursorMZ"];
        double precursor_intensity = single_List["PrecursorIntensity"];
        double rt = single_List["rt"];
        string index = single_List["num"];
        // 获取二级那个数据框
        Rcpp::List ms2_List = single_List["MS2mz"];
        vector<double> ms2_mz = ms2_List["da.temp.mz"];
        vector<double> ms2_intensity = ms2_List["da.temp.intensity"];

        vector<double> null_vector = {};

        // 二级
        ms2 this_ms2 = ms2(move(ms2_mz), move(ms2_intensity), null_vector, null_vector, "", true);

        // 加入ms1_vector中
        ms1_vector.emplace_back(ms1(precursor_mz, precursor_intensity, rt, move(this_ms2), move(index)));
    }

    // 载入二级向量
    for (int i = 0; i < library_list.size(); i++) {
        Rcpp::List single_List = library_list.at(i);
        string info = single_List["inf"];

        // 获取二级那个数据框
        Rcpp::List ms2_List = single_List["ms2.spe"];

        vector<double> ms2_mz = ms2_List["mz"];
        vector<double> ms2_intensity = ms2_List["intensity"];
        vector<double> ms2_min_mz = ms2_List["down"];
        vector<double> ms2_max_mz = ms2_List["up"];

        // 加入ms2_vector中
        ms2_vector.emplace_back(ms2(move(ms2_mz), move(ms2_intensity), move(ms2_min_mz), move(ms2_max_mz), info, false));
    }
    end = clock();
    Rcpp::Rcout << "copying cost:" << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();

    fastMatch* process = new fastMatch(ms1_vector, ms2_vector, library_mz, ppm, output_path);
    process->compare();

    end = clock();
    Rcpp::Rcout << "comparing cost:" << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
    return;
}
