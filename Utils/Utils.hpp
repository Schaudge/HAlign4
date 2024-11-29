#pragma once
//读入，写回，预处理，后处理，统一化
#include "Fasta.hpp"
#include "Pseudo.hpp"
#include "Insertion.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>


#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#include <process.h>
#include <io.h>
inline void EmptySet()
{
    EmptyWorkingSet(GetCurrentProcess());
}
void getFiles_win(std::string path, std::vector<std::string>& files);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/types.h>
#include <dirent.h>
#include <malloc.h>
#include <unistd.h>
#include <sys/resource.h>
#include <pthread.h>
inline void EmptySet()
{
    malloc_trim(0);
}
void getFiles_linux(std::string path, std::vector<std::string>& filenames);
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

#undef max
#undef min
void cout_cur_time();
size_t getPeakRSS();
//测试内存峰值
inline void GetMemoryUsage()
{
    int mem = getPeakRSS() / 1024.0 / 1024.0;
    std::cout << "****process mem****" << std::endl;
    std::cout << "current pid: " << getpid() << std::endl;
    std::cout << "memory usage: " << mem << "MB" << std::endl;
}

namespace utils
{ 
    struct MAF_info
    {
        std::string path;
        int thresh1; //长度100
        int thresh2; //条数1
        int thresh3; //分数95
    };
    struct block
    {
        int name;
        size_t start;
        size_t length;
        //unsigned char* seqi;
        //std::string sign;
        //size_t all_length;
        std::vector<unsigned char> seqi;
    };
    struct PSA_ni_block
    {
        size_t start[2];
        size_t end[2];
        size_t length[2];
        bool sign;
        std::vector<unsigned char> a_seq;
        std::vector<unsigned char> b_seq;
        PSA_ni_block* next=NULL;
    };
    struct in_block
    {
        float score;
        float score_100;
        size_t start;
        size_t end;
        std::vector<size_t> name;
        std::vector<size_t> length;
        std::vector<size_t> si;
        in_block* next=NULL;
    };
    struct m_block
    {
        size_t start1;
        size_t end1;
        size_t start2;
        size_t end2;
        std::vector<std::tuple<int, int>> gap1;
        std::vector<std::tuple<int, int>> gap2;
    };
    using more_block = std::vector<m_block>;
    struct MAF_block
    {
        float score;
        int tag_num;  //记录有tag的列数，连续性
        std::vector<block> seq;
    };
    std::string remove_white_spaces(const std::string &str); //去掉空格

    unsigned char to_pseudo(char c); //预处理，char->int

    std::vector<unsigned char> to_pseudo(const std::string &str);
    std::string from_pseudo(const std::vector<unsigned char> &pseu);

    template<typename InputIterator, typename OutputIterator>
    void transform_to_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)   //预处理，char->int
    {
        std::vector<unsigned char> (*op)(const std::string &) = &to_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator, typename OutputIterator>
    void transform_from_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des) //后处理，int->char
    {
        std::string (*op)(const std::vector<unsigned char> &) = &from_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator>
    InputIterator iter_of_max(InputIterator first, InputIterator last) //找到最大的元素
    {
        auto result = first;

        for (; first != last; ++first) if (*result < *first) result = first;
        return result;
    }

    std::vector<std::vector<unsigned char>> read_to_pseudo(std::istream& is, std::string& center_name, int& II, int& center_);
    unsigned char* copy_DNA(const std::vector<unsigned char>& sequence, unsigned char* A, size_t a_begin, size_t a_end);
    void insert_and_write(std::ostream &os, std::istream &is, const std::vector<std::vector<Insertion>> &insertions); //原文件写回比对结果
    void write_to_fasta(std::ostream& os, std::istream& is, std::vector<std::vector<Insertion>>& insertions, size_t& II);
    void insert_and_write_file(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign); //写回比对结果
    int* vector_insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions);
    void write_to_str(std::string& ans, std::string& each_sequence, std::vector<Insertion>& insertions);
    void insert_and_write_fasta(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, bool TU);//写回fasta
    template<typename InputIterator>
    static void cut_and_write(std::ostream &os, InputIterator first, InputIterator last) //一条长序列分多行写入
    {
        const size_t sequence_length = std::distance(first, last);

        for (size_t i = 0; i < sequence_length; i += Fasta::max_line_length)
        {
            if (i) os << '\n';

            size_t write_length = sequence_length - i;
            if (write_length > Fasta::max_line_length) write_length = Fasta::max_line_length;

            const auto begin = first;
            std::advance(first, write_length);
            std::copy(begin, first, std::ostream_iterator<decltype(*first)>(os));
        }
    }

}

int my_mk_dir(std::string output_dir);

template<typename Representation, typename Period>
std::ostream &operator<<(std::ostream &os, std::chrono::duration<Representation, Period> duration) //时间消耗
{
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms";
    return os;
}
