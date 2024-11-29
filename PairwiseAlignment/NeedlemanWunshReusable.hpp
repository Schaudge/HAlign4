#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <tuple>
#include <cmath>
#include "../SuffixArray/SuffixArray.hpp"
#include "../Utils/Utils.hpp"
#include <iostream>
#include <string>
#include "WFA2-lib/bindings/cpp/WFAligner.hpp"

const int thresh0 = 10001;

using RandomAccessIterator = std::vector<unsigned char>::const_iterator;
using gap_vector_type = std::vector<int>;
using char_vector_type = std::vector<unsigned char>;
using triple = std::array<int, 4>;
using quadra = std::array<int, 5>;
using insert = std::vector<std::tuple<int, int>>;
using in = std::tuple<int, int>;

class Kband {
public:
    int match;
    int mismatch;
    int d;
    int e;
    int my_INT_MIN;
    unsigned char* A;
    unsigned char* B;
    unsigned char* Aq;
    unsigned char* Bq;
    int** pm;
    int** pm2;
    int** pmt1;
    int** pmt2;
    int** pmt;

    unsigned char** bt;
    std::vector<unsigned char> seq_A;
    std::vector<unsigned char> seq_B;

public:
    Kband();

    ~Kband();
    inline int score(unsigned char xi, unsigned char yi);
    inline bool InsiderStrip(int i, int j, int k, int diff = 0);
    inline int index(int i, int l);
    inline int maxi(int a, int b);
    inline int maxi(int a, int b, int c);
    void Init(int m, int k, int diff);
    void InitTwo(int ii, int k, int diff);
    int ChooseWay(int p0, int p1, int p2, bool state = true);
    inline int parse(int b, int s);

    std::tuple<insert, insert>
        PSA_AGP_Kband3(const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end,
            int cmatch = 1, int cmismatch = -2, int cd = 3, int ce = 1);
};
std::tuple<insert, insert> mywfa(wfa::WFAlignerGapAffine& aligner, const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end);
std::tuple<insert, insert> parseCigar(const std::string& cigar);
std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin);
std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin, bool tag);
void insertGaps(std::string& str, const insert& insertions);