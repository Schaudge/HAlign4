#pragma once
#include "../Utils/Utils.hpp"
#include "divsufsort.h"
#include "../Utils/Arguments.hpp"

#include <iostream>
#include <vector>

#include <iterator>
#include <algorithm>
#include <cstring>
#include <array>
#include <limits>
#include <unordered_map>
#include <set>

namespace suffix_array
{
    template<size_t width> 
    class SuffixArray
    {
    public:
        using triple = std::array<size_t, 3>;

        template<typename InputIterator>
        SuffixArray(InputIterator first, InputIterator last, unsigned char end_mark)
            : length(last- first + 1)
            , dis(1)
            , reword(_copy_reword(first, last, end_mark))
            , SA(new int32_t[length])
        {
            divsufsort(reword, SA, length, 4);
            endp = build_b();
            delete[] reword;
        }

        ~SuffixArray()
        {
            delete[] O;
            delete[] SA;
            delete[] B;
            delete[] begin;
        }
        
        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
        {
            size_t common_prefix_length = 0;
            int lbegin, lend, start, end, len_sub = last - first;
            char sub = *(first++); //第0个
            start = begin[(int)sub - 1];
            end = begin[(int)sub] - 1;
            //std::cout << start << " " << end << " s-e?\n";
            while (first < last)
            {
                sub = *(first);
                lbegin = find(sub, start, end);
                lend = rfind(sub, start, end);
                //std::cout << lbegin << " " << lend << "?\n";
                if (lbegin == -1)
                    break;
                start = begin[(int)sub - 1] + O_index_num(lbegin, sub) - 1;
                end = begin[(int)sub - 1] + O_index_num(lend, sub) - 1;
                //std::cout << start << " " << end << " ??\n";
                first++;
            }

            common_prefix_length = last - first;
            common_prefix_length = len_sub - common_prefix_length;
            if (common_prefix_length < threshold)
                return std::vector<size_t>();
            std::vector<size_t> starts{ common_prefix_length };
            for (int i = start; i <= end; i++)
            {
                //std::cout << SA[i] << " ";
                starts.emplace_back(length - 1 - SA[i] - common_prefix_length);
            }
            
            return std::move(starts);
        }


        template<typename RandomAccessIterator>
        std::vector<triple> get_common_substrings(RandomAccessIterator first, RandomAccessIterator last, size_t threshold)const
        {
            std::vector<triple> common_substrings; 
            const size_t rhs_len = last - first;
            if (rhs_len < threshold)
                return common_substrings;
            for (size_t rhs_index = 0; rhs_index < rhs_len; )
            {
                auto found = search_for_prefix(first + rhs_index, last, threshold); 
                
                if (found.empty())// || (tag_center && found.size() > 2))
                {
                    ++rhs_index;
                }
                else
                {
                    //std::cout << found.size() << "\n";
                    for (size_t i = 1; i != found.size(); ++i)
                        common_substrings.emplace_back(triple({found[i], rhs_index,found[0] })); 
                    rhs_index += found[0]- threshold + 1;
                    //rhs_index += found[0];
                }
            }
            
            return std::move(common_substrings);
        }

        int find(char now, int start, int end) const
        {
            for (int i = start; i <= end; i++)
                if (B[i] == now)
                    return i;
            return -1;
        }

        int rfind(char now, int start, int end) const
        {
            for (int i = end; i >= start; i--)
                if (B[i] == now)
                    return i;
            return -1;
        }

        int O_index_num(int x, char now) const
        {
            int num, i, quotient = x / dis;
            if (((x - quotient * dis) <= (dis / 2)) || (quotient == (Osize - 1)))
            {
                num = O[quotient][(int)now - 1];
                for (i = quotient * dis + 1; i <= x; i++)
                    if (B[i] == now)
                        num++;
            }
            else            
            {
                num = O[quotient + 1][(int)now - 1];
                for (i = (quotient + 1) * dis; i > x; i--)
                    if (B[i] == now)
                        num--;
            }
            return num;
        }

    private:

        template<typename InputIterator>
        unsigned char* _copy_reword(InputIterator first, InputIterator last, unsigned char end_mark)
        {
            unsigned char* result = new unsigned char[length];
            int i = length - 2;
            while (i >= 0)
            {
                result[i] = *(first++);
                i--;
            }
            result[length - 1] = end_mark;
            return result;
        }

        inline bool leq(int a1, int a2, int b1, int b2) { // lexic. orderfor pairs
            return(a1 < b1 || a1 == b1 && a2 <= b2);
        }                                                  // and triples
        inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
            return(a1 < b1 || a1 == b1 && leq(a2, a3, b2, b3));
        }
        // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K fromr
        static void radixPass(int* a, int* b, int* r, int n, int K)
        {// count occurrences
            int* c = new int[K + 1];                          // counter array
            for (int i = 0; i <= K; i++) c[i] = 0;         // resetcounters
            for (int i = 0; i < n; i++) c[r[a[i]]]++;    // countoccurences
            for (int i = 0, sum = 0; i <= K; i++) { // exclusive prefix sums
                int t = c[i];  c[i] = sum;  sum += t;
            }
            for (int i = 0; i < n; i++) b[c[r[a[i]]]++] = a[i];      //sort
            delete[] c;
        }
        // find the suffix array SA of s[0..n-1] in {1..K}^n
        // require s[n]=s[n+1]=s[n+2]=0, n>=2
        void suffixArray(int* s, int* SA, int n, int K) {
            int n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n02 = n0 + n2;
        
            int* s12 = new int[n02 + 3]; s12[n02] = s12[n02 + 1] = s12[n02 + 2] = 0;
            int* SA12 = new int[n02 + 3]; SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
            int* s0 = new int[n0];
            int* SA0 = new int[n0];

            // generatepositions of mod 1 and mod  2 suffixes
            // the"+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
            for (int i = 0, j = 0; i < n + (n0 - n1); i++) if (i % 3 != 0) s12[j++] = i;
    

            // lsb radix sortthe mod 1 and mod 2 triples
            radixPass(s12, SA12, s + 2, n02, K);
            radixPass(SA12, s12, s + 1, n02, K);
            radixPass(s12, SA12, s, n02, K);


            // findlexicographic names of triples
            int name = 0, c0 = -1, c1 = -1, c2 = -1;
            for (int i = 0; i < n02; i++) {
                if (s[SA12[i]] != c0 || s[SA12[i] + 1] != c1 || s[SA12[i] + 2] != c2) {
                    name++; c0 = s[SA12[i]];  c1 = s[SA12[i] + 1];  c2 = s[SA12[i] + 2];
                }

                if (SA12[i] % 3 == 1) { s12[SA12[i] / 3] = name; }// left half
                else { s12[SA12[i] / 3 + n0] = name; } // right half
   
            }

            // recurse if namesare not yet unique
            if (name < n02) {
    
                suffixArray(s12, SA12, n02, name);
                // store uniquenames in s12 using the suffix array
                for (int i = 0; i < n02; i++) s12[SA12[i]] = i + 1;
            }
            else // generate the suffix array of s12 directly
                for (int i = 0; i < n02; i++) SA12[s12[i] - 1] = i;




            // stably sort themod 0 suffixes from SA12 by their first character
            for (int i = 0, j = 0; i < n02; i++) if (SA12[i] < n0) s0[j++] = 3 * SA12[i];

            radixPass(s0, SA0, s, n0, K);


            // merge sorted SA0suffixes and sorted SA12 suffixes
            for (int p = 0, t = n0 - n1, k = 0; k < n; k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) *3 + 2)
                int i = GetI(); // pos of current offset 12 suffix
                int j = SA0[p]; // pos of current offset 0  suffix
                if (SA12[t] < n0 ?
                    leq(s[i], s12[SA12[t] + n0], s[j], s12[j / 3]) :
                    leq(s[i], s[i + 1], s12[SA12[t] - n0 + 1], s[j], s[j + 1], s12[j / 3 + n0]))
                { // suffix fromSA12 is smaller
                    SA[k] = i;  t++;
                    if (t == n02) { // done --- only SA0 suffixes left
                        for (k++; p < n0; p++, k++) SA[k] = SA0[p];
                    }
                }
                else {
                    SA[k] = j;  p++;
                    if (p == n0) { // done--- only SA12 suffixes left
                        for (k++; t < n02; t++, k++) SA[k] = GetI();
                    }
                }
            }
            delete[]s12; delete[] SA12; delete[] SA0; delete[] s0;
        }

        int* build_sa()
        {
            int* s = new int[length + 3];
            int* sa = new int[length + 3];
            for (int i = 0; i < length; i++) s[i] = (int)reword[i];
            s[length] = s[length + 1] = s[length + 2] = sa[length] = sa[length + 1] = sa[length + 2] = 0;
            suffixArray(s, sa, length, 5);
            delete[] s;
            return sa;
        }

        int build_b()
        {
            quadra* o = new quadra[length / dis + 2]();
            Osize = length / dis + 2;
            int* _begin = new int[5];
            _begin[4] = length;
            unsigned char* b = new unsigned char[length];
            quadra num = { 0,0,0,0 };
            int  i, e = 0;
            for (i = 0; i < length; i++)
            {
                b[i] = reword[(SA[i] + length - 1) % length];
                if (((int)b[i]) == 0)
                    e = i;
                else
                    num[(int)b[i] - 1]++;
                if (i % dis == 0)
                    o[i / dis] = num;
            }
            _begin[0] = 1;
            _begin[1] = num[0] + 1;
            _begin[2] = num[0] + num[1] + 1;
            _begin[3] = num[0] + num[1] + num[2] + 1;
            B = b;  
            begin = _begin; 
            O = o;
            return e; 
        }
    private:
        using quadra = std::array<int32_t, 4>;
        const int dis;
        int endp;
        int Osize;
    public:
        const size_t length;                        // including the ending '$'
        const unsigned char* reword;          
        int32_t* SA;                              
        const unsigned char* B;            
        const quadra* O;
        const int* begin;
    };
}
