#pragma once
#include "../SuffixArray/SuffixArray.hpp"  //利用后缀树
#include "../Utils/Utils.hpp"
#include "../multi-thread/multi.hpp"

#include <vector>
#include <array>
#include <string>
#include<algorithm>


namespace star_alignment //星比对命名空间
{

    class StarAligner//星比对类
    {
    private:
        using triple = std::array<size_t, 3>;
        using quadra = std::array<size_t, 4>;
        using sequence_type = std::vector<unsigned char>;

    public:
        static std::vector<sequence_type> align(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);
        static void get_gaps(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);

        static std::vector<triple> _optimal_path(const std::vector<triple>& common_substrings);
        static std::vector<int> _trace_back_bp(const std::vector<triple>& common_substrings, int* p);
        static std::vector<triple> _optimal_path_bp(const std::vector<triple>& optimal_common_substrings);

        
    private:
        StarAligner(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center);

        std::vector<sequence_type> _align() const;
        void _get_gaps() const;

        std::vector<size_t> _set_lengths() const;
        size_t _set_centre() const;

        // main steps
        std::vector<std::array<std::vector<utils::Insertion>, 2>> _pairwise_align() const;
        std::vector<std::vector<utils::Insertion>> _merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>>& pairwise_gaps) const;
        std::vector<sequence_type> _insert_gaps(const std::vector<std::vector<utils::Insertion>>& gaps) const;

        void mul_pairwise_align() const;
        void mul_fasta_func(int i, const suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER>& st,
            std::vector<std::array<std::vector<utils::Insertion>, 2>>& all_pairwise_gaps, int threshold1) const;
        
        // support
        static void _append(const std::vector<size_t>& src_gaps, std::vector<utils::Insertion>& des_gaps, size_t start);
        std::vector<std::vector<utils::Insertion>>& Insertions;
        std::vector<sequence_type>& _sequences;
        const size_t _row;
        std::vector<size_t> _lengths;
        size_t thresh1;
        size_t _centre;
        size_t _centre_len;
    };

}
