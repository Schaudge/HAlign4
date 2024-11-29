#include "StarAligner.hpp"
#include "../Utils/Pseudo.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Graph.hpp"
#include "../PairwiseAlignment/NeedlemanWunshReusable.hpp"
//std::map<std::thread::id, wfa::WFAlignerGapAffine*> threadMAP;
std::vector<std::vector<unsigned char>> star_alignment::StarAligner::align(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh,int center)
{
    return StarAligner(insertions, sequences, thresh, center)._align();
}

void star_alignment::StarAligner::get_gaps(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center)
{
    StarAligner(insertions, sequences, thresh, center)._get_gaps();
}

star_alignment::StarAligner::StarAligner(std::vector<std::vector<utils::Insertion>>& insertions, std::vector<sequence_type>& sequences, size_t thresh, int center)
    : thresh1(thresh)
    , Insertions(insertions)
    , _sequences(sequences)
    , _row(_sequences.size())
    , _lengths(_set_lengths())
    , _centre(_set_centre())
    , _centre_len(_sequences[_centre].size())
{
    if (center != -1)
    {
        _centre = center;
        _centre_len = _sequences[_centre].size();
    }
}

std::vector<size_t> star_alignment::StarAligner::_set_lengths() const
{
    std::vector<size_t> lengths(_row);

    for (size_t i = 0; i != _row; ++i) lengths[i] = _sequences[i].size();
    return lengths;
}

size_t star_alignment::StarAligner::_set_centre() const
{
    size_t centre_index = 0;
    for (size_t i = 1; i != _row; ++i)
        if (_lengths[i] > _lengths[centre_index])
            centre_index = i;
    return centre_index;
}

std::vector<std::vector<unsigned char>> star_alignment::StarAligner::_align() const
{
    return _insert_gaps(_merge_results(_pairwise_align()));
}

void star_alignment::StarAligner::_get_gaps() const
{
    //if (threadPool0->Thread_num == 1)
    //    return _merge_results(_pairwise_align());
    mul_pairwise_align();
    //_merge_results(mul_pairwise_align());
}

auto star_alignment::StarAligner::_pairwise_align() const -> std::vector<std::array<std::vector<utils::Insertion>, 2>>
{

    suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER> st(_sequences[_centre].cbegin(), _sequences[_centre].cend(),nucleic_acid_pseudo::end_mark);//实例化后缀树对象st，比对同源区域
    std::vector<std::array<std::vector<utils::Insertion>, 2>> all_pairwise_gaps;
    //!!!
    wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    
    for (size_t i = 0; i != _row; ++i)
    {
        //auto common_substrings = _optimal_path_bp(st.get_common_substrings(_sequences[i].cbegin(), _sequences[i].cend(), thresh1));
        auto common_substrings = _optimal_path(st.get_common_substrings(_sequences[i].cbegin(), _sequences[i].cend(), thresh1));
        
        std::vector<quadra> intervals;
        intervals.reserve(common_substrings.size() + 1);

        if (common_substrings.empty())
        {
            intervals.emplace_back(quadra({ 0, _centre_len, 0, _sequences[i].size() }));
        }
        else
        {
            if (common_substrings[0][0] || common_substrings[0][1])
                intervals.emplace_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

            for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                    common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                    intervals.emplace_back(quadra
                    ({
                        common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                        common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                        }));

            if (common_substrings.back()[0] + common_substrings.back()[2] != _centre_len ||
                common_substrings.back()[1] + common_substrings.back()[2] != _lengths[i])
                intervals.emplace_back(quadra
                ({
                    common_substrings.back()[0] + common_substrings.back()[2], _centre_len,
                    common_substrings.back()[1] + common_substrings.back()[2], _lengths[i]
                    }));
        }

        std::array<std::vector<utils::Insertion>, 2> pairwise_gaps;
        for (size_t j = 0; j != intervals.size(); ++j)
        {
            const size_t centre_begin = intervals[j][0];
            const size_t centre_end = intervals[j][1];
            const size_t sequence_begin = intervals[j][2];
            const size_t sequence_end = intervals[j][3];
            //std::cout << i<<"\t" << centre_begin << "\t" << centre_end << "\t" << sequence_begin << "\t" << sequence_end << "\n";
            //std::cout << i << "\t" << centre_end-centre_begin <<  "\t"  << sequence_end-sequence_begin << "\n";

            auto [lhs_gaps, rhs_gaps] = mywfa(aligner, _sequences[_centre], centre_begin, centre_end,  //非同源区域，动态规划比对
                _sequences[i], sequence_begin, sequence_end);

            int an = 0, bn = 0;
            for (int ii = 0; ii < lhs_gaps.size(); ii++)
            {
                an += std::get<1>(lhs_gaps[ii]);
                if ((!pairwise_gaps[0].empty()) && pairwise_gaps[0].back().index == std::get<0>(lhs_gaps[ii]))
                    pairwise_gaps[0].back().number += std::get<1>(lhs_gaps[ii]);
                else
                    pairwise_gaps[0].emplace_back(utils::Insertion({ (size_t)std::get<0>(lhs_gaps[ii]),(size_t)std::get<1>(lhs_gaps[ii]) }));
            }
            for (int ii = 0; ii < rhs_gaps.size(); ii++)
            {
                bn += std::get<1>(rhs_gaps[ii]);
                if ((!pairwise_gaps[1].empty()) && pairwise_gaps[1].back().index == std::get<0>(rhs_gaps[ii]))
                    pairwise_gaps[1].back().number += std::get<1>(rhs_gaps[ii]);
                else
                    pairwise_gaps[1].emplace_back(utils::Insertion({ (size_t)std::get<0>(rhs_gaps[ii]),(size_t)std::get<1>(rhs_gaps[ii]) }));
            }
        }

        size_t sequence_gap_num = 0;
        for (auto gap : pairwise_gaps[0]) sequence_gap_num += gap.number;

        all_pairwise_gaps.emplace_back(pairwise_gaps);
        sequence_type().swap(_sequences[i]);
    }

    return all_pairwise_gaps;
}

auto star_alignment::StarAligner::_optimal_path(const std::vector<triple>& common_substrings) 
-> std::vector<triple>
{
    std::vector<triple> optimal_common_substrings;
    if (common_substrings.empty()) return optimal_common_substrings;

    const size_t pair_num = common_substrings.size();
    utils::AdjacencyList graph(pair_num + 1);

    for (size_t i = 0; i != pair_num; ++i)
        for (size_t j = 0; j != pair_num; ++j)
            if (i != j && common_substrings[i][0] + common_substrings[i][2] < common_substrings[j][0] + common_substrings[j][2]
                && common_substrings[i][1] + common_substrings[i][2] < common_substrings[j][1] + common_substrings[j][2])
            {
                const int possible_overlap = std::max(
                    static_cast<int>(common_substrings[i][0] + common_substrings[i][2]) - static_cast<int>(common_substrings[j][0]),
                    static_cast<int>(common_substrings[i][1] + common_substrings[i][2]) - static_cast<int>(common_substrings[j][1]));

                unsigned weight = common_substrings[j][2];
                if (possible_overlap > 0) weight -= possible_overlap;
                graph.add_edge(i + 1, j + 1, weight);
            }

    for (size_t i = 0; i != pair_num; ++i)
        graph.add_edge(0, i + 1, common_substrings[i][2]);

    const auto optimal_path = graph.get_longest_path();

    optimal_common_substrings.reserve(optimal_path.size());
    optimal_common_substrings.emplace_back(triple({ common_substrings[optimal_path[0] - 1][0],
                                                     common_substrings[optimal_path[0] - 1][1],
                                                     common_substrings[optimal_path[0] - 1][2] }));

    for (size_t i = 0; i < optimal_path.size() - 1; ++i)
    {
        size_t new_len = graph.get_weight(optimal_path[i], optimal_path[i + 1]);
        size_t old_len = common_substrings[optimal_path[i + 1] - 1][2];
        int difference = static_cast<int>(old_len) - static_cast<int>(new_len);

        size_t lhs_first = common_substrings[optimal_path[i + 1] - 1][0];
        size_t rhs_first = common_substrings[optimal_path[i + 1] - 1][1];
        if (difference > 0)
        {
            lhs_first += difference; rhs_first += difference;
        }

        optimal_common_substrings.emplace_back(triple({ lhs_first, rhs_first, new_len }));
    }

    return optimal_common_substrings;
}

void star_alignment::StarAligner::_append(const std::vector<size_t>& src_gaps, std::vector<utils::Insertion>& des_gaps, size_t start)
{
    for (size_t i = 0; i != src_gaps.size(); ++i)
        if (src_gaps[i])
        {
            if (des_gaps.size() && des_gaps.back().index == start + i)
                des_gaps.back().number += src_gaps[i];
            else
                des_gaps.emplace_back(utils::Insertion({ start + i, src_gaps[i] }));
        }
}

auto star_alignment::StarAligner::_merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>>& pairwise_gaps) const
-> std::vector<std::vector<utils::Insertion>>
{
    std::vector<utils::Insertion> final_centre_gaps;
    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        for (size_t lhs_pointer = 0, rhs_pointer = 0; rhs_pointer != curr_centre_gaps.size(); )
        {
            if (lhs_pointer == final_centre_gaps.size())
            {
                final_centre_gaps.insert(final_centre_gaps.cend(), curr_centre_gaps.cbegin() + rhs_pointer, curr_centre_gaps.cend());
                break;
            }

            if (final_centre_gaps[lhs_pointer].index == curr_centre_gaps[rhs_pointer].index)
            {
                if (final_centre_gaps[lhs_pointer].number < curr_centre_gaps[rhs_pointer].number)
                    final_centre_gaps[lhs_pointer].number = curr_centre_gaps[rhs_pointer].number;
                ++lhs_pointer;
                ++rhs_pointer;
            }
            else if (final_centre_gaps[lhs_pointer].index < curr_centre_gaps[rhs_pointer].index)
            {
                ++lhs_pointer;
            }
            else
            {
                final_centre_gaps.insert(final_centre_gaps.cbegin() + lhs_pointer, curr_centre_gaps[rhs_pointer]);
                ++lhs_pointer; // because of the insert above
                ++rhs_pointer;
            }
        }
    }

    std::vector<std::vector<utils::Insertion>> final_sequence_gaps;
    final_sequence_gaps.reserve(_row);
    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        const auto& curr_sequence_gaps = pairwise_gaps[i][1];

        std::vector<utils::Insertion> centre_addition;
        centre_addition.reserve(final_centre_gaps.size());
        utils::Insertion::minus(final_centre_gaps.cbegin(), final_centre_gaps.cend(),
            curr_centre_gaps.cbegin(), curr_centre_gaps.cend(),
            std::back_inserter(centre_addition));

        std::vector<utils::Insertion> sequence_addition;
        for (size_t centre_index = 0, sequence_index = 0, centre_gaps_index = 0, sequence_gaps_index = 0, centre_addition_index = 0;
            centre_addition_index != centre_addition.size(); ++centre_addition_index)
        {
            const auto curr_addition = centre_addition[centre_addition_index]; // current addition pending process

            while (centre_index < curr_addition.index)
            {
                size_t centre_distance = centre_gaps_index < curr_centre_gaps.size() ?
                    curr_centre_gaps[centre_gaps_index].index - centre_index : std::numeric_limits<size_t>::max();
                size_t sequence_distance = sequence_gaps_index < curr_sequence_gaps.size() ?
                    curr_sequence_gaps[sequence_gaps_index].index - sequence_index : std::numeric_limits<size_t>::max();

                size_t step = std::min({ sequence_distance, centre_distance, curr_addition.index - centre_index }); // assure centre_index <= curr_addtion.index
                centre_index += step;
                sequence_index += step;

                if (centre_gaps_index < curr_centre_gaps.size() && curr_centre_gaps[centre_gaps_index].index == centre_index)
                    sequence_index += curr_centre_gaps[centre_gaps_index++].number;

                else if (sequence_gaps_index < curr_sequence_gaps.size() && curr_sequence_gaps[sequence_gaps_index].index == sequence_index)
                    centre_index += curr_sequence_gaps[sequence_gaps_index++].number;
            }

            if (sequence_addition.size() && sequence_index == sequence_addition.back().index)
                sequence_addition.back().number += curr_addition.number;
            else
                sequence_addition.emplace_back(utils::Insertion({ sequence_index, curr_addition.number }));
        }

        std::vector<utils::Insertion> indels_of_current_sequence;
        indels_of_current_sequence.reserve(curr_sequence_gaps.size() + sequence_addition.size());
        utils::Insertion::plus(curr_sequence_gaps.cbegin(), curr_sequence_gaps.cend(),
            sequence_addition.cbegin(), sequence_addition.cend(),
            std::back_inserter(indels_of_current_sequence));
        final_sequence_gaps.emplace_back(indels_of_current_sequence);
    }

    return final_sequence_gaps;
}

auto star_alignment::StarAligner::_insert_gaps(const std::vector<std::vector<utils::Insertion>>& gaps) const
-> std::vector<sequence_type>
{
    size_t length = _lengths[0];
    for (const auto gap : gaps[0]) length += gap.number;

    std::vector<sequence_type> aligned;
    aligned.reserve(_row);

    for (size_t i = 0; i != _row; ++i)
    {
        aligned.emplace_back(length);
        utils::Insertion::insert_gaps(_sequences[i].cbegin(), _sequences[i].cend(),
            gaps[i].cbegin(), gaps[i].cend(), aligned.back().begin(), nucleic_acid_pseudo::GAP);
    }
    return aligned;
}

void star_alignment::StarAligner::mul_pairwise_align() const
{
    cout_cur_time();
    const auto align_1 = std::chrono::high_resolution_clock::now(); //记录比对起始时间
    std::cout << "Start: build Suffix Array No." << _centre << "\n";
    suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER> st(_sequences[_centre].cbegin(), _sequences[_centre].cend(), nucleic_acid_pseudo::end_mark);//实例化后缀树对象st，比对同源区域
    //size_t peakMem2 = get_peak_memory(); // 获取被测试代码执行后当前进程的内存占用峰值
    cout_cur_time();
    std::cout << "End  : build Suffix Array " << (std::chrono::high_resolution_clock::now() - align_1)<<"\n";
    std::vector<std::array<std::vector<utils::Insertion>, 2>> pairwise_gaps(_row); //定义变量-接收结果

    //for (int i = 0; i < threadPool0->Thread_num; i++)
    //    threadMAP[threadPool0->workers[i].get_id()] = new wfa::WFAlignerGapAffine(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
    cout_cur_time();
    std::cout << "Start: pairwise sequence alignment\n";
    for (int i = 0; i != _row; ++i)
        threadPool0->execute(&star_alignment::StarAligner::mul_fasta_func, this, i, std::ref(st), std::ref(pairwise_gaps), thresh1);
    threadPool0->waitFinished();

    _sequences.clear();
 
    delete threadPool0;

    EmptySet();
    cout_cur_time();
    std::cout << "End  : pairwise sequence alignment\n";
    //ns.close();
    //for (int i = 0; i < threadPool0->Thread_num; i++)
    //    delete threadMAP[threadPool0->workers[i].get_id()];

    std::vector<utils::Insertion> final_centre_gaps;
    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        for (size_t lhs_pointer = 0, rhs_pointer = 0; rhs_pointer != curr_centre_gaps.size(); )
        {
            if (lhs_pointer == final_centre_gaps.size())
            {
                final_centre_gaps.insert(final_centre_gaps.cend(), curr_centre_gaps.cbegin() + rhs_pointer, curr_centre_gaps.cend());
                break;
            }

            if (final_centre_gaps[lhs_pointer].index == curr_centre_gaps[rhs_pointer].index)
            {
                if (final_centre_gaps[lhs_pointer].number < curr_centre_gaps[rhs_pointer].number)
                    final_centre_gaps[lhs_pointer].number = curr_centre_gaps[rhs_pointer].number;
                ++lhs_pointer;
                ++rhs_pointer;
            }
            else if (final_centre_gaps[lhs_pointer].index < curr_centre_gaps[rhs_pointer].index)
            {
                ++lhs_pointer;
            }
            else
            {
                final_centre_gaps.insert(final_centre_gaps.cbegin() + lhs_pointer, curr_centre_gaps[rhs_pointer]);
                ++lhs_pointer; // because of the insert above
                ++rhs_pointer;
            }
        }
    }

    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        const auto& curr_sequence_gaps = pairwise_gaps[i][1];

        std::vector<utils::Insertion> centre_addition;
        centre_addition.reserve(final_centre_gaps.size());
        utils::Insertion::minus(final_centre_gaps.cbegin(), final_centre_gaps.cend(),
            curr_centre_gaps.cbegin(), curr_centre_gaps.cend(),
            std::back_inserter(centre_addition));

        std::vector<utils::Insertion> sequence_addition;
        for (size_t centre_index = 0, sequence_index = 0, centre_gaps_index = 0, sequence_gaps_index = 0, centre_addition_index = 0;
            centre_addition_index != centre_addition.size(); ++centre_addition_index)
        {
            const auto curr_addition = centre_addition[centre_addition_index]; // current addition pending process

            while (centre_index < curr_addition.index)
            {
                size_t centre_distance = centre_gaps_index < curr_centre_gaps.size() ?
                    curr_centre_gaps[centre_gaps_index].index - centre_index : std::numeric_limits<size_t>::max();
                size_t sequence_distance = sequence_gaps_index < curr_sequence_gaps.size() ?
                    curr_sequence_gaps[sequence_gaps_index].index - sequence_index : std::numeric_limits<size_t>::max();

                size_t step = std::min({ sequence_distance, centre_distance, curr_addition.index - centre_index }); // assure centre_index <= curr_addtion.index
                centre_index += step;
                sequence_index += step;

                if (centre_gaps_index < curr_centre_gaps.size() && curr_centre_gaps[centre_gaps_index].index == centre_index)
                    sequence_index += curr_centre_gaps[centre_gaps_index++].number;

                else if (sequence_gaps_index < curr_sequence_gaps.size() && curr_sequence_gaps[sequence_gaps_index].index == sequence_index)
                    centre_index += curr_sequence_gaps[sequence_gaps_index++].number;
            }

            if (sequence_addition.size() && sequence_index == sequence_addition.back().index)
                sequence_addition.back().number += curr_addition.number;
            else
                sequence_addition.emplace_back(utils::Insertion({ sequence_index, curr_addition.number }));
        }

        std::vector<utils::Insertion>& indels_of_current_sequence = Insertions[i];
        indels_of_current_sequence.reserve(curr_sequence_gaps.size() + sequence_addition.size());
        utils::Insertion::plus(curr_sequence_gaps.cbegin(), curr_sequence_gaps.cend(),
            sequence_addition.cbegin(), sequence_addition.cend(),
            std::back_inserter(indels_of_current_sequence));
    }
    //return std::move(all_pairwise_gaps); //双序列比对得到的两两gap，长度为n的vector，每个元素有长度为2的array，每个元素是有若干个Insertion的vector
}

void star_alignment::StarAligner::mul_fasta_func(int i, const suffix_array::SuffixArray<nucleic_acid_pseudo::NUMBER>& st,
    std::vector<std::array<std::vector<utils::Insertion>, 2>>& all_pairwise_gaps, int threshold1) const
{
    if (i == _centre)
        return;
    wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    //Kband* kband = new Kband();
    auto common_substrings = _optimal_path(st.get_common_substrings(_sequences[i].cbegin(), _sequences[i].cend(), thresh1));

    std::vector<quadra> intervals;
    intervals.reserve(common_substrings.size() + 1);

    if (common_substrings.empty())
    {
        intervals.emplace_back(quadra({ 0, _centre_len, 0, _sequences[i].size() }));
    }
    else
    {
        if (common_substrings[0][0] || common_substrings[0][1])
            intervals.emplace_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

        for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
            if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                intervals.emplace_back(quadra
                ({
                    common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                    common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                    }));

        if (common_substrings.back()[0] + common_substrings.back()[2] != _centre_len ||
            common_substrings.back()[1] + common_substrings.back()[2] != _lengths[i])
            intervals.emplace_back(quadra
            ({
                common_substrings.back()[0] + common_substrings.back()[2], _centre_len,
                common_substrings.back()[1] + common_substrings.back()[2], _lengths[i]
                }));
    }

    std::array<std::vector<utils::Insertion>, 2> pairwise_gaps;
    for (size_t j = 0; j != intervals.size(); ++j)
    {
        const size_t centre_begin = intervals[j][0];
        const size_t centre_end = intervals[j][1];
        const size_t sequence_begin = intervals[j][2];
        const size_t sequence_end = intervals[j][3];
        //std::cout << i<<"\t" << centre_begin << "\t" << centre_end << "\t" << sequence_begin << "\t" << sequence_end << "\n";
        //std::cout << i << "\t" << centre_end - centre_begin << "\t" << sequence_end - sequence_begin << "\n";

        auto [lhs_gaps, rhs_gaps] = mywfa(aligner, _sequences[_centre], centre_begin, centre_end,  //非同源区域，动态规划比对
            _sequences[i], sequence_begin, sequence_end);

        //auto [lhs_gaps, rhs_gaps] = kband->PSA_AGP_Kband3(_sequences[_centre], centre_begin, centre_end,  //非同源区域，动态规划比对
        //    _sequences[i], sequence_begin, sequence_end); //分治到比thresh小，然后k-band

        int an = 0, bn = 0;
        for (int ii = 0; ii < lhs_gaps.size(); ii++)
        {
            an += std::get<1>(lhs_gaps[ii]);
            if ((!pairwise_gaps[0].empty()) && pairwise_gaps[0].back().index == std::get<0>(lhs_gaps[ii]))
                pairwise_gaps[0].back().number += std::get<1>(lhs_gaps[ii]);
            else
                pairwise_gaps[0].emplace_back(utils::Insertion({ (size_t)std::get<0>(lhs_gaps[ii]),(size_t)std::get<1>(lhs_gaps[ii]) }));
        }
        for (int ii = 0; ii < rhs_gaps.size(); ii++)
        {
            bn += std::get<1>(rhs_gaps[ii]);
            if ((!pairwise_gaps[1].empty()) && pairwise_gaps[1].back().index == std::get<0>(rhs_gaps[ii]))
                pairwise_gaps[1].back().number += std::get<1>(rhs_gaps[ii]);
            else
                pairwise_gaps[1].emplace_back(utils::Insertion({ (size_t)std::get<0>(rhs_gaps[ii]),(size_t)std::get<1>(rhs_gaps[ii]) }));
        }
    }
    //delete kband;
    all_pairwise_gaps[i].swap(pairwise_gaps);
    _sequences[i].clear();
    return;
}

std::vector<int> star_alignment::StarAligner::_trace_back_bp(const std::vector<triple>& common_substrings, int* p)
{
    std::vector<int> ansi;
    int* tmp = p;
    for (int i = 0; i < common_substrings.size(); i++)
    {
        if (*p < tmp[i])
            p = &tmp[i];
    }
    int j, i = p - tmp;
    ansi.emplace_back(i);
    p = tmp;
    while (i > 0)
    {
        j = i - 1;
        if (p[i] == common_substrings[i][2]) break;
        while (j >= 0)
        {

            if ((common_substrings[i][0] >= (common_substrings[j][0] + common_substrings[j][2])) && (p[i] == (p[j] + common_substrings[i][2])))
            {
                ansi.emplace_back(j);
                i = j;
                break;
            }
            j--;
        }
    }
    reverse(ansi.begin(), ansi.end());//反转
    return ansi;
}
//第二步，依据动态规划，选出合适的不重叠的同源区段
auto star_alignment::StarAligner::_optimal_path_bp(const std::vector<triple>& optimal_common_substrings) //最优路径
-> std::vector<triple>
{
    std::vector<triple> ans_common_substrings;
    int m = optimal_common_substrings.size();
    if (m <= 1)
    {
        if (m == 1)
            ans_common_substrings.emplace_back(optimal_common_substrings[0]);
        return ans_common_substrings;
    }
    int* p = new int[m];
    for (int i = 0; i < m; i++) p[i] = optimal_common_substrings[i][2];
    for (int i = 1; i < m; i++)
        for (int j = 0; j < i; j++)
            if (optimal_common_substrings[i][0] >= (optimal_common_substrings[j][0] + optimal_common_substrings[j][2]))
                p[i] = (p[i] > (p[j] + optimal_common_substrings[i][2])) ? p[i] : (p[j] + optimal_common_substrings[i][2]);
    std::vector<int> ansi = _trace_back_bp(optimal_common_substrings, p);


    for (int i = 0; i < ansi.size(); i++) ans_common_substrings.emplace_back(optimal_common_substrings[ansi[i]]);
    std::vector<int>().swap(ansi);
    delete[] p;
    return ans_common_substrings;
}
