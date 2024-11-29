#include "NeedlemanWunshReusable.hpp"
//kband
Kband::Kband() {
    match = 1;
    mismatch = -2;
    d = 3;
    e = 1;
    my_INT_MIN = 0; // 初始化为合适的值
    A = new unsigned char[thresh0 + 2];
    B = new unsigned char[thresh0 + 2];
    pm = new int* [3];
    pm2 = new int* [3];
    for (int i = 0; i < 3; ++i) {
        pm[i] = new int[thresh0];
        pm2[i] = new int[thresh0];
    }
    //pm = new int[3][thresh0];
    //pm2 = new int[3][thresh0];
    pmt1 = pm;
    pmt2 = pm2;
    pmt = pm;
    bt = new unsigned char* [thresh0];
    for (int i = 0; i < thresh0; ++i) {
        bt[i] = new unsigned char[thresh0];
    }
    //bt = new unsigned char[thresh0][thresh0];
    seq_A.resize(thresh0);
    seq_B.resize(thresh0);
}
Kband::~Kband() {
    delete[] A;
    delete[] B;
    for (int i = 0; i < 3; ++i) {
        delete[] pmt1[i];
        delete[] pmt2[i];
    }
    delete[] pmt1;
    delete[] pmt2;
    for (int i = 0; i < thresh0; ++i) 
        delete[] bt[i];
    delete[] bt;
    char_vector_type().swap(seq_A);
    char_vector_type().swap(seq_B);
}
inline int Kband::score(unsigned char xi, unsigned char yi)
{
    if (xi == yi)
        return match;
    else
        return mismatch;
}
inline bool Kband::InsiderStrip(int i, int j, int k, int diff)
{
    return ((-k <= (j - i)) && ((j - i) <= (k + diff)));
}
inline int Kband::index(int i, int l)
{
    return (i + l) % l;
}
inline int Kband::maxi(int a, int b)
{
    if (a > b)return a;
    else return b;
}
inline int Kband::maxi(int a, int b, int c)
{
    int max; if (a > b)max = a; else max = b;
    if (max > c)return max; else return c;
}
void Kband::Init(int m, int k, int diff)
{
    for (int i = 0; i < (m + 1); i++)
    {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            bt[i][j] = '\0';
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            pm[i][j] = my_INT_MIN;
    }
    pm[0][k] = 0;
    bt[0][k] = '\16';
    for (int j = 0; j < (diff + k + 1); j++)
    {
        pm[1][j + k] = -d - e * (j - 1);
        bt[0][j + k] = (char)8;
    }
    for (int i = 0; i < (k + 1); i++)
        bt[i][k - i] = '\3';
}
void Kband::InitTwo(int ii, int k, int diff)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < (diff + 2 * k + 1); j++)
            pm2[i][j] = my_INT_MIN;
    }
    if (ii < k + 1)
        pm2[2][index(k - ii, diff + 2 * k + 1)] = -d - e * (ii - 1);
}
int Kband::ChooseWay(int p0, int p1, int p2, bool state)
{
    if (p0 >= p1)
    {
        if (p0 >= p2)
            return state ? 16 : 0;
        else
            return state ? 48 : 2;
    }
    else if (p1 >= p2)
        return state ? 32 : 1;
    else
        return state ? 48 : 2;
}
inline int Kband::parse(int b, int s)
{
    //b = (int)b;
    b = (b >> (4 - s * 2));
    return (b & 3) - 1;
}

std::tuple<insert, insert>
    Kband::PSA_AGP_Kband3(const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end,
        int cmatch, int cmismatch, int cd, int ce)
{

    match = cmatch;
    mismatch = cmismatch;
    d = cd;
    e = ce;
    pm = pmt1;
    pm2 = pmt2;
    int i = 0, j, z, diff, k = 1, m, n, b_j, l, old, bt1, bt2, bt3, newi, l2, channel;
    bool state_ex = false;//交换标识
    int a_len = a_end - a_begin;
    int b_len = b_end - b_begin;
    //std::cout << a_len << " " << b_len << "     1\n";
    insert a_gap;
    insert b_gap;
    if (a_len == 0)                             //a为0
    {
        a_gap.emplace_back(in(a_begin, b_len));
        return std::make_tuple(std::move(a_gap), std::move(b_gap));
    }
    else if (b_len == 0)                        //b为0
    {
        b_gap.emplace_back(in(b_begin, a_len));
        return std::make_tuple(std::move(a_gap), std::move(b_gap));
    }

    A = utils::copy_DNA(sequence1, A, a_begin, a_end);
    B = utils::copy_DNA(sequence2, B, b_begin, b_end);

    if (a_len > b_len)  //保证  B长，A短
    {
        auto tmp = A; A = B; B = tmp;
        diff = a_len - b_len;
        i = a_len; a_len = b_len; b_len = i;
        state_ex = true; //交换标识
    }
    else
        diff = b_len - a_len;
    m = a_len, n = b_len;
    my_INT_MIN = old = -d * n;

    while (k <= m)
    {// init
        //Init(bt, pm, m, k, diff);
        pm = pmt1;
        pm2 = pmt2;
        for (i = 0; i < (m + 1); i++)
        {
            for (int j = 0; j < (diff + 2 * k + 1); j++)
                bt[i][j] = '\0';
        }
        for (i = 0; i < 3; i++)
        {
            for (int j = 0; j < (diff + 2 * k + 1); j++)
                pm[i][j] = my_INT_MIN;
        }
        pm[0][k] = 0;
        bt[0][k] = '\16';
        for (j = 0; j < (diff + k + 1); j++)
        {
            pm[1][j + k] = -d - e * (j - 1);
            bt[0][j + k] = (char)8;
        }
        for (i = 0; i < (k + 1); i++)
            bt[i][k - i] = '\3';
        l = diff + 2 * k + 1;
        //end-init
        for (i = 1; i < (m + 1); i++)
        {
            //InitTwo(pm2, i, k, diff);
            for (int q = 0; q < 3; q++)
            {
                for (j = 0; j < (diff + 2 * k + 1); j++)
                    pm2[q][j] = my_INT_MIN;
            }
            if (i < k + 1)
                pm2[2][index(k - i, diff + 2 * k + 1)] = -d - e * (i - 1);
            l2 = diff + 2 * k + 1;
            //end-init
            for (int z = -k; z < (diff + k + 1); z++)
            {
                j = z;
                if ((1 <= (j + i)) && ((j + i) <= n))
                {
                    j = j + k;
                    bt1 = bt2 = bt3 = 0;
                    bt1 = ChooseWay(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]);
                    pm2[0][index(j, l2)] = maxi(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]) + score(A[i - 1], B[j + i - k - 1]);
                    if (InsiderStrip(i, j + i - k - 1, k, diff))// x : B[j] ~_
                    {
                        pm2[1][index(j, l2)] = maxi(pm2[0][index(j - 1, l2)] - d, pm2[1][index(j - 1, l2)] - e);
                        if ((pm2[0][index(j - 1, l2)] - d) > (pm2[1][index(j - 1, l2)] - e)) bt2 = 4;
                        else bt2 = 8;
                    }
                    if (InsiderStrip(i - 1, j + i - k, k, diff))// y : A[i] ~_
                    {
                        pm2[2][index(j, l2)] = maxi(pm[0][index(j + 1, l)] - d, pm[2][index(j + 1, l)] - e);
                        if ((pm[0][index(j + 1, l)] - d) > (pm[2][index(j + 1, l)] - e)) bt3 = 1;
                        else bt3 = 3;
                    }
                    bt[i][index(j, l)] = (char)(bt1 + bt2 + bt3);
                }
            }
            pmt = pm;
            pm = pm2;
            pm2 = pmt;
        }
        newi = maxi(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k]);
        if (old == newi || (k * 2) > m) break;
        else { old = newi; k *= 2; }
    }
    channel = ChooseWay(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k], false);

    //traceback
    i = m;
    b_j = n;
    j = diff + k;

    seq_A.clear();
    seq_B.clear();

    while (i > 0 || j > k)
    {
        if (channel == 0)
        {
            channel = parse(bt[i][j], 0);
            seq_A.emplace_back(A[--i]);
            seq_B.emplace_back(B[--b_j]);
        }
        else if (channel == 1)
        {
            channel = parse(bt[i][j], 1);
            seq_A.emplace_back('\7');
            seq_B.emplace_back(B[b_j - 1]);
            b_j--;
            j--;
        }
        else if (channel == 2)
        {
            channel = parse(bt[i][j], 2);
            seq_A.emplace_back(A[i - 1]);
            seq_B.emplace_back('\7');
            i--;
            j++;
        }
        else
        {
            std::cout << "channel error!\n";
            exit(-1);
        }
    }
    //std::cout << seq_A.size() << " " << seq_B.size() << "\n";
    int j1 = 0, j2 = 0, num1 = 0, num2 = 0, match_num = 0;
    if (state_ex)
    {
        for (i = seq_A.size() - 1; i >= 0; i--)
        {
            if (seq_A[i] == '\7') { num1++; }
            else { if (num1 != 0)b_gap.emplace_back(in(b_begin + j1, num1)); num1 = 0; j1++; }
            if (seq_B[i] == '\7') { num2++; }
            else
            {
                if (num2 != 0)a_gap.emplace_back(in(a_begin + j2, num2)); num2 = 0; j2++;
                if (seq_A[i] == seq_B[i]) match_num++;
            }
        }
        if (num1 != 0)b_gap.emplace_back(in(b_begin + j1, num1));
        if (num2 != 0)a_gap.emplace_back(in(a_begin + j2, num2));
    }
    else
    {
        for (i = seq_A.size() - 1; i >= 0; i--)
        {
            if (seq_A[i] == '\7') { num1++; }
            else { if (num1 != 0)a_gap.emplace_back(in(a_begin + j1, num1)); num1 = 0; j1++; }
            if (seq_B[i] == '\7') { num2++; }
            else
            {
                if (num2 != 0)b_gap.emplace_back(in(b_begin + j2, num2)); num2 = 0; j2++;
                if (seq_A[i] == seq_B[i]) match_num++;
            }
        }
        if (num1 != 0)a_gap.emplace_back(in(a_begin + j1, num1));
        if (num2 != 0)b_gap.emplace_back(in(b_begin + j2, num2));

    }

    return std::make_tuple(std::move(a_gap), std::move(b_gap));

}

std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin, bool tag) {
    insert insertions1, insertions2;
    int index1 = a_begin, index2 = b_begin;
    int Dnum = 0;
    int Inum = 0;

    for (char c : cigar) {
        if (c == 'I') {
            index2++;
            Inum++;
        }
        else if (c == 'D') {
            index1++;
            Dnum++;
        }
        else {
            if (Dnum > 0) {
                insertions2.emplace_back(std::make_tuple(index2, Dnum));
                Dnum = 0;
            }
            else if (Inum > 0) {
                insertions1.emplace_back(std::make_tuple(index1, Inum));
                Inum = 0;
            }
            index1++;
            index2++;
        }
    }
    //最后可能为 D I
    if (Dnum > 0) {
        insertions2.emplace_back(std::make_tuple(index2, Dnum));
        Dnum = 0;
    }
    else if (Inum > 0) {
        insertions1.emplace_back(std::make_tuple(index1, Inum));
        Inum = 0;
    }
    if(tag)
        return make_tuple(insertions2, insertions1);
    return make_tuple(insertions1, insertions2);
}

std::tuple<insert, insert> parseCigar(const std::string& cigar, size_t a_begin, size_t b_begin) {
    insert insertions1, insertions2;
    int index1 = a_begin, index2 = b_begin;
    int Dnum = 0;
    int Inum = 0;

    for (char c : cigar) {
        if (c == 'I') {
            index2++;
            Inum++;
        }
        else if (c == 'D') {
            index1++;
            Dnum++;
        }
        else {
            if (Dnum > 0) {
                insertions2.emplace_back(std::make_tuple(index2, Dnum));
                Dnum = 0;
            }
            else if (Inum > 0) {
                insertions1.emplace_back(std::make_tuple(index1, Inum));
                Inum = 0;
            }
            index1++;
            index2++;
        }
    }
    //最后可能为 D I
    if (Dnum > 0) {
        insertions2.emplace_back(std::make_tuple(index2, Dnum));
        Dnum = 0;
    }
    else if (Inum > 0) {
        insertions1.emplace_back(std::make_tuple(index1, Inum));
        Inum = 0;
    }

    return make_tuple(insertions1, insertions2);
}


std::tuple<insert, insert> parseCigar(const std::string& cigar) {
    insert insertions1, insertions2;
    int index1 = 0, index2 = 0;
    int Dnum = 0;
    int Inum = 0;

    for (char c : cigar) {
        if (c == 'I') {
            index2++;
            Inum++;
        }
        else if (c == 'D') {
            index1++;
            Dnum++;
        }
        else {
            if (Dnum > 0) {
                insertions2.emplace_back(std::make_tuple(index2, Dnum));
                Dnum = 0;
            }
            else if (Inum > 0) {
                insertions1.emplace_back(std::make_tuple(index1, Inum));
                Inum = 0;
            }
            index1++;
            index2++;
        }
    }
    //最后可能为 D I
    if (Dnum > 0) {
        insertions2.emplace_back(std::make_tuple(index2, Dnum));
        Dnum = 0;
    }
    else if (Inum > 0) {
        insertions1.emplace_back(std::make_tuple(index1, Inum));
        Inum = 0;
    }

    return make_tuple(insertions1, insertions2);
}


// 在字符串中插入空插入
void insertGaps(std::string& str, const insert& insertions) {
    for (auto it = insertions.rbegin(); it != insertions.rend(); ++it) {
        int index = std::get<0>(*it);
        int num = std::get<1>(*it);
        str.insert(index, num, '-');
    }
}

std::tuple<insert, insert> mywfa(wfa::WFAlignerGapAffine& aligner, const std::vector<unsigned char>& sequence1, size_t a_begin, size_t a_end, const std::vector<unsigned char>& sequence2, int b_begin, int b_end) {
    int a_len = a_end - a_begin;
    int b_len = b_end - b_begin;
    insert a_gap;
    insert b_gap;
    if (a_len == 0)                             //a为0
    {
        a_gap.emplace_back(in(a_begin, b_len));
        return std::make_tuple(std::move(a_gap), std::move(b_gap));
    }
    else if (b_len == 0)                        //b为0
    {
        b_gap.emplace_back(in(b_begin, a_len));
        return std::make_tuple(std::move(a_gap), std::move(b_gap));
    }
    /*if (b_len > a_len)
    {
        aligner.alignEnd2End(&sequence2[b_begin], b_len, &sequence1[a_begin], a_len);
        return parseCigar(aligner.getAlignment(), b_begin, a_begin, 1);//1 交换
    }
    aligner.alignEnd2End(&sequence1[a_begin], a_len, &sequence2[b_begin], b_len);
    return parseCigar(aligner.getAlignment(), a_begin, b_begin, 0);
    
    aligner.alignEnd2End(&sequence1[a_begin], a_len, &sequence2[b_begin], b_len);
    return parseCigar(aligner.getAlignment(), a_begin, b_begin);
*/
    std::string A = std::string(sequence1.begin() + a_begin, sequence1.begin() + a_end);
    std::string B = std::string(sequence2.begin() + b_begin, sequence2.begin() + b_end);
    
    
    
    if (b_len > a_len)
    {
        aligner.alignEnd2End(B, A);
        return parseCigar(aligner.getAlignment(), b_begin, a_begin, 1);//1 交换
    }



    
    aligner.alignEnd2End(A, B);
    return parseCigar(aligner.getAlignment(), a_begin, b_begin);
}
