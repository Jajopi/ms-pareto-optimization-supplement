#pragma once

#include <vector>
#include <algorithm>
#include <string>

#include "../kmers.h"
#include "../parser.h"
#include "objective.h"

/// Sort k-mers with key_bits significant bits
template <typename kmer_t>
void ParallelSort(typename std::vector<kmer_t>::iterator begin, typename std::vector<kmer_t>::iterator end, int key_bits){
    std::sort(begin, end);
}

/// Merge sorted runs [0, middle) and [middle, size) of the vector (backward, only the second run is copied)
template <typename kmer_t>
void MergeSortedRuns(std::vector<kmer_t> &kMers, size_t middle){
    std::vector<kmer_t> buffer(kMers.begin() + middle, kMers.end());
    size_t i = middle, j = buffer.size(), out = kMers.size();
    while (j > 0){
        if (i > 0 && buffer[j - 1] < kMers[i - 1]) kMers[--out] = kMers[--i];
        else kMers[--out] = buffer[--j];
    }
}

/// Add reverse complements to sorted canonical k-mers, keeping them sorted (self-complementary k-mers end up twice)
template <typename kmer_t>
void AddComplementsSorted(std::vector<kmer_t> &kMers, size_k_max k){
    size_t n = kMers.size();
    if (!std::is_sorted(kMers.begin(), kMers.end())) ParallelSort<kmer_t>(kMers.begin(), kMers.end(), 2 * k);
    kMers.resize(2 * n);
    for (size_t i = 0; i < n; ++i) kMers[n + i] = ReverseComplement(kMers[i], k);
    ParallelSort<kmer_t>(kMers.begin() + n, kMers.end(), 2 * k);
    MergeSortedRuns(kMers, n);
}

/// Read distinct k-mers of a FASTA file into a sorted vector (canonical ones if complements is true).
/// Occurrences are buffered and periodically sorted, deduplicated and merged, a small table of recently seen k-mers
/// filters most repetitions (e.g. overlaps of unitigs). The result has capacity for the reverse complements.
template <typename kmer_t>
std::vector<kmer_t> ReadKMersSorted(std::string &path, int k, bool complements){
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    std::vector<kmer_t> kMers;
    size_t sorted_size = 0;
    const size_t MIN_BATCH = size_t(1) << 22;

    /// The filter grows with the number of distinct k-mers, so that small inputs stay small
    const int MIN_FILTER_BITS = 16, MAX_FILTER_BITS = 22;
    const kmer_t NO_KMER = ((kmer_t) 1) << (2 * k); /// Not a valid k-mer (k-mer types always have more than 2k bits)
    int filter_bits = MIN_FILTER_BITS;
    std::vector<kmer_t> recently_seen(size_t(1) << filter_bits, NO_KMER);
    auto compact = [&](){
        ParallelSort<kmer_t>(kMers.begin() + sorted_size, kMers.end(), 2 * k);
        kMers.erase(std::unique(kMers.begin() + sorted_size, kMers.end()), kMers.end());
        MergeSortedRuns(kMers, sorted_size);
        kMers.erase(std::unique(kMers.begin(), kMers.end()), kMers.end());
        sorted_size = kMers.size();
        int wanted_bits = MIN_FILTER_BITS;
        while (wanted_bits < MAX_FILTER_BITS && (size_t(1) << wanted_bits) < sorted_size) ++wanted_bits;
        if (wanted_bits > filter_bits && !recently_seen.empty()){
            filter_bits = wanted_bits;
            recently_seen.assign(size_t(1) << filter_bits, NO_KMER);
        }
        size_t wanted = sorted_size + std::max(MIN_BATCH, sorted_size / 2);
        if (kMers.capacity() < wanted) kMers.reserve(wanted);
    };
    kMers.reserve(MIN_BATCH);

    kmer_t mask = (((kmer_t) 1) << (2 * k)) - 1;
    kmer_t shift = 2 * (k - 1);
    while (kseq_read(seq) >= 0){
        int64_t currentLength = 0;
        kmer_t currentKMer = 0, reverseComplement = 0;
        for (size_t i = 0; i < seq->seq.l; ++i){
            auto data = nucleotideToInt[(uint8_t)seq->seq.s[i]];
            if (data >= 4){ /// Restart if "N"-like nucleotide
                currentKMer = reverseComplement = 0;
                currentLength = 0;
                continue;
            }
            currentKMer = ((currentKMer << 2) | data) & mask;
            reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
            if (++currentLength < k) continue;
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            size_t slot = (uint64_t(canonical) * 0x9E3779B97F4A7C15ULL) >> (64 - filter_bits);
            if (recently_seen[slot] == canonical) continue;
            recently_seen[slot] = canonical;
            kMers.push_back(canonical);
            if (kMers.size() == kMers.capacity()) compact();
        }
    }
    std::vector<kmer_t>().swap(recently_seen);
    compact();
    kseq_destroy(seq);
    gzclose(fp);

    std::vector<kmer_t> result;
    result.reserve(complements ? 2 * kMers.size() : kMers.size());
    result.assign(kMers.begin(), kMers.end());
    return result;
}
