// Find many orfs faster by putting all the data into Cpp directly
// Adapted from code by Haakon Tjeldnes et al., originally part of ORFik
// License: MIT



#include <Rcpp.h>
#include <algorithm>
#include <locale>
#include <vector>
#include "findORFsHelpers.h"

using vi = std::vector<int>;
using namespace Rcpp;

Function GRangesC("GRanges", Environment::namespace_env("GenomicRanges"));
Function IRangesC("IRanges", Environment::namespace_env("IRanges"));
Function namesC("names", Environment::namespace_env("base"));

char complement(char n) {
    switch (n) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
    case 'N': return 'N';
    case 'a': return 't';
    case 't': return 'a';
    case 'g': return 'c';
    case 'c': return 'g';
    default: return ' ';
    }
}

// [[Rcpp::export]]
S4 findORFsFastaCpp(CharacterVector &fastaSeqs,
                    std::string startCodon,
                    std::string stopCodon,
                    int minimumLength,
                    bool isCircular,
                    bool plusStrandOnly) {
    vi all_orfs;
    std::vector<std::string> Seqnames;
    std::vector<int> strands;
    CharacterVector headers = namesC(fastaSeqs);

    for (int i = 0; i < fastaSeqs.size(); i++) {
        std::string fastaSeq = static_cast<std::string>(fastaSeqs[i]);
        std::string header = static_cast<std::string>(headers[i]);
        header = header.substr(0, header.find(' '));

        int length = fastaSeq.length();
        int chromoLength = length + 1;

        // Handle lowercase sequences
        std::locale loc;
        if (std::islower(fastaSeq.at(0), loc)) {
            std::transform(startCodon.begin(), startCodon.end(), startCodon.begin(), ::tolower);
            std::transform(stopCodon.begin(), stopCodon.end(), stopCodon.begin(), ::tolower);
        }

        // Forward strand ORFs
        vi ORFdef = orfs_as_vector(fastaSeq, startCodon, stopCodon, minimumLength);
        all_orfs.insert(all_orfs.end(), ORFdef.begin(), ORFdef.end());
        Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, header);
        strands.insert(strands.end(), ORFdef.size() / 2, 1);

        // Circular genome handling for forward strand
        if (isCircular) {
            std::string startStopBoundary = fastaSeq + fastaSeq; // full wrap
            vi ORFdefBoundary = orfs_as_vector(startStopBoundary, startCodon, stopCodon, minimumLength);
            vi ORFdefOverlap;

            for (size_t j = 0; j < ORFdefBoundary.size() / 2; j++) {
                int startPos = ORFdefBoundary[2 * j];
                int endPos = ORFdefBoundary[2 * j + 1];

                if (endPos >= length) { // crosses boundary
                    int startWrapped = startPos % length;
                    int endWrapped = endPos % length;
                    if (endWrapped < startWrapped) endWrapped += length; // keep order
                    ORFdefOverlap.push_back(startWrapped);
                    ORFdefOverlap.push_back(endWrapped);
                }
            }

            all_orfs.insert(all_orfs.end(), ORFdefOverlap.begin(), ORFdefOverlap.end());
            Seqnames.insert(Seqnames.end(), ORFdefOverlap.size() / 2, header);
            strands.insert(strands.end(), ORFdefOverlap.size() / 2, 1);
        }

        // Reverse complement strand if allowed
        if (!plusStrandOnly) {
            std::reverse(fastaSeq.begin(), fastaSeq.end());
            std::transform(fastaSeq.begin(), fastaSeq.end(), fastaSeq.begin(), complement);

            ORFdef = orfs_as_vector(fastaSeq, startCodon, stopCodon, minimumLength);
            for (size_t j = 0; j < ORFdef.size(); j++) ORFdef[j] = chromoLength - ORFdef[j];
            all_orfs.insert(all_orfs.end(), ORFdef.rbegin(), ORFdef.rend());
            Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, header);
            strands.insert(strands.end(), ORFdef.size() / 2, -1);

            if (isCircular) {
                std::string startStopBoundary = fastaSeq + fastaSeq;
                vi ORFdefBoundary = orfs_as_vector(startStopBoundary, startCodon, stopCodon, minimumLength);
                vi ORFdefOverlapMin;

                for (size_t j = 0; j < ORFdefBoundary.size() / 2; j++) {
                    int startPos = ORFdefBoundary[2 * j];
                    int endPos = ORFdefBoundary[2 * j + 1];

                    if (endPos >= length) {
                        int startWrapped = startPos % length;
                        int endWrapped = endPos % length;
                        if (endWrapped < startWrapped) endWrapped += length;
                        ORFdefOverlapMin.push_back(startWrapped);
                        ORFdefOverlapMin.push_back(endWrapped);
                    }
                }

                for (size_t j = 0; j < ORFdefOverlapMin.size(); j++) {
                    ORFdefOverlapMin[j] = chromoLength - ORFdefOverlapMin[j];
                }

                all_orfs.insert(all_orfs.end(), ORFdefOverlapMin.rbegin(), ORFdefOverlapMin.rend());
                Seqnames.insert(Seqnames.end(), ORFdefOverlapMin.size() / 2, header);
                strands.insert(strands.end(), ORFdefOverlapMin.size() / 2, -1);
            }
        }
    }

    // De-interlace ORFs
    std::vector<vi> result_value(2);
    result_value[0].resize(all_orfs.size() / 2);
    result_value[1].resize(all_orfs.size() / 2);
    for (size_t i = 0; i < all_orfs.size() / 2; i++) {
        result_value[0][i] = all_orfs[2 * i];
        result_value[1][i] = all_orfs[2 * i + 1];
    }

    return GRangesC(Seqnames, IRangesC(wrap(result_value[0]), wrap(result_value[1])), strands);
}
