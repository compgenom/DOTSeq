// Adapted from code by Haakon Tjeldnes et al., originally part of ORFik
// License: MIT

#ifndef PKG_findORFs_H
#define PKG_findORFs_H

std::vector<int> orfs_as_vector(
    std::string& main_string,
    std::string s, std::string e,
    int minimumLength);

#endif
