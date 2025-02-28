#pragma once

namespace seq_utils {

static inline constexpr unsigned char complementLookup[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', // 0-15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', // 16-31
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    '-', 'N', // 32-47
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', // 48-63
    'N', 'T', 'B', 'G', 'D', 'N', 'N', 'C', 'H', 'N', 'N', 'K', 'M', 'N',
    'N', 'N', // 64-79  (@,A-O)
    'N', 'N', 'Y', 'S', 'A', 'U', 'V', 'W', 'N', 'R', 'N', 'N', 'N', 'N',
    'N', 'N', // 80-95  (P-Z)
    'N', 't', 'b', 'g', 'd', 'n', 'n', 'c', 'h', 'n', 'n', 'k', 'm', 'n',
    'n', 'n', // 96-111  (`,a-o)
    'n', 'n', 'y', 's', 'a', 'u', 'v', 'w', 'n', 'r', 'n', 'N', 'N', 'N',
    'N', 'N', // 112-127 (p-z)
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', // 128-255 (extended ASCII)
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

static inline constexpr char getComplementCharacter(char nuc) {
  return complementLookup[static_cast<unsigned char>(nuc)];
}

} // namespace seq_utils