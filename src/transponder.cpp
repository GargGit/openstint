#include "transponder.hpp"

#include <algorithm>
#include <bit>
#include <string>

extern "C" {
#include <fec.h>
}
#include <liquid/liquid.h>

static void *viterbi_decoder;
static crc_scheme crc8_scheme = LIQUID_CRC_8;

void init_transponders() {
    viterbi_decoder = create_viterbi29(32);
}

int decode_openstint(const uint8_t *softbits, uint32_t *transponder_id) {
    uint8_t decoded[4];

    init_viterbi29(viterbi_decoder, 0);
    update_viterbi29_blk(viterbi_decoder, const_cast<uint8_t*>(softbits), 32+8); // khm...
    chainback_viterbi29(viterbi_decoder, decoded, 32, 0);
    
    *transponder_id = (static_cast<uint32_t>(decoded[0]) << 16) | (static_cast<uint32_t>(decoded[1]) << 8) | static_cast<uint32_t>(decoded[2]);
    return crc_validate_message(crc8_scheme, decoded, 3, decoded[3]);
}

// Split the 32 bit message into the transponder id and the status code.
//
// Example: your transponder id is "1234567", or "0b00010010_11010110_10000111".
// To get the pre-encoded message, reverse the binary word and break into
// chunks of 3. You'll end up with 24/3=8 chunks:
// 111 000 010 110 101 101 001 000
// Suffix each chunk with a bit from a "status code",
// in my RC4 hybrid transponder, it was 00000101
// 1110 0000 0100 1100 1010 1011 0010 0001
static void split_rc3_message(uint32_t message, uint32_t *transponder_id, uint8_t *status_code) {
    uint32_t tid = 0; // transponder_id
    uint8_t status = 0; // status code
    for (int i=0; i<32; i++) {
        uint32_t bitmask = (1 << i);
        uint32_t bit = (message & bitmask) ? 1 : 0;

        if (i % 4 != 0) { // every 4th is status bit
            tid = (tid << 1) | bit;
        } else {
            status = (status << 1) | bit;
        }
    }
    *transponder_id = tid;
    *status_code = status;
}

// Strict RC3 decode: read the message off the symbols and reject anything that
// is not an exact codeword. This is the fast path; it decodes every frame that
// arrived without a single symbol error.
//
// RC3 use a K=24, r=1/2 convolutional encoder with polynoms 0xEEC20F and 0xEEC20D.
// Note, the generating polynoms differ only in one bit (0xf=1111 vs 0xd=1101).
// This reduces the complexity to:
// bit0 = parity(SHREG & 0xEEC20C) ^ SHREG[1] ^ SHREG[0]
// bit1 = parity(SHREG & 0xEEC20C) ^            SHREG[0]
// As such, we know:
// bit0 ^ bit1 = SHREG[1]
// SHREG[0] = bit1 ^ parity(SHREG & 0xEEC20C)
// so every pair carries one redundant bit: bit0 ^ bit1 must match SHREG[1].
// Together with the zero trail that checks all 48 redundant bits of the frame.
static bool decode_rc3_strict(const uint8_t *softbits, uint32_t *message) {
    uint64_t shreg = 0; // shift register

    // Before encoding, the 24 bit transponder id is scrambled with extra 8 bits,
    // resulting in 32 bits. This is further appended with 0x00 so viterbi-decoder (???)
    // can process it. The rate=1/2 encoder generates 2x40=80 bits in total.
    int prev_sym = 0; // differential decoder
    for (int i=0; i<80; i+=2) {
        int p = std::popcount(shreg & 0xEEC20C) % 2; // parity bit from SHREG

        // differential-BPSK is decoded here as well (sym ^ prev_sym magic):
        int sym = (softbits[i+0] > 127) ? 1 : 0;
        int b0 = sym ^ prev_sym;
        prev_sym = (softbits[i+1] > 127) ? 1 : 0;
        int b1 = prev_sym ^ sym; // decode softbit bit1

        int shreg1 = (shreg & 2) ? 1 : 0; // shift register last-1 bit
        if ((b0 ^ b1) != shreg1) { // the redundant bit disagrees: not a codeword
            return false;
        }
        shreg |= p ^ b1;
        shreg <<= 1;
    }

    shreg >>= 1;
    *message = static_cast<uint32_t>((shreg>>8) & 0xffffffff);

    // the last byte must be zero (tail==0 error check)
    return (shreg & 0xff) == 0;
}

// Beam (M-algorithm) decoder: search for the codeword the received symbols are
// most likely to be, instead of reading bits off them and hoping they are right.
// Worth roughly 3 dB over decode_rc3_strict() at a 90% decode rate.
//
// A Viterbi decoder is out of reach here -- the K=24 shift register has 8M
// states -- so the search keeps only the RC3_BEAM_PATHS best paths per step.
// Paths are scored on the soft symbol values, which is where most of the gain
// over a hard-decision search comes from.
//
// Every candidate is also tracked in hard symbol errors, and a path that is
// already beyond RC3_BEAM_MAX_DIST is dropped: a longer path can only differ in
// more symbols, so this cannot discard a codeword we would have accepted. It
// makes the common case cheap -- a preamble triggered by noise runs out of
// paths within a few steps instead of walking the whole trellis.
//
// The accept radius trades sensitivity against wrongly decoded ids
// increase RC3_BEAM_MAX_DIST by +1, runtime increase by ca x3

static constexpr int RC3_BEAM_PATHS = 64;   // paths kept per step
static constexpr int RC3_BEAM_MAX_DIST = 4; // accept a codeword within this many symbols

struct Rc3Path {
    uint64_t history; // message bits fed so far, newest in bit 0
    int last;         // last channel symbol emitted (differential encoder state)
    int mismatch;     // symbols differing from the received ones so far
    int metric;       // soft agreement with the received symbols, higher is better
};

static bool decode_rc3_beam(const uint8_t *softbits, uint32_t *message) {
    int hard[80];   // hard decision per symbol
    int weight[80]; // how much the symbol is worth, 0 when the symbol is a coin flip
    for (int i=0; i<80; i++) {
        int w = static_cast<int>(softbits[i]) - 128; // >=0 is a one, same as (softbits[i] > 127)
        hard[i] = (w >= 0) ? 1 : 0;
        weight[i] = w;
    }

    Rc3Path paths[2*RC3_BEAM_PATHS], next[2*RC3_BEAM_PATHS];
    int path_count = 1;
    paths[0] = {0, 0, 0, 0}; // the differential decoder's reference symbol is 0

    for (int step=0; step<40; step++) {
        const int branches = (step < 32) ? 2 : 1; // the trailing 8 bits are zero
        int next_count = 0;
        for (int i=0; i<path_count; i++) {
            const Rc3Path &p = paths[i];
            for (int u=0; u<branches; u++) {
                uint64_t history = (p.history << 1) | static_cast<uint64_t>(u);
                uint64_t shreg = history & 0xFFFFFF;
                // encode this candidate bit, differentially encoded as transmitted
                int s0 = p.last ^ (std::popcount(shreg & 0xEEC20F) % 2);
                int s1 = s0 ^ (std::popcount(shreg & 0xEEC20D) % 2);

                int mismatch = p.mismatch + (s0 != hard[2*step]) + (s1 != hard[2*step+1]);
                if (mismatch > RC3_BEAM_MAX_DIST) {
                    continue; // already too far off, and it cannot get closer
                }
                int metric = p.metric
                           + (s0 ? weight[2*step] : -weight[2*step])
                           + (s1 ? weight[2*step+1] : -weight[2*step+1]);
                next[next_count++] = {history, s1, mismatch, metric};
            }
        }
        if (next_count == 0) { // nothing plausible left
            return false;
        }
        if (next_count > RC3_BEAM_PATHS) { // keep the best paths only
            std::nth_element(next, next+RC3_BEAM_PATHS, next+next_count,
                [](const Rc3Path &a, const Rc3Path &b) { return a.metric > b.metric; });
            next_count = RC3_BEAM_PATHS;
        }
        std::copy(next, next+next_count, paths);
        path_count = next_count;
    }

    const Rc3Path *best = std::max_element(paths, paths+path_count,
        [](const Rc3Path &a, const Rc3Path &b) { return a.metric < b.metric; });

    *message = static_cast<uint32_t>(best->history >> 8); // drop the zero tail
    return true;
}

int decode_rc3(const uint8_t *softbits, uint32_t *transponder_id, uint8_t *status_code) {
    uint32_t message;
    if (!decode_rc3_strict(softbits, &message) && !decode_rc3_beam(softbits, &message)) {
        return 0;
    }

    split_rc3_message(message, transponder_id, status_code);
    return 1;
}

int decode_vostok(const uint8_t *softbits, uint32_t *transponder_id) {
    // Vostok transponders reuse the RC3 preamble, but the payload is not
    // convolutionally encoded at all. Once the differential-BPSK is undone, an
    // 80 bit frame is plain: a fixed header, the transponder id sent twice, and
    // two running XOR checksums.
    //
    //   byte:   0    1  |  2  3  4  |  5  |  6  7  8  |  9
    //          0x63 0x1A|    id     | c1  |    id     | c2
    //
    // The id is a 24 bit big-endian word, so 0x4A6151 is transponder 4874577.
    // Each checksum is the XOR of every byte before it:
    //
    //   c1 = 0x63 ^ 0x1A ^ id[0] ^ id[1] ^ id[2]   (XOR of bytes 0..4)
    //   c2 = c1 ^ id[0] ^ id[1] ^ id[2] = c1 ^ 0x79  (XOR of bytes 0..8)
    //
    // so c2 is fully determined by c1, and c1 ^ c2 == 0x79 in every frame.

    // differential-decode into 10 bytes, MSB first. prev=0 as in decode_rc3.
    uint8_t frame[10] = {0};
    int prev = 0;
    for (int i = 0; i < 80; i++) {
        int raw = (softbits[i] > 127) ? 1 : 0;
        frame[i / 8] |= static_cast<uint8_t>((raw ^ prev) << (7 - (i % 8)));
        prev = raw;
    }

    if (frame[0] != 0x63 || frame[1] != 0x1a) { // fixed header
        return 0;
    }
    // both copies of the transponder id must agree
    if (frame[2] != frame[6] || frame[3] != frame[7] || frame[4] != frame[8]) {
        return 0;
    }
    uint8_t c1 = 0;
    for (int i = 0; i < 5; i++) {
        c1 ^= frame[i];
    }
    if (c1 != frame[5] || static_cast<uint8_t>(c1 ^ 0x79) != frame[9]) {
        return 0;
    }

    *transponder_id = (static_cast<uint32_t>(frame[2]) << 16)
                    | (static_cast<uint32_t>(frame[3]) << 8)
                    | static_cast<uint32_t>(frame[4]);

    return 1;
}

void AmbRcBlacklist::process(uint64_t timestamp, uint8_t status_code, uint32_t transponder_id) {
    // not a candidate status/validation message:
    if ((status_code & 0xf8) != 0xf8) return; // not an AmbRc message
    if ((status_code & 0x07) == 0) return; // no counter set

    uint8_t msb8 = static_cast<uint8_t>((transponder_id >> 16) & 0xff);
    if (status_code != 0xff) { // status/validation message for sure
        msb8_timestamps[msb8] = timestamp;
    } else { // might be a status/validation message
        auto it = msb8_timestamps.find(msb8);
        if (it != msb8_timestamps.end() && timestamp <= (it->second + 250000ul)) {
            banned_transponders.insert(transponder_id);
        }
    }
}

bool AmbRcBlacklist::check_banned(uint32_t transponder_id) const {
    return banned_transponders.contains(transponder_id);
}
