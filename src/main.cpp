#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__aarch64__) || defined(__arm64__)
#include <arm_neon.h>
#endif

namespace {

struct Sequence {
    std::string name;
    std::string bases;
};

struct Scoring {
    int match = 2;
    int mismatch = -1;
    int gap = -2;
};

struct Alignment {
    int score = 0;
    std::size_t reference_start = 0;
    std::size_t reference_end = 0;
    std::size_t read_start = 0;
    std::size_t read_end = 0;
};

struct Cell {
    int score = 0;
    std::size_t read_start = 0;
    std::size_t reference_start = 0;
};

struct Options {
    std::string reference_path;
    std::string reads_path;
    std::string mode = "rolling";
    Scoring scoring;
};

std::string trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<Sequence> read_fasta(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Ne mogu otvoriti FASTA fajl: " + path);
    }

    std::vector<Sequence> sequences;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            sequences.push_back({trim(line.substr(1)), {}});
            if (sequences.back().name.empty()) {
                throw std::runtime_error("FASTA zaglavlje ne sadrži ime sekvence");
            }
            continue;
        }
        if (sequences.empty()) {
            throw std::runtime_error("FASTA sekvenca mora početi zaglavljem '>'");
        }
        for (const unsigned char ch : line) {
            if (!std::isspace(ch)) {
                sequences.back().bases.push_back(
                    static_cast<char>(std::toupper(ch)));
            }
        }
    }

    if (sequences.empty()) {
        throw std::runtime_error("FASTA fajl ne sadrži sekvence: " + path);
    }
    for (const auto& sequence : sequences) {
        if (sequence.bases.empty()) {
            throw std::runtime_error("Prazna sekvenca u FASTA fajlu: " + sequence.name);
        }
    }
    return sequences;
}

int substitution_score(char read_base, char reference_base,
                       const Scoring& scoring) {
    return read_base == reference_base ? scoring.match : scoring.mismatch;
}

Cell choose_cell(const Cell& diagonal, const Cell& up, const Cell& left,
                 int diagonal_score, int up_score, int left_score,
                 std::size_t read_index, std::size_t reference_index) {
    Cell result;
    result.score = std::max({0, diagonal_score, up_score, left_score});
    if (result.score == 0) {
        // Coordinates are boundaries between bases. A positive diagonal
        // starting here will convert these to zero-based inclusive starts.
        result.read_start = read_index;
        result.reference_start = reference_index;
    } else if (diagonal_score == result.score) {
        if (diagonal.score > 0) {
            result.read_start = diagonal.read_start;
            result.reference_start = diagonal.reference_start;
        } else {
            result.read_start = read_index - 1;
            result.reference_start = reference_index - 1;
        }
    } else if (up_score == result.score) {
        result.read_start = up.read_start;
        result.reference_start = up.reference_start;
    } else {
        result.read_start = left.read_start;
        result.reference_start = left.reference_start;
    }
    return result;
}

void consider_best(const Cell& cell, std::size_t read_end,
                   std::size_t reference_end, Alignment& best) {
    if (cell.score > best.score) {
        best.score = cell.score;
        best.read_start = cell.read_start;
        best.reference_start = cell.reference_start;
        best.read_end = read_end;
        best.reference_end = reference_end;
    }
}

Alignment smith_waterman_full(const std::string& read,
                              const std::string& reference,
                              const Scoring& scoring,
                              double& dp_seconds) {
    const std::size_t rows = read.size() + 1;
    const std::size_t columns = reference.size() + 1;
    if (columns > std::numeric_limits<std::size_t>::max() / rows) {
        throw std::length_error("Matrica je prevelika");
    }
    std::vector<Cell> matrix(rows * columns);
    const auto at = [&matrix, columns](std::size_t i, std::size_t j) -> Cell& {
        return matrix[i * columns + j];
    };

    Alignment best;
    const auto dp_start = std::chrono::steady_clock::now();
    for (std::size_t i = 1; i < rows; ++i) {
        for (std::size_t j = 1; j < columns; ++j) {
            const Cell& diagonal = at(i - 1, j - 1);
            const Cell& up = at(i - 1, j);
            const Cell& left = at(i, j - 1);
            Cell& current = at(i, j);
            current = choose_cell(
                diagonal, up, left,
                diagonal.score + substitution_score(read[i - 1], reference[j - 1], scoring),
                up.score + scoring.gap,
                left.score + scoring.gap,
                i, j);
            consider_best(current, i, j, best);
        }
    }
    dp_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dp_start).count();
    return best;
}

Alignment smith_waterman_rolling(const std::string& read,
                                 const std::string& reference,
                                 const Scoring& scoring,
                                 double& dp_seconds) {
    std::vector<Cell> previous(reference.size() + 1);
    std::vector<Cell> current(reference.size() + 1);
    Alignment best;

    const auto dp_start = std::chrono::steady_clock::now();
    for (std::size_t i = 1; i <= read.size(); ++i) {
        current[0] = Cell{0, i, 0};
        for (std::size_t j = 1; j <= reference.size(); ++j) {
            const Cell& diagonal = previous[j - 1];
            const Cell& up = previous[j];
            const Cell& left = current[j - 1];
            current[j] = choose_cell(
                diagonal, up, left,
                diagonal.score + substitution_score(read[i - 1], reference[j - 1], scoring),
                up.score + scoring.gap,
                left.score + scoring.gap,
                i, j);
            consider_best(current[j], i, j, best);
        }
        std::swap(previous, current);
    }
    dp_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dp_start).count();
    return best;
}

Alignment smith_waterman_trace(const std::string& read,
                               const std::string& reference,
                               const Scoring& scoring,
                               double& dp_seconds) {
    const std::size_t rows = read.size() + 1;
    const std::size_t columns = reference.size() + 1;
    const std::size_t direction_bytes_per_row =
        reference.size() / 4 + (reference.size() % 4 != 0);
    if (direction_bytes_per_row >
        std::numeric_limits<std::size_t>::max() / read.size()) {
        throw std::length_error("Traceback matrica je prevelika");
    }
    const std::size_t direction_bytes = read.size() * direction_bytes_per_row;
    // Every packed byte is assigned during the DP pass, so avoid zeroing the
    // entire traceback buffer before writing it.
    std::unique_ptr<std::uint8_t[]> directions(new std::uint8_t[direction_bytes]);
    std::uint8_t* direction_data = directions.get();
    const char* read_bases = read.data();
    const char* reference_bases = reference.data();
    std::vector<int> previous(columns, 0), current(columns, 0);
    Alignment best;

    const auto dp_start = std::chrono::steady_clock::now();
    const auto process_cell = [&](std::size_t j, std::size_t i,
                                  const char read_base, int* previous_scores,
                                  int* current_scores) {
        const int substitution = read_base == reference_bases[j - 1]
            ? scoring.match : scoring.mismatch;
        const int diagonal_score = previous_scores[j - 1] + substitution;
        const int up_score = previous_scores[j] + scoring.gap;
        const int left_score = current_scores[j - 1] + scoring.gap;
        int score = diagonal_score;
        if (up_score > score) score = up_score;
        if (left_score > score) score = left_score;
        if (score < 0) score = 0;

        std::uint8_t direction = 0;
        if (score > 0) {
            if (diagonal_score == score) direction = 1;
            else if (up_score == score) direction = 2;
            else direction = 3;
        }
        current_scores[j] = score;
        if (score > best.score) {
            best.score = score;
            best.read_end = i;
            best.reference_end = j;
        }
        return direction;
    };

    for (std::size_t i = 1; i < rows; ++i) {
        int* previous_scores = previous.data();
        int* current_scores = current.data();
        const char read_base = read_bases[i - 1];
        current_scores[0] = 0;
        const std::size_t direction_row_offset = (i - 1) * direction_bytes_per_row;
        std::size_t j = 1;
        for (; j + 3 < columns; j += 4) {
            std::uint8_t packed_directions = process_cell(
                j, i, read_base, previous_scores, current_scores);
            packed_directions |= static_cast<std::uint8_t>(
                process_cell(j + 1, i, read_base, previous_scores, current_scores) << 2);
            packed_directions |= static_cast<std::uint8_t>(
                process_cell(j + 2, i, read_base, previous_scores, current_scores) << 4);
            packed_directions |= static_cast<std::uint8_t>(
                process_cell(j + 3, i, read_base, previous_scores, current_scores) << 6);
            direction_data[direction_row_offset + (j - 1) / 4] = packed_directions;
        }
        if (j < columns) {
            std::uint8_t packed_directions = 0;
            unsigned int shift = 0;
            const std::size_t first_column = j;
            for (; j < columns; ++j, shift += 2) {
                packed_directions |= static_cast<std::uint8_t>(
                    process_cell(j, i, read_base, previous_scores, current_scores) << shift);
            }
            direction_data[direction_row_offset + (first_column - 1) / 4] =
                packed_directions;
        }
        std::swap(previous, current);
    }

    std::size_t read_index = best.read_end;
    std::size_t reference_index = best.reference_end;
    while (read_index > 0 && reference_index > 0) {
        const std::size_t direction_row_offset =
            (read_index - 1) * direction_bytes_per_row;
        const std::size_t direction_index = reference_index - 1;
        const unsigned int shift = static_cast<unsigned int>((direction_index % 4) * 2);
        const std::uint8_t direction =
            static_cast<std::uint8_t>((direction_data[
                direction_row_offset + direction_index / 4] >> shift) & 0x3);
        if (direction == 0) break;
        if (direction == 1) {
            --read_index;
            --reference_index;
        } else if (direction == 2) {
            --read_index;
        } else {
            --reference_index;
        }
    }
    best.read_start = read_index;
    best.reference_start = reference_index;
    dp_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dp_start).count();
    return best;
}

void smith_waterman_simd_batch(const std::vector<Sequence>& reads,
                               const std::string& reference,
                               const Scoring& scoring,
                               std::vector<Alignment>& alignments,
                               double& dp_seconds) {
#if defined(__aarch64__) || defined(__arm64__)
    constexpr std::size_t lanes = 4;
    const int32x4_t zero = vdupq_n_s32(0);
    const std::size_t columns = reference.size() + 1;

    for (std::size_t base = 0; base + lanes <= reads.size(); base += lanes) {
        std::size_t max_read_length = 0;
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            max_read_length = std::max(max_read_length, reads[base + lane].bases.size());
        }

        std::vector<int32x4_t> previous(columns, zero), current(columns, zero);
        std::vector<int32x4_t> previous_read_start(columns, zero);
        std::vector<int32x4_t> previous_reference_start(columns, zero);
        std::vector<int32x4_t> current_read_start(columns, zero);
        std::vector<int32x4_t> current_reference_start(columns, zero);
        int32x4_t best_score = zero;
        int32x4_t best_read_start = zero;
        int32x4_t best_reference_start = zero;
        int32x4_t best_read_end = zero;
        int32x4_t best_reference_end = zero;

        const auto dp_start = std::chrono::steady_clock::now();
        for (std::size_t i = 1; i <= max_read_length; ++i) {
            current[0] = zero;
            current_read_start[0] = vdupq_n_s32(static_cast<int32_t>(i));
            current_reference_start[0] = zero;

            for (std::size_t j = 1; j < columns; ++j) {
                int32_t substitutions[lanes];
                int32_t active_lanes[lanes];
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    const auto& read = reads[base + lane].bases;
                    active_lanes[lane] = i <= read.size() ? -1 : 0;
                    substitutions[lane] = i <= read.size()
                        ? substitution_score(read[i - 1], reference[j - 1], scoring)
                        : 0;
                }

                const int32x4_t diagonal = previous[j - 1];
                const int32x4_t up = previous[j];
                const int32x4_t left = current[j - 1];
                const int32x4_t diagonal_score = vaddq_s32(
                    diagonal, vld1q_s32(substitutions));
                const int32x4_t up_score = vaddq_s32(up, vdupq_n_s32(scoring.gap));
                const int32x4_t left_score = vaddq_s32(left, vdupq_n_s32(scoring.gap));
                int32x4_t cell_score = vmaxq_s32(zero,
                    vmaxq_s32(diagonal_score, vmaxq_s32(up_score, left_score)));
                const int32x4_t active = vld1q_s32(active_lanes);
                cell_score = vandq_s32(cell_score, active);

                const int32x4_t positive = vcgtq_s32(cell_score, zero);
                const int32x4_t choose_diagonal = vandq_s32(
                    vceqq_s32(diagonal_score, cell_score), positive);
                const int32x4_t choose_up = vandq_s32(
                    vandq_s32(vceqq_s32(up_score, cell_score), positive),
                    vmvnq_s32(choose_diagonal));
                const int32x4_t choose_left = vandq_s32(
                    vandq_s32(vceqq_s32(left_score, cell_score), positive),
                    vmvnq_s32(vorrq_s32(choose_diagonal, choose_up)));

                int32x4_t read_start = vdupq_n_s32(static_cast<int32_t>(i));
                int32x4_t reference_start = vdupq_n_s32(static_cast<int32_t>(j));
                const int32x4_t diagonal_read_start = vbslq_s32(
                    vcgtq_s32(diagonal, zero), previous_read_start[j - 1],
                    vdupq_n_s32(static_cast<int32_t>(i - 1)));
                const int32x4_t diagonal_reference_start = vbslq_s32(
                    vcgtq_s32(diagonal, zero), previous_reference_start[j - 1],
                    vdupq_n_s32(static_cast<int32_t>(j - 1)));
                read_start = vbslq_s32(choose_diagonal, diagonal_read_start, read_start);
                reference_start = vbslq_s32(
                    choose_diagonal, diagonal_reference_start, reference_start);
                read_start = vbslq_s32(
                    choose_up, previous_read_start[j], read_start);
                reference_start = vbslq_s32(
                    choose_up, previous_reference_start[j], reference_start);
                read_start = vbslq_s32(choose_left, current_read_start[j - 1], read_start);
                reference_start = vbslq_s32(
                    choose_left, current_reference_start[j - 1], reference_start);

                current[j] = cell_score;
                current_read_start[j] = read_start;
                current_reference_start[j] = reference_start;

                const int32x4_t new_best = vcgtq_s32(cell_score, best_score);
                best_score = vbslq_s32(new_best, cell_score, best_score);
                best_read_start = vbslq_s32(new_best, read_start, best_read_start);
                best_reference_start = vbslq_s32(
                    new_best, reference_start, best_reference_start);
                best_read_end = vbslq_s32(
                    new_best, vdupq_n_s32(static_cast<int32_t>(i)), best_read_end);
                best_reference_end = vbslq_s32(
                    new_best, vdupq_n_s32(static_cast<int32_t>(j)), best_reference_end);
            }
            previous.swap(current);
            previous_read_start.swap(current_read_start);
            previous_reference_start.swap(current_reference_start);
        }
        dp_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - dp_start).count();

        int32_t scores[lanes], read_starts[lanes], reference_starts[lanes];
        int32_t read_ends[lanes], reference_ends[lanes];
        vst1q_s32(scores, best_score);
        vst1q_s32(read_starts, best_read_start);
        vst1q_s32(reference_starts, best_reference_start);
        vst1q_s32(read_ends, best_read_end);
        vst1q_s32(reference_ends, best_reference_end);
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            alignments[base + lane] = Alignment{
                scores[lane],
                static_cast<std::size_t>(reference_starts[lane]),
                static_cast<std::size_t>(reference_ends[lane]),
                static_cast<std::size_t>(read_starts[lane]),
                static_cast<std::size_t>(read_ends[lane])};
        }
    }

    // Process fewer than four remaining reads with the scalar reference.
    const std::size_t remainder = reads.size() - reads.size() % lanes;
    for (std::size_t i = remainder; i < reads.size(); ++i) {
        alignments[i] = smith_waterman_rolling(
            reads[i].bases, reference, scoring, dp_seconds);
    }
#else
    (void)reads;
    (void)reference;
    (void)scoring;
    (void)alignments;
    (void)dp_seconds;
    throw std::runtime_error("--mode simd je dostupan samo na ARM64 procesorima sa NEON podrškom");
#endif
}

std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\n\r") == std::string::npos) {
        return value;
    }
    std::string escaped = "\"";
    for (const char ch : value) {
        if (ch == '\"') {
            escaped += "\"\"";
        } else {
            escaped += ch;
        }
    }
    escaped += '\"';
    return escaped;
}

void print_usage(std::ostream& out) {
    out << "Upotreba: sw_mapper --reference ref.fa --reads reads.fa [opcije]\n"
           "Opcije:\n"
           "  --mode full|rolling|trace|simd   sekvencijalno, compact trace ili batch SIMD\n"
           "  --match N                    bodovi za podudaranje (default 2)\n"
           "  --mismatch N                 kazna za nepodudaranje (default -1)\n"
           "  --gap N                      linearna kazna za gap (default -2)\n"
           "  --help                       prikaz ove poruke\n";
}

int parse_int(const std::string& value, const std::string& option) {
    std::size_t parsed = 0;
    int result = 0;
    try {
        result = std::stoi(value, &parsed);
    } catch (const std::exception&) {
        throw std::runtime_error("Neispravan broj za " + option + ": " + value);
    }
    if (parsed != value.size()) {
        throw std::runtime_error("Neispravan broj za " + option + ": " + value);
    }
    return result;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(std::cout);
            std::exit(0);
        }
        if (i + 1 >= argc) {
            throw std::runtime_error("Nedostaje vrijednost za opciju " + arg);
        }
        const std::string value = argv[++i];
        if (arg == "--reference") options.reference_path = value;
        else if (arg == "--reads") options.reads_path = value;
        else if (arg == "--mode") options.mode = value;
        else if (arg == "--match") options.scoring.match = parse_int(value, arg);
        else if (arg == "--mismatch") options.scoring.mismatch = parse_int(value, arg);
        else if (arg == "--gap") options.scoring.gap = parse_int(value, arg);
        else throw std::runtime_error("Nepoznata opcija: " + arg);
    }

    if (options.reference_path.empty() || options.reads_path.empty()) {
        throw std::runtime_error("Opcije --reference i --reads su obavezne");
    }
    if (options.mode != "full" && options.mode != "rolling" &&
        options.mode != "trace" && options.mode != "simd") {
        throw std::runtime_error("--mode mora biti full, rolling, trace ili simd");
    }
    if (options.scoring.match <= 0 || options.scoring.mismatch > 0 ||
        options.scoring.gap >= 0) {
        throw std::runtime_error("Očekuje se match > 0, mismatch <= 0 i gap < 0");
    }
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const auto reference_records = read_fasta(options.reference_path);
        if (reference_records.size() != 1) {
            throw std::runtime_error("Referentni FASTA mora sadržavati tačno jednu sekvencu");
        }
        const auto reads = read_fasta(options.reads_path);
        const std::string& reference = reference_records.front().bases;
        std::vector<Alignment> alignments(reads.size());
        double dp_seconds = 0.0;

        const auto start = std::chrono::steady_clock::now();
        if (options.mode == "simd") {
            smith_waterman_simd_batch(
                reads, reference, options.scoring, alignments, dp_seconds);
        } else {
            for (std::size_t i = 0; i < reads.size(); ++i) {
                if (options.mode == "full") {
                    alignments[i] = smith_waterman_full(
                        reads[i].bases, reference, options.scoring, dp_seconds);
                } else if (options.mode == "trace") {
                    alignments[i] = smith_waterman_trace(
                        reads[i].bases, reference, options.scoring, dp_seconds);
                } else {
                    alignments[i] = smith_waterman_rolling(
                        reads[i].bases, reference, options.scoring, dp_seconds);
                }
            }
        }
        const auto end = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(end - start).count();

        std::cout << "read_id,score,reference_start,reference_end,read_start,read_end\n";
        for (std::size_t i = 0; i < reads.size(); ++i) {
            const auto& a = alignments[i];
            std::cout << csv_escape(reads[i].name) << ',' << a.score << ','
                      << a.reference_start << ',' << a.reference_end << ','
                      << a.read_start << ',' << a.read_end << '\n';
        }

        long double cells = 0;
        for (const auto& read : reads) {
            cells += static_cast<long double>(reference.size()) * read.bases.size();
        }
        const double gcups = seconds > 0.0
            ? static_cast<double>(cells / seconds / 1'000'000'000.0L) : 0.0;
        std::cerr << std::fixed << std::setprecision(6)
                  << "mode=" << options.mode
                  << " reads=" << reads.size()
                  << " reference_length=" << reference.size()
                  << " seconds=" << seconds
                  << " dp_seconds=" << dp_seconds
                  << " GCUPS=" << gcups;
        std::cerr << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Greška: " << error.what() << '\n';
        print_usage(std::cerr);
        return 1;
    }
    return 0;
}
