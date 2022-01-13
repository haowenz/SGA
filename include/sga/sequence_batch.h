#ifndef SGA_SEQUENCEBATCH_H_
#define SGA_SEQUENCEBATCH_H_

#include <assert.h>

#include <string>
#include <vector>

#include "kseq.h"
#include "sequence.h"
#include "utils.h"

KSEQ_INIT(int, read)

namespace sga {

class SequenceBatch {
 public:
  SequenceBatch(const uint32_t max_batch_size) {
    max_batch_size_ = max_batch_size;
    num_loaded_sequences_ = 0;
    sequence_batch_.reserve(max_batch_size_);
    // Construct once and use update methods when loading each batch
    sequence_batch_.assign(max_batch_size_, sga::Sequence(0, "", ""));
    sequence_file_ = nullptr;
    sequence_kseq_ = nullptr;
  }

  // Only string and vector, thus nothing to free
  ~SequenceBatch() {}

  uint32_t GetMaxBatchSize() const { return max_batch_size_; }

  uint32_t GetNumLoadedSequences() const { return num_loaded_sequences_; }

  std::vector<sga::Sequence> const &GetSequenceBatch() const {
    return sequence_batch_;
  }

  sga::Sequence const &GetSequence(uint32_t sequence_index) const {
    return sequence_batch_.at(sequence_index);
  }

  void InitializeLoading(const std::string &sequence_file_path) {
    sequence_file_path_ = sequence_file_path;
    sequence_file_ = fopen(sequence_file_path_.c_str(), "r");
    assert(sequence_file_ != NULL);
    sequence_kseq_ = kseq_init(fileno(sequence_file_));
  }

  // Return the number of sequences loaded into the batch
  // and return 0 if there is no more sequences
  uint32_t LoadBatch() {
    double real_start_time = sga::GetRealTime();
    num_loaded_sequences_ = 0;
    for (uint32_t sequence_index = 0; sequence_index < max_batch_size_;
         ++sequence_index) {
      int length = kseq_read(sequence_kseq_);
      // Skip the sequences of length 0
      if (length == 0)
        continue;
      else if (length > 0) {
        if (sequence_kseq_->qual.l == 0) {  // fasta file
          sequence_batch_[sequence_index].Update(length, sequence_kseq_->name.s,
                                                 sequence_kseq_->seq.s);
        } else {  // fastq file
          sequence_batch_[sequence_index].Update(length, sequence_kseq_->name.s,
                                                 sequence_kseq_->seq.s,
                                                 sequence_kseq_->qual.s);
        }
        ++num_loaded_sequences_;
      } else {
        assert(length == -1);  // make sure to reach the end of file rather than
                               // meet an error
        break;
      }
    }
    std::cerr << "Number of sequences: " << num_loaded_sequences_ << "."
              << std::endl;
    std::cerr << "Loaded sequence batch successfully in "
              << sga::GetRealTime() - real_start_time << "s." << std::endl;
    return num_loaded_sequences_;
  }

  void FinalizeLoading() {
    if (sequence_file_ != nullptr) {
      fclose(sequence_file_);
      sequence_file_ = nullptr;
    }

    if (sequence_kseq_ != nullptr) {
      kseq_destroy(sequence_kseq_);
      sequence_kseq_ = nullptr;
    }
  }

 protected:
  uint32_t max_batch_size_;
  uint32_t num_loaded_sequences_;
  std::vector<sga::Sequence> sequence_batch_;
  std::string sequence_file_path_;
  FILE *sequence_file_;
  kseq_t *sequence_kseq_;
};

}  // namespace sga

#endif  // SGA_SEQUENCEBATCH_H_
