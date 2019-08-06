#ifndef SGA_SEQUENCE_H_
#define SGA_SEQUENCE_H_

#include <cstring>
#include <iostream>
#include <string>
#include <unistd.h>

namespace sga {
class Sequence {
 public:
  Sequence(){}
  Sequence(const uint32_t length, const char *name, const char *sequence) {
    Sequence();
    name_ = name;
    sequence_ = sequence;
  }
  Sequence(const uint32_t length, const char *name, const char *sequence, const char *quality_scores) {
    Sequence(length, name, sequence);
    quality_scores_ = quality_scores;
  }
  ~Sequence() {
    name_.clear();
    sequence_.clear();  
    quality_scores_.clear();
  }

  std::string const & GetName() const {
    return name_;
  }
  std::string const & GetSequence() const {
    return sequence_;
  }
  uint32_t GetLength() const {
    return sequence_.length();
  }
  std::string const & GetQualityScores() const {
    return quality_scores_;
  }

  inline void Update(const uint32_t length, const char *name, const char *sequence) { // for fasta
    name_.clear();
    name_ = name;
    sequence_.clear();
    sequence_ = sequence;
  }
  inline void Update(const uint32_t length, const char *name, const char *sequence, const char *quality_scores) { // for fastq
    Update(length, name, sequence);
    quality_scores_.clear();
    quality_scores_ = quality_scores;
  }

 protected:
  std::string name_; 
  std::string sequence_;
  std::string quality_scores_;
};

} // namespace sga

#endif // SGA_SEQUENCE_H_
