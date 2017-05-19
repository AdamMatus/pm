#ifndef _SPEAKER_H_
#define _SPEAKER_H_

#include <vector>
#include <array>
#include <string>

template <int C, int K>
using codebook_array = typename std::array<std::array<double, K>, C>; 

template <int C, int K>
class speaker {
  public:
    speaker(const std::string& name)
      : speaker_info{name} {}  

    std::string name() {return speaker_info.name;};

    struct Code{
      Code(const std::string& text, codebook_array<C,K>&& code)
        : centroids{code}, text{text} {}
      std::array<std::array<double, K>, C> centroids;
      std::string text;
    };

    void add_code(speaker::Code&& code)
    {
      codebook.push_back(code); 
    }

    std::vector<Code> codebook;
  private:

    struct speaker_info{
      speaker_info(const std::string& name)
        : name{name} {}
      std::string name;
    } speaker_info;
};

#endif
