#ifndef NUCSEQUENCE_HPP
#define NUCSEQUENCE_HPP

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <divsufsort.h>
#include <map>
#include <list>
using namespace std;

#define MYBLOCKSIZE 18


namespace Nuc {
    string complementary(const string & seq);

    inline map<char,unsigned char> index()
    {
        map<char,unsigned char> res;
        res[ 0 ] = 0; // Needed character (end of sequence)
        res['A'] = 1;
        res['C'] = 2;
        res['G'] = 3;
        res['N'] = 4;
        res['T'] = 5;
        return res;
    }
}


// Class for queries
class NucQuery
{
  protected:
    string _name;
    string _sequence;
    bool   _sense;
    vector<saidx_t> _positions;

  public:
    // Only used when looking for submatches using a linked list
    NucQuery * next;

  public:
    // Constructor
    NucQuery() : next(0) {}

    // Destructor
    ~NucQuery();

    // Getters
    inline const string & name()     { return _name;             }
    inline const string & sequence() { return _sequence;         }
    inline const bool   & sense()    { return _sense;            }
    inline       int      count()    { return _positions.size(); }

    // Setters
    inline void name    ( const string& name )      { _name = name;         }
    inline void sequence( const string& sequence )  { _sequence = sequence; }
    inline void sense   ( const bool& sense )       { _sense = sense;       }

    // Positions handling
    inline void    addPosition(saidx_t pos)    { _positions.push_back(pos);                }
    inline void    removePosition(saidx_t pos) { _positions.erase(_positions.begin()+pos); }
    inline saidx_t position(int k)             { return _positions[k];                     }
};


class NucSequence
{
  protected:
    string _name;
    string _sequence;

  protected:
    map<char,unsigned char> _nuc;
    short _nchar;
    saidx_t _seqsize;
    vector<saidx_t> _C;
    vector<saidx_t> _SA;
    vector<saidx_t> _occ; // Used as 2D vector, but cleaner with 1D-vector
    short _blocksize;
    saidx_t _pidx;

  // No default constructor (no empty object)
  private:
    NucSequence();

  // Constructors
  public:
    NucSequence(string & seqfilename);

    NucSequence(string name, string sequence) :
      _name(name), _sequence(sequence), _nuc(Nuc::index()), _nchar(_nuc.size()), _C(_nchar),
      _blocksize(MYBLOCKSIZE)
    {
      lowercasename();
      if(!check()) throw invalid_argument( "Invalid characters in the sequence." );
    }

    // Comparison operators (only names are taken into account)
    friend bool operator== (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name == seq2._name; }
    friend bool operator!= (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name != seq2._name; }
    friend bool operator<  (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name < seq2._name; }
    friend bool operator>= (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name >= seq2._name; }
    friend bool operator>  (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name > seq2._name; }
    friend bool operator<= (const NucSequence &seq1, const NucSequence &seq2) { return seq1._name <= seq2._name; }

    // Const Getters
    const string & name()          { return _name;         }
    const string & sequence()      { return _sequence;     }

    // Setters
    void name(const string & name) { _name = name; }

    // BWT
    void bwt();
    void inverse_bwt();

    template <bool SUBMATCHES,bool MISMATCHES, bool BWT>
    void search(NucQuery & query, const int & mismatches, const int & submatches);

    template <bool MISMATCHES, bool BWT>
    void search(NucQuery & query, const int & mismatches);

  protected:
    // Puts the name in lower case
    void lowercasename() { transform(_name.begin(), _name.end(), _name.begin(), (int (*)(int))tolower); }
    // Checks the sequence
    bool check();
};


class NucSequences : public vector<NucSequence>
{
  // FASTA loading constructor
  public:
    NucSequences() {}
    NucSequences(string & fastafilename);
};


struct Candidates
{
    int count;
    saidx_t low;
    saidx_t high;
    Candidates(int c, saidx_t l, saidx_t h, Candidates * n) : count(c), low(l), high(h), next(n) {}

    Candidates * next;
};

#include "nucsequences.hxx"

#endif
