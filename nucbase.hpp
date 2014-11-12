#ifndef NUCBASE_HPP
#define NUCBASE_HPP

#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include "nucsequences.hpp"
using namespace std;

// Utility functions

enum fq_encoding { SANGER=0, SOLEXA=1, IL13=2, IL15=3, IL18=4 };

void fastq2txt(string inputname, string & adapter3, string & adapter5,
               fq_encoding encoding, int minsize=0, int maxsize=0, int score=0);

void fasta2txt(string inputname, string & adapter3, string & adapter5,
               int minsize=0, int maxsize=0);


class NucBase
{
  // Class attributes (obviously protected)
  protected:
    string         _inputname;
    string         _dataname;
    string         _outputfolder;
    vector<string> _labels;
    bool           _labelled;
    int            _colmapnum;
    int            _colname;
    int            _nlines;
    int            _maxsize;
  
  
  
  // Private default constructors (to prevent misuse)
  private:
    NucBase();
    NucBase(const NucBase & data);
  
  
  
  // Public methods
  public:
  
    // Constructor
    NucBase( string inputname, 
             string outputfolder = "./Results/" );
    
    // Searches words from the database in the sequence
    bool search( NucSequences & sequences,
                 const vector<int> & columns,
                 int mismatch,
                 int submatch,
                 bool absent,
                 bool seqfile,
                 bool mapnum,
                 vector<int> &progress) const;

    // Gives the columns names
    void getLabels(vector<string> & labels) { labels.clear(); labels = _labels; }

    // Puts labels in the output files
    void labelOutputs(ofstream * output, const vector<int> & columns, const vector<string> & labels, string seqname) const;

    // Checks the desired columns
    void checkColumns(vector<int> & columns);

    // Gets the number of lines
    int getNlines() const { return _nlines; }
  
  
  
  // Protected methods
  protected:

    // Consensus-expansion
    void expand();
    
    // Saves the sequences without matching parts
    void saveChangedSequences(const vector<int> & columns, NucSequences &sequences,
                              int mismatches, int submatches, bool absent) const;
  
  
  
  // Protected Template const methods
  protected:
  
    // Opens and browses the input file
    template <bool MAPNUM, bool MISMATCHES, bool SUBMATCHES, bool BWT>
    bool processDatabase(NucSequences & sequences, const vector<int> & columns, int mismatch, int submatch, bool absent, vector<int> & progress) const;
    
    // Outputs one search result
    template <bool GFF3, bool SUBMATCHES>
    void writeOutput(ofstream & out, NucQuery & query, NucSequence & sequence, bool absent, const string & val, const string & mapnum) const;
};

#include "nucbase.hxx"

#endif
