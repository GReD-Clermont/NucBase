#include "nucbase.hpp"
#include <errno.h>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
using namespace std;

// Macro to manage mkdir call (Windows vs Unix)
#ifdef _WIN32
#include <direct.h>
#define MKDIR(PATH) mkdir(PATH)
#else
#include <sys/stat.h>
#define MKDIR(PATH) mkdir(PATH, 0775)
#endif

// Utility functions
void fastq2txt(string inputname, string & adapter3, string & adapter5, fq_encoding encoding, int minsize, int maxsize, int score)
{
  map<string, int> data;
  char base_score = 33;
  bool nomaxsize = (maxsize == 0);
  bool ad3 = !adapter3.empty();
  bool ad5 = !adapter5.empty();

  switch(encoding)
  {
    case SANGER : base_score = 33; break;
    case SOLEXA : base_score = 59; break;
    case IL13 : base_score = 64; break;
    case IL15 : base_score = 64; break;
    case IL18 : base_score = 33; break;
  }

  string outputname = inputname;
  size_t extpos = outputname.rfind('.');
  size_t slashpos = outputname.rfind('/');
  outputname.replace(extpos+1,outputname.size()-extpos,"txt");

  ifstream fastq(inputname.c_str());
  if(fastq.is_open())
  {
    string line;
    while(getline(fastq,line))
    {
      string seqid;
      string seq;
      string scores;
      char minscore;

      if(line[0] == '@')
      {
        seqid = line.substr(1);
        if(getline(fastq,line))
        {
          seq = line;
          if(getline(fastq,line) && line[0] == '+')
          {
            if(getline(fastq,line))
            {
              scores = line;

              if(ad5)
              {
                size_t ad5pos = seq.find(adapter5);
                scores = scores.substr(ad5pos);
                seq = seq.substr(ad5pos);
              }

              if(ad3)
              {
                size_t ad3pos = 0;
                size_t size = adapter3.size();
                do
                {
                  ad3pos = seq.find(adapter3.substr(0,size));
                  --size;
                } while(ad3pos == string::npos && size >= 8);

                if(size == 8 || ad3pos == string::npos)
                {
                  int mm = 0;
                  size = adapter3.size();
                  ad3pos = seq.find(adapter3.substr(0,5));
                  if(ad3pos != string::npos)
                  {
                    for(size_t i=5; i<size; ++i)
                      if(seq[ad3pos+i] != adapter3[i])
                        ++mm;

                    if( (float)mm/(float)(size-5) < 0.2 )
                    {
                      scores = scores.substr(0,ad3pos);
                      seq = seq.substr(0,ad3pos);
                    }
                  }
                }
                else
                {
                  scores = scores.substr(0,ad3pos);
                  seq = seq.substr(0,ad3pos);
                }
              }

              if(encoding == IL15)
              {
                size_t posB = scores.find('B');
                scores = scores.substr(0, posB);
                seq = seq.substr(0, posB);
              }
            }
            else
              throw ios::failure( "Error while reading the fastq file ! (4)" );
          }
          else
            throw ios::failure( "Error while reading the fastq file ! (3)" );
        }
        else
          throw ios::failure( "Error while reading the fastq file ! (2)" );

        minscore = (char)127;
        for(string::iterator it=scores.begin(); it != scores.end(); ++it)
          if(minscore > *it)
            minscore = *it;

        if((int)(minscore-base_score) >= score)
          if(nomaxsize || (size_t)maxsize >= seq.size())
            if((size_t)minsize <= seq.size())
              data[seq]++;
      }
    }
    fastq.close();

    ofstream txt(outputname.c_str());
    if(txt.is_open())
    {
      txt << "labels\t" << inputname.substr(slashpos+1, extpos-slashpos-1) << endl;

      for(map<string,int>::iterator it=data.begin(); it!=data.end(); ++it)
        txt << it->first << "\t" << it->second << endl;

      txt.close();
    }
  }
  else
    throw ios::failure( "Error while reading the fastq file !" );

}

void fasta2txt(string inputname, string & adapter3, string & adapter5, int minsize, int maxsize)
{
  map<string, int> data;
  bool nomaxsize = (maxsize == 0);
  bool ad3 = !adapter3.empty();
  bool ad5 = !adapter5.empty();

  string outputname = inputname;
  size_t extpos = outputname.rfind('.');
  size_t slashpos = outputname.rfind('/');
  outputname.replace(extpos+1,outputname.size()-extpos,"txt");

  ifstream fasta(inputname.c_str());
  if(fasta.is_open())
  {
    string line;
    string seqname;
    string seq = "";
    while(getline(fasta,line))
    {
      if(line[0] == '>')
      {
        if(ad5)
          seq = seq.substr(seq.find(adapter5));

        if(ad3)
        {
          size_t ad3pos;
          size_t size = adapter3.size();
          do
          {
            ad3pos = seq.find(adapter3.substr(0,size));
            --size;
          } while(ad3pos == string::npos && size >= 8);

          if(size == 8 || ad3pos == string::npos)
          {
            int mm = 0;
            size = adapter3.size();
            ad3pos = seq.find(adapter3.substr(0,5));
            if(ad3pos != string::npos)
            {
              for(size_t i=5; i<size; ++i)
                if(seq[ad3pos+i] != adapter3[i])
                  ++mm;

              if( (float)mm/(float)(size-5) < 0.2 )
                seq = seq.substr(0,ad3pos);
            }
            else
              seq = seq.substr(0,ad3pos);
          }
        }

        if(nomaxsize || (size_t)maxsize >= seq.size())
          if((size_t)minsize <= seq.size())
            data[seq]++;

        seq = "";
        seqname = line.substr(1);
      }
      else
        seq += line;
    }
    fasta.close();

    ofstream txt(outputname.c_str());
    if(txt.is_open())
    {
      txt << "labels\t" << inputname.substr(slashpos+1, extpos-slashpos-1) << endl;

      for(map<string,int>::iterator it=data.begin(); it!=data.end(); ++it)
        txt << it->first << "\t" << it->second << endl;

      txt.close();
    }
  }
}

NucBase::NucBase( string inputname, string outputfolder ) : 
  _inputname(inputname), _dataname("data"), _outputfolder(outputfolder),
  _labelled(false), _colmapnum(0), _colname(0),
  _nlines(1), _maxsize(0)
{
  bool invalid = false;
  bool consensus = false;
  string cons_char = "UKSYMWRBDHVN";
  string accepted_chars = "ACGTUKSYMWRBDHVN";

  // We give the database a name
  size_t length = string::npos;
  size_t startname = _inputname.rfind('/');
  size_t endname = _inputname.rfind('.');
  if(endname != string::npos)
    length = endname - startname - 1;

  if(length != 0)
    _dataname = _inputname.substr(startname+1, length);

  int mkd = MKDIR(outputfolder.c_str());

  if(mkd == 0 || errno == EEXIST)
  {
    // We now open the database to get some information
    ifstream input(_inputname.c_str());
    if(input.is_open())
    {
      vector<string> words;
      string line;
      string word;
      size_t found;
      size_t bad_char;

      // We parse the first line
      getline(input, line);
      if (line[line.size() - 1] == '\r')
        line.resize(line.size() - 1);

      stringstream strstr(line);
      while (getline(strstr, word, '\t'))
        words.push_back(word);

      // Number of columns in input files
      int nbcol = words.size();

      // We check that we have at least 1 column
      if(nbcol > 0)
      {
        // We save the first word for later
        word = words[0];

        // We convert the first line words to lower case
        for(int i=0; i<nbcol; ++i)
          std::transform(words[i].begin(), words[i].end(), words[i].begin(), (int (*)(int))tolower);

        // If the first line contains the labels
        if(words[0] == "labels")
        {
          _labelled = true;

          if(nbcol > 1)
          {
            _labels = words;

            // We don't count the first line
            _nlines--;

            for(int i=0; i<nbcol; ++i)
              if(words[i] == "mapnum")
                _colmapnum = i;

            for(int i=0; i<nbcol; ++i)
              if(words[i] == "name")
                _colname = i;
          }
          else
            _labels.push_back(_dataname);
        }
        else
        {
          // We create generic names for the columns
          if(nbcol>2)
          {
            for(int i=0; i<nbcol; ++i)
            {
              std::ostringstream oss;
              oss << "Column " << i;
              _labels.push_back(oss.str());
            }
          }
          else
          {
            if(nbcol > 1)
              _labels.push_back("labels");

            _labels.push_back(_dataname);
          }

          // We check if the string is valid
          bad_char = word.find_first_not_of(accepted_chars);
          invalid |= bad_char != string::npos;

          // We check for special notation characters in the first line
          found = word.find_first_of(cons_char);
          consensus |= found != string::npos;

          // We look for the largest read
          if(_maxsize < (int)word.size())
            _maxsize = (int)word.size();
        }

        // We count the number of lines (and keep checking)
        while(getline(input, line) && !invalid)
        {
          ++_nlines;

          stringstream strstr(line);
          getline(strstr, word, '\t');

          // We check if the word is valid
          bad_char = word.find_first_not_of(accepted_chars);
          invalid |= bad_char != string::npos;

          // We look for degenerate bases
          size_t tab = line.find('\t');
          found = line.find_first_of(cons_char);

          // We look for the largest read
          if(_maxsize < (int)word.size())
            _maxsize = (int)word.size();

          if(!consensus)
            consensus = found < tab;
        }

        // We close the file
        input.close();

        if(invalid)
          throw invalid_argument("Invalid characters in the database.");
        else
          if(consensus)
            expand();
      }
      else
        invalid = true;

    }
    else
      throw ios::failure( "Error opening database file !" );
  }
  else
    throw ios::failure("Could not create \"Results\" folder.");
}


void NucBase::expand()
{
  string cons_char = "UKSYMWRBDHVN";

  // Map
  typedef multimap<char, char> MyMap;
  MyMap mapping;
  mapping.insert(pair<char,char>('A','A'));
  mapping.insert(pair<char,char>('C','C'));
  mapping.insert(pair<char,char>('G','G'));
  mapping.insert(pair<char,char>('T','T'));
  mapping.insert(pair<char,char>('U','T'));
  mapping.insert(pair<char,char>('K','G'));
  mapping.insert(pair<char,char>('K','T'));
  mapping.insert(pair<char,char>('S','G'));
  mapping.insert(pair<char,char>('S','C'));
  mapping.insert(pair<char,char>('Y','C'));
  mapping.insert(pair<char,char>('Y','T'));
  mapping.insert(pair<char,char>('M','A'));
  mapping.insert(pair<char,char>('M','C'));
  mapping.insert(pair<char,char>('W','A'));
  mapping.insert(pair<char,char>('W','T'));
  mapping.insert(pair<char,char>('R','A'));
  mapping.insert(pair<char,char>('R','G'));
  mapping.insert(pair<char,char>('B','C'));
  mapping.insert(pair<char,char>('B','G'));
  mapping.insert(pair<char,char>('B','T'));
  mapping.insert(pair<char,char>('D','A'));
  mapping.insert(pair<char,char>('D','G'));
  mapping.insert(pair<char,char>('D','T'));
  mapping.insert(pair<char,char>('H','A'));
  mapping.insert(pair<char,char>('H','C'));
  mapping.insert(pair<char,char>('H','T'));
  mapping.insert(pair<char,char>('V','A'));
  mapping.insert(pair<char,char>('V','C'));
  mapping.insert(pair<char,char>('V','G'));
  mapping.insert(pair<char,char>('N','A'));
  mapping.insert(pair<char,char>('N','C'));
  mapping.insert(pair<char,char>('N','G'));
  mapping.insert(pair<char,char>('N','T'));

  // Map iterator
  pair<MyMap::iterator,MyMap::iterator> mapit;

  ifstream input(_inputname.c_str());
  ofstream output("tmp.txt");

  if(input.is_open() && output.is_open())
  {
    _nlines = 0;
    string line;

    if(_labelled)
    {
      getline(input,line);
      output << line << endl;
    }

    // We go through the input file
    while(getline(input, line))
    {
      vector<string> words;
      string word;
      string * expansion;
      int times = 1;

      stringstream strstr(line);
      while (getline(strstr, word, '\t'))
        words.push_back(word);
      size_t wordsize = words[0].size();

      // We count how many words there will be instead
      size_t found = words[0].find_first_of(cons_char);
      while(found != string::npos)
      {
        times *= (int)mapping.count(words[0][found]);
        found = words[0].find_first_of(cons_char, found+1);
      }
      int n = times;

      // Array allocation;
      expansion = new string[times];

      // Array initialization
      for(int i=0; i<times; ++i)
        expansion[i] = words[0];

      // We do the replacement
      for(size_t j=0; j<wordsize; ++j)
      {
        vector<char> repl;
        int count = mapping.count(words[0][j]);
        mapit = mapping.equal_range(words[0][j]);
        n /= count;

        multimap<char,char>::iterator it;
        for(it=mapit.first; it!=mapit.second; ++it)
          repl.push_back(it->second);

        for(int i=0; i<times; ++i)
          expansion[i][j] = repl[(i/n)%count];
      }

      // We write the output file
      for(int i=0; i<times; ++i)
      {
        output << expansion[i];
        for(size_t j=1; j<words.size(); ++j)
          output << "\t" << words[j];
        output << endl;
        ++_nlines;
      }

      // We delete the string array
      delete [] expansion;
    }
    output.close();
    input.close();

    string old(_inputname);
    old += ".old";

    // We remove the old "old file" if it exists
    remove(old.c_str());

    // We rename the temp file
    int res1 = rename(_inputname.c_str(), old.c_str());
    int res2 = rename("tmp.txt", _inputname.c_str());
    if(res1 != 0 || res2 != 0)
      throw ios::failure( "Could not create new database file (extended notation)...");
  }
  else
    throw ios::failure( "Error expanding the database file !" );
}


bool NucBase::search( NucSequences & sequences, const vector<int> & columns, int mismatch, int submatch, bool absent, bool seqfile, bool mapnum, vector<int> & progress ) const
{
  bool ok = false;
  bool bwt = true;
  int nseq = sequences.size();
  vector<int> ind2remove;

  // Statistics to evaluate cost
  int maxseqsize = 0;
  double meanseqsize = 0;

  // For each sequence
  for(int i=0; i<nseq; ++i)
  {
    // We create the sequence output folder, in case it doesn't exist
    string seqresdir(_outputfolder);
    seqresdir += sequences[i].name();
    int mkd = MKDIR(seqresdir.c_str());

    // If name forbidden by OS we try with an underscore
    if(mkd != 0 && errno != EEXIST)
    {
      string newseqname = sequences[i].name();
      newseqname.insert(0,"_");
      seqresdir = _outputfolder;
      seqresdir += newseqname;
      mkd = MKDIR(seqresdir.c_str());

      // If we still can't create it, we throw an error and remember the sequence to remove
      if(mkd != 0 && errno != EEXIST)
      {
        ind2remove.push_back(i);
        throw ios::failure("Could not create sequence output folder for : "+sequences[i].name()+"\n");
      }
      else
        // We keep this name
        sequences[i].name(newseqname);
    }

    // We get the size of the largest sequence
    if((int)sequences[i].sequence().size() > maxseqsize)
      maxseqsize = sequences[i].sequence().size();

    // We compute the mean size
    meanseqsize += sequences[i].sequence().size();
  }

  // Mean size of the sequences
  meanseqsize /= nseq;

  // We remove previously marked sequences
  for(vector<int>::reverse_iterator it=ind2remove.rbegin(); it!=ind2remove.rend(); --it)
    sequences.erase(sequences.begin()+*it);

  // We estimate the costs of the bwt and "naive" methods
  long long bwt_cost = _maxsize;
  long long std_cost = meanseqsize;
  if(mismatch > 0)
  {
    std_cost = maxseqsize*_maxsize;
    bwt_cost = _maxsize*(_maxsize-mismatch)*pow(4,mismatch)*(log(maxseqsize)-log(4));
    if(bwt_cost > _maxsize*maxseqsize*(log(maxseqsize)-log(4)))
      bwt_cost = _maxsize*maxseqsize*(log(maxseqsize)-log(4));
  }

  // We see if bwt is more interesting
  if(_nlines*bwt_cost + meanseqsize*log(meanseqsize) > _nlines*std_cost)
    bwt = false;

  // We create a variable to sum up the options
  char options = 0;
  if(mapnum)       options += 1;
  if(mismatch > 0) options += 2;
  if(submatch > 9) options += 4; //Strings of 9 nucleotids will give too many results
  if(bwt == true ) options += 8;

  // The way we browse the database depends on the options
  switch(options)
  {
    case 0 : ok = processDatabase<false, false, false, false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 1 : ok = processDatabase<true , false, false, false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 2 : ok = processDatabase<false, true , false, false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 3 : ok = processDatabase<true , true , false, false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 4 : ok = processDatabase<false, false, true , false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 5 : ok = processDatabase<true , false, true , false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 6 : ok = processDatabase<false, true , true , false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 7 : ok = processDatabase<true , true , true , false>(sequences, columns, mismatch, submatch, absent, progress); break;
    case 8 : ok = processDatabase<false, false, false, true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 9 : ok = processDatabase<true , false, false, true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 10: ok = processDatabase<false, true , false, true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 11: ok = processDatabase<true , true , false, true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 12: ok = processDatabase<false, false, true , true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 13: ok = processDatabase<true , false, true , true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 14: ok = processDatabase<false, true , true , true >(sequences, columns, mismatch, submatch, absent, progress); break;
    case 15: ok = processDatabase<true , true , true , true >(sequences, columns, mismatch, submatch, absent, progress); break;
    default: ok = processDatabase<false, false, false, true >(sequences, columns, mismatch, submatch, absent, progress); break;
  }

  if(seqfile)
      saveChangedSequences(columns, sequences, mismatch, submatch, absent);

  return ok;
}


void NucBase::saveChangedSequences(const vector<int> & columns, NucSequences & sequences, int mismatches, int submatches, bool absent) const
{
  int  nseq = sequences.size();
  int  ncol = columns.size();

  for(int i=0; i<ncol; ++i)
  {
    for(int j=0; j<nseq; ++j)
    {
      ostringstream oss;
      oss << _outputfolder << sequences[j].name() << "/" << sequences[j].name() << "_" << _labels[columns[i]];

      if(absent)
        oss << "_absent";

      if(mismatches > 0)
        oss << "_" << mismatches << "mm";

      if(submatches > 9)
        oss << "_" << submatches << "minblock";

      string base = oss.str();
      string gff3 = base + ".gff3";
      string out = base + ".txt";

      ifstream input(gff3.c_str());
      ofstream output(out.c_str());

      if(input.is_open() && output.is_open())
      {
        string line = "";
        string tmp_s(sequences[j].sequence());
        string tmp_a(Nuc::complementary(sequences[j].sequence()));

        int seqsize = tmp_s.size();

        // We skip the first three lines of the GFF3
        getline(input,line); // "##gff_version 3"
        getline(input,line); // "##Index_subfeatures 1"
        getline(input,line); // ""

        // We go through the GFF3 file
        while(getline(input,line))
        {
          // We get and parse the line.
          string word;
          vector<string> words;
          stringstream strstr(line);
          while (getline(strstr, word, '\t'))
            words.push_back(word);

          bool sense = words[6] == "+";
          int  pos1  = atoi(words[3].c_str())-1;
          int  pos2  = atoi(words[4].c_str())-1;
          int  size  = pos2 - pos1 + 1;

          size_t posname = words[8].find("Name=");
          word = words[8].substr(posname+5,size);
          if(!sense)
            word = Nuc::complementary(word);

          if(sense)
          {
            for(int i=0; i<size; ++i)
              if(tmp_s[pos1+i] == word[i])
                tmp_s[pos1+i] = '*';
          }
          else
          {
            for(int i=0; i<size; ++i)
              if(tmp_a[seqsize-1-pos2+i] == word[i])
                tmp_a[seqsize-1-pos2+i] = '*';
          }
        }

        // Output sense
        output << ">" << sequences[j].name() << " (sense)" << endl;
        for(int k=0; k<=(int)(seqsize/80); ++k)
          output << tmp_s.substr(k*80,80) << endl;

        // Output antisense
        output << ">" << sequences[j].name() << " (antisense)" << endl;
        for(int k=0; k<=(int)(seqsize/80); ++k)
          output << tmp_a.substr(k*80,80) << endl;

        output.close();
        input.close();
      }
      else
        throw ios::failure( "ProcessDatabase : error creating \"unmatched sequences\" files !" );
    }
  }
}


void NucBase::checkColumns(vector<int> & columns)
{
  int  nlabels = _labels.size();

  // We check that each column number exists or we delete it
  for(vector<int>::iterator it = columns.begin(); it<columns.end(); ++it)
    if(*it >= nlabels || *it <= 0)
      columns.erase(it--);
}


void NucBase::labelOutputs(ofstream * output, const vector<int> & columns, const vector<string> & labels, string seqname) const
{
  int ncol = columns.size();
  string mapnum = "";

  if(_colmapnum > 0)
    mapnum = "\tmapnum";

  for(int i=0; i<ncol; ++i)
  {
      // GFF3
      output[i+0*ncol] << "##gff_version 3" << endl;
      output[i+0*ncol] << "##Index_subfeatures 1" << endl;
      output[i+0*ncol] << endl;

      // Sense
      output[i+1*ncol] << "labels\t" << labels[columns[i]] << "_on_" << seqname << "\t" << _labels[columns[i]] << mapnum << endl;

      // Antisense
      output[i+2*ncol] << "labels\t" << labels[columns[i]] << "_on_" << seqname << "\t" << _labels[columns[i]] << mapnum << endl;
  }
}
