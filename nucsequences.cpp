#include "nucsequences.hpp"
#include <fstream>
using namespace std;


namespace Nuc
{
  string complementary(const string &seq)
  {
    string res;
    for(string::const_reverse_iterator rit=seq.rbegin(); rit<seq.rend(); ++rit)
    {
      char temp = 'N';
      switch(*rit)
      {
        case 'T': temp = 'A'; break;
        case 'A': temp = 'T'; break;
        case 'G': temp = 'C'; break;
        case 'C': temp = 'G'; break;
        // Unused
        //case 'K': temp = 'M'; break;
        //case 'S': temp = 'W'; break;
        //case 'Y': temp = 'R'; break;
        //case 'M': temp = 'K'; break;
        //case 'W': temp = 'S'; break;
        //case 'R': temp = 'Y'; break;
        //case 'B': temp = 'V'; break;
        //case 'D': temp = 'H'; break;
        //case 'H': temp = 'D'; break;
        //case 'V': temp = 'B'; break;
        //case 'U': temp = 'A'; break;
      }
      res.push_back(temp);
    }
    return res;
  }
}


NucQuery::~NucQuery()
{
  if(next != 0)
    delete next;
  next = 0;
}


NucSequence::NucSequence(string & seqfilename) : _nuc(Nuc::index()), _nchar(_nuc.size()), _C(_nchar), _blocksize(MYBLOCKSIZE)
{
  // We get the sequence name
  size_t length = string::npos;
  size_t startname = seqfilename.rfind('/');
  size_t endname = seqfilename.rfind('.');
  if(endname != string::npos)
    length = endname - startname - 1;
  
  if(length == 0)
    throw ios::failure( "Error opening sequence file : empty sequence name..." );
  
  _name = seqfilename.substr(startname+1, length);

  // We convert the name to lower case
  lowercasename();

  // We open the file
  ifstream sequencefile(seqfilename.c_str());
  if(sequencefile.is_open())
  {
    string line;

    // We get the sequence
    while(getline(sequencefile, line))
      _sequence += line;

    sequencefile.close();

    if(_sequence.size() == 0)
      throw invalid_argument( "Empty sequence." );

    // We check if the sequence is correct
    if(!check())
      throw invalid_argument( "Invalid characters in the sequence." );
  }
  else
    throw ios::failure( "Error opening sequence file !" );
}


bool NucSequence::check()
{
  bool ok = true;

  // We convert the sequence to upper case
  transform(_sequence.begin(), _sequence.end(), _sequence.begin(), (int (*)(int))toupper);

  // We check if the characters are known
  for(string::reverse_iterator rit=_sequence.rbegin(); rit<_sequence.rend(); ++rit)
    if(_nuc.count(*rit)<1)
      ok = false;

  return ok;
}


void NucSequence::bwt()
{
  // Append end character
  _sequence.append(1,(char)0);

  // Sizes
  _seqsize = _sequence.size();
  saidx_t nb = (_seqsize+1)/_blocksize;

  // Resize Suffix Array
  _SA.resize(_seqsize,0);

  // Temporarily use arrays for compatibility with libdivsufsort
  // WARNING : we use the fact that data in vectors and strings is contiguous
  saidx_t * SA = &_SA[0];
  sauchar_t * str = (sauchar_t *)&_sequence[0];

  // Suffix array computation
  divsufsort(str, SA, _seqsize);

  // Index construction
  map<char,unsigned char>::iterator it;
  for(it = _nuc.begin(); it != _nuc.end(); ++it)
    sa_simplesearch(str, _seqsize, SA, _seqsize, it->first, &_C[it->second]);

  // BWT
  _pidx = divbwt(str, str, NULL, _seqsize);

  // Bug correction : end-character always at the start of BWT
  saidx_t end = 0;
  while(SA[end] != 0 && end < _seqsize)
    ++end;

  for(saidx_t o=0; o<end; ++o)
    _sequence[o] = _sequence[o+1];
  _sequence[end] = (char) 0;
  // End correction

  // Occurrences table (divided in blocks to save space)
  _occ = vector<saidx_t>(_nchar*(nb+1),0);
  // For each block
  for(saidx_t ind=1; ind<=nb; ++ind)
  {
    // We initialize the block value to the previous block
    for(short a=0; a<_nchar; ++a)
      _occ[a + ind*_nchar] = _occ[a + (ind-1)*_nchar];

    // We then add the characters occurrences between them
    for(short a=0; a<_blocksize; ++a)
    {
      // Character in the BWT string
      char c = _sequence[a + (ind-1)*_blocksize];
      // Increment counter for this block
      ++_occ[_nuc[c] + ind*_nchar];
    }
  }
}


void NucSequence::inverse_bwt()
{
  // Bug correction : re-place end-character at the start of BWT
  saidx_t end = 0;
  while(_SA[end] != 0 && end < _seqsize)
    ++end;

  for(saidx_t o=end; o>0; --o)
    _sequence[o] = _sequence[o-1];
  _sequence[0] = (char) 0;
  // End correction

  if(_seqsize > 0)
  {
    // WARNING : We use data contiguity in strings
    sauchar_t * str = (sauchar_t *)&_sequence[0];
    inverse_bw_transform(str, str, NULL, _seqsize, _pidx);
  }

  _seqsize = _sequence.size()-1;
  _sequence = _sequence.substr(0,_seqsize);

  _occ.clear();
  _SA.clear();
}


NucSequences::NucSequences(string & filename)
{
  // We open the file
  ifstream file(filename.c_str());
  if(file.is_open())
  {
    string line;
    string name = "";
    string sense = "";

    bool noend = getline(file, line);

    //We check the first line (FASTA or TXT)
    if(line[0] == '>')
    {
      while(noend)
      {
        if(line[0] == '>')
        {
          size_t startname = 1;
          size_t endname = line.find_first_of(" |/\\;.,",startname);

          size_t posname = line.find("name=");
          if(posname != string::npos)
          {
            startname = posname+5;
            endname = line.find_first_of("; |/\\.,",startname);
          }

          name = line.substr(startname,endname - startname);

          noend = getline(file, line);
        }
        else
          throw invalid_argument( "Unexpected format (FASTA file)." );

        while(noend && line[0] != '>')
        {
          sense += line;
          noend = getline(file, line);
        }

        if(name != "" && sense != "")
          this->push_back(NucSequence(name,sense));
        else
          throw invalid_argument( "Empty name or sequence (FASTA file)." );

        name.clear();
        sense.clear();
      }
      file.close();
    }
    else
    {
      file.close();

      this->push_back(NucSequence(filename));
    }
  }
}
