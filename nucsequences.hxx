#ifndef NUCSEQUENCES_HXX
#define NUCSEQUENCES_HXX

#include "nucsequences.hpp"


template <bool SUBMATCHES, bool MISMATCHES, bool BWT>
void NucSequence::search(NucQuery & query, const int & mismatches, const int & submatches)
{
  // Alias to the sequence we are looking for
  const string & word = query.sequence();
  // Size of the word
  saidx_t size = word.size();

  search<MISMATCHES,BWT>(query,mismatches);

  if(SUBMATCHES)
  {
    if(size >= submatches)
    {
      NucQuery * last = &query;

      for(int i=size-submatches; i>=0; --i)
      {
        last->next = new NucQuery;

        NucQuery * elt = last->next;
        elt->name(query.name());
        elt->sequence(word.substr(i, submatches));
        elt->sense(query.sense());

        search<MISMATCHES,BWT>(*elt,mismatches);

        last = elt;
      }

      last->next = 0;

      // We merge adjacent submatches
      NucQuery * elt = query.next;
      NucQuery * nxt = elt->next;
      while(nxt != 0)
      {
        string str1(elt->sequence(),0,elt->sequence().size()-1);
        string str2(nxt->sequence(),1,nxt->sequence().size()-1);

        NucQuery * lastprev = last;
        if(str1 == str2)
        {
          vector<int> ind;
          for(int i=0;i<elt->count(); ++i)
            for(int j=0;j<nxt->count(); ++j)
              if(elt->position(i) == (nxt->position(j)+1))
                ind.push_back(i);

          if(ind.size() > 0)
          {
            string tmp(nxt->sequence(),0,1);
            tmp += elt->sequence();

            last->next = new NucQuery;
            last = last->next;

            last->name(query.name());
            last->sequence(tmp);
            last->sense(query.sense());
            last->next = 0;

            for(int i=ind.size()-1;i>=0; --i)
            {
              last->addPosition(elt->position(ind[i])-1);
              elt->removePosition(ind[i]);
            }
          }
        }

        if(str1 == str2 || elt->sequence().size() < nxt->sequence().size())
        {
          vector<int> ind2rm;
          for(int i=0;i<elt->count(); ++i)
            for(int j=0;j<lastprev->count(); ++j)
              if(elt->position(i) == lastprev->position(j))
                ind2rm.push_back(i);

          for(int i=ind2rm.size()-1;i>=0; --i)
            elt->removePosition(ind2rm[i]);
        }
        elt = nxt;
        nxt = elt->next;
      }

      // List cleaning
      elt = &query;
      nxt = elt->next;
      while(nxt != last)
      {
        if(nxt->count() == 0)
        {
          elt->next = nxt->next;
          nxt->next = 0;
          delete nxt;
        }
        else
          elt = nxt;
        nxt = elt->next;
      }
      if(last->sequence() == query.sequence() || last->count() == 0)
      {
        elt->next = 0;
        delete last;
      }
    }
  }
}


template <bool MISMATCHES, bool BWT>
void NucSequence::search(NucQuery & query, const int & mismatches)
{
  // Alias to the sequence we are looking for
  const string & word = query.sequence();
  // Size of the word
  saidx_t size = word.size();

  if(BWT)
  {
    if(!MISMATCHES)
    {
      // We initialize the loop variables
      // Low index
      saidx_t low = 0;
      // High index
      saidx_t high = _seqsize+1;

      // We search for character in ith position
      // with consideration to the previous character treated
      for(saidx_t i=size-1; i>=0 && low < high; --i)
      {
        // ith character in word
        char c = word[i];
        // Corresponding index
        short ic = _nuc[c];

        // Occurrences (table divided in blocks)
        // Blocks indexes
        saidx_t lowb = low/_blocksize;
        saidx_t highb = high/_blocksize;

        // Blocks values
        saidx_t lowocc = _occ[ic+lowb*_nchar];
        saidx_t highocc = _occ[ic+highb*_nchar];

        // Remaining characters to browse in BWT
        short lowmodb = low%_blocksize;
        short highmodb = high%_blocksize;

        // Counting remaining characters in BWT (low)
        for(short a=0; a<lowmodb; ++a)
          if(_sequence[a+lowb*_blocksize] == c)
            ++lowocc;

        // Counting remaining characters in BWT (high)
        for(short a=0; a<highmodb; ++a)
          if(_sequence[a+highb*_blocksize] == c)
            ++highocc;

        // New low and high indexes
        low = _C[ic] + lowocc;
        high = _C[ic] + highocc;
      }

      // We store their positions
      for(saidx_t k=low; k<high; ++k)
        query.addPosition(_SA[k]);
    }
    else
    {
      // We initialize the loop variables
      // Character being treated
      saidx_t i = size-1;

      //list<Candidate> candidates;
      //candidates.push_front(Candidate(0,0,_seqsize+1));
      Candidates * init = 0;
      Candidates * lst = new Candidates(0,(saidx_t)0,_seqsize+1,init);

      // We search for character in ith position
      // with consideration to the previous character treated
      while(lst != 0 && i >= 0)
      {
        Candidates * ptr = lst;
        Candidates ** prev = &lst;

        while(ptr != 0)
        {
          // We get the candidate triplet
          int count = ptr->count;
          saidx_t clow = ptr->low;
          saidx_t chigh = ptr->high;

          // Map iteration
          map<char,unsigned char>::iterator it = _nuc.begin();
          ++it;
          while(it!=_nuc.end())
          {
            char c = it->first;
            char ic = it->second;

            // Local low & high
            saidx_t high = chigh;
            saidx_t low = clow;

            // Occurrences (table divided in blocks)
            // Blocks indexes
            saidx_t lowb = low/_blocksize;
            saidx_t highb = high/_blocksize;

            // Blocks values
            saidx_t lowocc = _occ[ic+lowb*_nchar];
            saidx_t highocc = _occ[ic+highb*_nchar];

            // Remaining characters to browse in BWT
            short lowmodb = low%_blocksize;
            short highmodb = high%_blocksize;

            // Counting remaining characters in BWT (low)
            for(short a=0; a<lowmodb; ++a)
              if(_sequence[a+lowb*_blocksize] == c)
                ++lowocc;

            // Counting remaining characters in BWT (high)
            for(short a=0; a<highmodb; ++a)
              if(_sequence[a+highb*_blocksize] == c)
                ++highocc;

            // New low and high indexes
            low = _C[ic] + lowocc;
            high = _C[ic] + highocc;

            if(c != word[i] && low < high && count<mismatches)
            {
              lst = new Candidates(count+1,low,high,lst);
              if(lst->next == ptr)
                prev = &(lst->next);
            }
            else if(c == word[i] && low < high)
            {
              lst = new Candidates(count,low,high,lst);
              if(lst->next == ptr)
                prev = &(lst->next);
            }

            ++it;
          }

          *prev = ptr->next;

          delete ptr;
          ptr = *prev;
        }

        // Previous character
        --i;
      }

      // We store their positions
      Candidates * ptr = lst;
      while(ptr != 0)
      {
        for(saidx_t k=ptr->low; k<ptr->high; ++k)
          query.addPosition(_SA[k]);

        Candidates * prev = ptr;
        ptr = ptr->next;
        delete prev;
      }
    }
  }
  else
  {
    if(!MISMATCHES)
    {
      size_t pos = string::npos;

      do {
        pos = _sequence.find(word,pos+1);

        if(pos != string::npos)
          query.addPosition((saidx_t)pos);

      } while(pos != string::npos);
    }
    else
    {
      saidx_t seqsize = _sequence.size();
      if(seqsize >= size)
      {
        saidx_t searchedseqsize = seqsize - size;

        // For each position, we check the number of mismatches
        for(saidx_t pos=0; pos <= searchedseqsize; ++pos)
        {
          int miss = 0;
          for(saidx_t j=0; j<size; ++j)
            if(word[j] != _sequence[j+pos])
              ++miss;

          // If miss <= mismatch, we found one more.
          if(miss <= mismatches)
            query.addPosition(pos);
        }
      }
    }
  }
}


#endif // NUCSEQUENCES_HXX
