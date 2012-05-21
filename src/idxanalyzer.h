#ifndef __idxanalyzer_h__
#define __idxanalyzer_h__

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include "indexpb.h"

using namespace std;

#define off_t long long int

template <class T> class SigStack;
template <class T> class PatternStack;
class IdxSigEntryList;

//used to describe a single pattern that found
//This will be saved in the stack
class PatternUnit {
    public:
        vector<off_t> seq;
        int cnt; //count of repeatition
        
        PatternUnit() {}
        PatternUnit( vector<off_t> sq, int ct )
            :seq(sq),cnt(ct)
        {}
        //return number of elements in total
        int size() const 
        {
            if ( cnt == 0 ) {
                return 1; //not repetition
            } else {
                return seq.size()*cnt;
            }
        }
        
        int memsize() 
        {
            return sizeof(cnt) + sizeof(off_t)*seq.size();
        }

        virtual void show() const
        {
            vector<off_t>::const_iterator iter;
            cout << "( " ;
            for (iter = seq.begin();
                    iter != seq.end();
                    iter++ )
            {
                cout << *iter << " ";
            }
            cout << ") ^" << cnt << endl;
        }
};

//This is used to describ a single repeating
//pattern, but with starting value
class IdxSigUnit: public PatternUnit {
    public:
        off_t init; // the initial value of 
                    // logical offset, length, or physical offset
        void show() const
        {
            cout << init << " ... ";
            PatternUnit::show();
        }

        int memsize() {
            return sizeof(init) + PatternUnit::memsize();
        }
};

template <class T> // T can be PatternUnit or IdxSigUnit
class PatternStack {
    public:
        PatternStack() {}
        void push( T pu ) 
        {
            the_stack.push_back(pu);
        }
        
        void clear() 
        {
            the_stack.clear();
        }

        //if popping out t elem breaks any patterns
        bool isPopSafe( int t ) 
        {
            typename vector<T>::reverse_iterator rit;
            
            int total = 0;
            rit = the_stack.rbegin();
            while ( rit != the_stack.rend()
                    && total < t )
            {
                total += rit->size();
                rit++;
            }
            return total == t;
        }

        //return false if it is not safe
        //t is number of elements, not pattern unit
        bool popElem ( int t )
        {
            if ( !isPopSafe(t) ) {
                return false;
            }

            int total = 0; // the number of elem already popped out
            while ( !the_stack.empty() && total < t ) {
                total += top().size();
                the_stack.pop_back();
            }
            assert( total == t );

            return true;
        }

        //pop out one pattern
        void popPattern () 
        {
            the_stack.pop_back();
        }
        
        //make sure the stack is not empty before using this
        T top () 
        {
            assert( the_stack.size() > 0 );
            return the_stack.back();
        }

        typename vector<T>::const_iterator
            begin() const
        {
            return the_stack.begin();
        }
        
        typename vector<T>::const_iterator
            end() const
        {
            return the_stack.end();
        }
        
        virtual void show()
        {
            typename vector<T>::const_iterator iter;
            
            for ( iter = the_stack.begin();
                    iter != the_stack.end();
                    iter++ )
            {
                iter->show();
                /*
                vector<off_t>::const_iterator off_iter;
                for ( off_iter = (iter->seq).begin();
                        off_iter != (iter->seq).end();
                        off_iter++ )
                {
                    cout << *off_iter << ", ";
                }
                cout << "^" << iter->cnt << endl;
                */
            }
        }
    
    protected:
        vector<T> the_stack;
};

//I really should not call it a stack since I use
//it in many ways...
template <class T>
class SigStack: public PatternStack <T> 
{
    public:
        virtual void show()
        {
            typename vector<T>::const_iterator iter;

            for ( iter = this->the_stack.begin();
                    iter != this->the_stack.end();
                    iter++ )
            {
                vector<off_t>::const_iterator off_iter;
                cout << iter->init << "- " ;
                for ( off_iter = (iter->seq).begin();
                        off_iter != (iter->seq).end();
                        off_iter++ )
                {
                    cout << *off_iter << ", ";
                }
                cout << "^" << iter->cnt << endl;
            }
        }

        int memsize() 
        {
            typename vector<T>::const_iterator iter;
            int totalsize = 0;
            for ( iter = PatternStack<T>::the_stack.begin();
                    iter != PatternStack<T>::the_stack.end();
                    iter++)
            {
                totalsize += ((T)(*iter)).memsize();
            }
            return totalsize;
        }

};

class Tuple {
    public:
        int offset; //note that this is not the 
        // offset when accessing file. But
        // the offset in LZ77 algorithm
        int length; //concept in LZ77
        off_t next_symbol;

        Tuple() {}
        Tuple(int o, int l, off_t n) {
            offset = o;
            length = l;
            next_symbol = n;
        }

        void put(int o, int l, off_t n) {
            offset = o;
            length = l;
            next_symbol = n;
        }

        bool operator== (const Tuple other) {
            if (offset == other.offset 
                    && length == other.length
                    && next_symbol == other.next_symbol)
            {
                return true;
            } else {
                return false;
            }
        }

        // Tell if the repeating sequences are next to each other
        bool isRepeatingNeighbor() {
            return (offset == length && offset > 0);
        }

        void show() {
            cout << "(" << offset 
                << ", " << length
                << ", " << next_symbol << ")" << endl;
        }
};


class IdxEntry {
    public:
        int Proc;
#define ID_WRITE 0
#define ID_READ  1
        int ID; //either ID_WRITE or ID_READ
        off_t Logical_offset;
        off_t Length;
        double Begin_timestamp;
        double End_timestamp;
        off_t Logical_tail;
        int ID_2;
        off_t Physical_offset;
};


// Each index has its own signature
class IdxSignature {
    public:
        IdxSignature():win_size(6) {} 
        void discoverPattern( vector<off_t> const &seq );
        SigStack<IdxSigUnit> discoverSigPattern( vector<off_t> const &seq,
                vector<off_t> const &orig );
        //It takes in a entry buffer like in PLFS,
        //analyzes it and generate Index Signature Entries
        IdxSigEntryList generateIdxSignature(vector<IdxEntry> &entry_buf, int proc);
    private:
        vector<IdxEntry> entry_buf;
        int win_size; //window size
        Tuple searchNeighbor( vector<off_t> const &seq,
                vector<off_t>::const_iterator p_lookahead_win ); 
};

//This is the new entry for the new index
//file using I/O signature. It corresponds to
//HostEntry in old PLFS.
//
//Damn, where can I put the time stamp :(
class IdxSigEntry {
    public:
        int memsize() 
        {
            return sizeof(proc) + logical_offset.memsize() 
                + length.memsize() + physical_offset.memsize();
        }

    public:
        int proc;
        IdxSigUnit logical_offset;
        SigStack<IdxSigUnit> length;
        SigStack<IdxSigUnit> physical_offset;
};


class IdxSigEntryList {
    public:
        vector<IdxSigEntry> list;
        idxfile::EntryList pb_list;

    public:
        void append(IdxSigEntryList other);
        void append(vector<IdxSigEntry> &other);
        void show();
        //TODO:
        void saveToFile(char *filename);
        void readFromFile(char *filename);
        void siglistToPblist(vector<IdxSigEntry> &slist,
                idxfile::EntryList &pblist);
};

void printIdxEntries( vector<IdxSigEntry> &idx_entry_list );
vector<off_t> buildDeltas( vector<off_t> seq );
#endif

