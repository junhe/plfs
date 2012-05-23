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

template <class T> class SigStack;
template <class T> class PatternStack;
class IdxSigEntryList;
class HostEntry;

typedef int32_t header_t; //the type to hold the body size in serialization

void appendToBuffer( string &to, const void *from, const int size );

//Note that this function will increase start
void readFromBuf( string &from, void *to, int &start, const int size );



//used to describe a single pattern that found
//This will be saved in the stack
class PatternUnit {
    public:
        vector<off_t> seq;
        int64_t cnt; //count of repeatition
        
        PatternUnit() {}
        PatternUnit( vector<off_t> sq, int ct )
            :seq(sq),cnt(ct)
        {}
        
        int size() const;
        virtual void show() const;
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

        header_t bodySize();
        string serialize();
        void deSerialize(string buf);
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
        int bodySize();
        string serialize();    
        void deSerialize( string buf );
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

/* Not longer in use
 *
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
*/

// this is the class that represents the records that get written into the
// index file for each host.
// TODO:
// Move this class somewhere proper. It is impelmented in Index.cpp
class HostEntry
{
    public:
        HostEntry();
        HostEntry( off_t o, size_t s, pid_t p );
        HostEntry( const HostEntry& copy );
        bool overlap( const HostEntry& );
        bool contains ( off_t ) const;
        bool splittable ( off_t ) const;
        bool abut   ( const HostEntry& );
        off_t logical_tail( ) const;
        bool follows(const HostEntry&);
        bool preceeds(const HostEntry&);

    protected:
        off_t  logical_offset;
        off_t  physical_offset;  // I tried so hard to not put this in here
        // to save some bytes in the index entries
        // on disk.  But truncate breaks it all.
        // we assume that each write makes one entry
        // in the data file and one entry in the index
        // file.  But when we remove entries and
        // rewrite the index, then we break this
        // assumption.  blech.
        size_t length;
        double begin_timestamp;
        double end_timestamp;
        pid_t  id;      // needs to be last so no padding

        friend class Index;
        friend class IdxSignature;
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
        IdxSigEntryList generateIdxSignature(vector<HostEntry> &entry_buf, int proc);
    private:
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
        pid_t original_chunk;  //used only when entry is in global
                               //complex index. 
        pid_t new_chunk_id;    //This is not serialized yet.
                               //it should only be serialized in
                               //the context of global complex index
        IdxSigUnit logical_offset;
        SigStack<IdxSigUnit> length;
        SigStack<IdxSigUnit> physical_offset;
        string serialize();
        void deSerialize(string buf);
        int bodySize();
};


class IdxSigEntryList {
    public:
        vector<IdxSigEntry> list;
        idxfile::EntryList pb_list;

    public:
        void append(IdxSigEntryList other);
        void append(vector<IdxSigEntry> &other);
        void show();
        void saveToFile(const char *filename);
        void saveToFile(const int fd);
        void readFromFile(char *filename);
        void siglistToPblist(vector<IdxSigEntry> &slist,
                idxfile::EntryList &pblist);
        void siglistToPblist();
        void clear();
        string serialize();
        void deSerialize(string buf);
        int bodySize();
};

void printIdxEntries( vector<IdxSigEntry> &idx_entry_list );
vector<off_t> buildDeltas( vector<off_t> seq );

template <class T>
string 
PatternStack<T>::serialize()
{
    header_t bodysize = 0;
    string buf;
    typename vector<T>::iterator iter;
    header_t realtotalsize = 0;

    bodysize = bodySize();
    //cout << "data size put in: " << bodysize << endl;
    appendToBuffer(buf, &bodysize, sizeof(bodysize));
    for ( iter = the_stack.begin() ; 
            iter != the_stack.end() ;
            iter++ )
    {
        string unit = iter->serialize();
        appendToBuffer(buf, &unit[0], unit.size());
        realtotalsize += unit.size();
        //to test if it serilized right
        //IdxSigUnit tmp;
        //tmp.deSerialize(unit);
        //cout << "test show.\n";
        //tmp.show();
    }
    assert(realtotalsize == bodysize);
    return buf;
}

template <class T>
void
PatternStack<T>::deSerialize( string buf )
{
    header_t bodysize, bufsize;
    int cur_start = 0;
    
    clear(); 

    readFromBuf(buf, &bodysize, cur_start, sizeof(bodysize));
    
    bufsize = buf.size();
    assert(bufsize == bodysize + sizeof(bodysize));
    while ( cur_start < bodysize ) {
        header_t unitbodysize;
        string unitbuf;
        T sigunit;

        readFromBuf(buf, &unitbodysize, cur_start, sizeof(unitbodysize));
        //cout << "Unitbodysize:" << unitbodysize << endl;
        int sizeofheadandunit = sizeof(unitbodysize) + unitbodysize;
        unitbuf.resize(sizeofheadandunit);
        if ( unitbodysize > 0 ) {
            //TODO:make this more delegate
            cur_start -= sizeof(unitbodysize);
            readFromBuf(buf, &unitbuf[0], cur_start, sizeofheadandunit); 
        }
        sigunit.deSerialize(unitbuf);
        push(sigunit);
    }

}

template <class T>
int
PatternStack<T>::bodySize()
{
    int totalsize = 0;
    typename vector<T>::iterator iter;
    for ( iter = the_stack.begin() ; 
            iter != the_stack.end() ;
            iter++ )
    {
        //IdxSigUnit body size and its header
        totalsize += (iter->bodySize() + sizeof(header_t));
    }
    return totalsize;
}

#endif

