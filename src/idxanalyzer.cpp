#include "idxanalyzer.h"
#include "Util.h"
#include "plfs_private.h"

#include <algorithm>

//for debugging
string printIdxEntries( vector<IdxSigEntry> &idx_entry_list )
{
    vector<IdxSigEntry>::iterator iter;

    //cout << "this is printIdxEntries" << endl;

    ostringstream showstr;
    for ( iter = idx_entry_list.begin();
            iter != idx_entry_list.end();
            iter++ )
    {
        showstr << iter->show();
    }
    return showstr.str();
}

string printVector( vector<off_t> vec )
{
    ostringstream oss;
    oss << "PrintVector: " ;
    vector<off_t>::const_iterator it;
    for ( it =  vec.begin() ;
          it != vec.end() ;
          it++ )
    {
        oss << *it << "," ;
    }
    oss << endl;
    return oss.str();
}

vector<off_t> buildDeltas( vector<off_t> seq ) 
{
    vector<off_t>::iterator it;
    vector<off_t> deltas;
    //ostringstream oss;
    for ( it = seq.begin() ; it != seq.end() ; it++ )
    {
        if ( it > seq.begin() ) {
            deltas.push_back( *it - *(it-1) );
            //oss << (*it - *(it-1)) << "-";
        }
    }
    //mlog(IDX_WARN, "deltas:%s", oss.str().c_str());
    //cout << "in builddeltas: " << seq.size() << " " << deltas.size() << endl; 
    return deltas;
}

//IdxEntry is designed to mimic the HostEntry in PLFS,
//so later I can easily copy and paste these code to PLFS
//and make it work.
//
//It gets signatures for a proc. 
//TODO:
//But do we really need to separate entries by original_chunk at first?
//Also, have to handle the case that there is only one entry
IdxSigEntryList IdxSignature::generateIdxSignature(
        vector<HostEntry> &entry_buf, 
        int proc) 
{
    vector<off_t> logical_offset, length, physical_offset; 
    vector<off_t> logical_offset_delta, 
                    length_delta, 
                    physical_offset_delta;
    IdxSigEntryList entrylist;
    vector<HostEntry>::const_iterator iter;
    
    for ( iter = entry_buf.begin() ; 
            iter != entry_buf.end() ;
            iter++ )
    {
        if ( iter->id != proc ) {
            continue;
        }

        logical_offset.push_back(iter->logical_offset);
        length.push_back(iter->length);
        physical_offset.push_back(iter->physical_offset);
    }
    
    if ( !(logical_offset.size() == length.size() &&
            length.size() == physical_offset.size()) ) {
        ostringstream oss;
        oss << "logical_offset.size():" << logical_offset.size() 
            << "length.size():" << length.size()
            << "physical_offset.size():" << physical_offset.size() << endl;
        mlog(IDX_WARN, "sizes should be equal. %s", oss.str().c_str());
        exit(-1);
    }

    SigStack<IdxSigUnit> offset_sig = generateComplexPatterns(logical_offset);
    //offset_sig.show();
    //cout << "Just showed offset_sig" << endl;
    //Now, go through offset_sig one by one and build the IdxSigEntry s
    vector<IdxSigEntry>idx_entry_list;
    vector<IdxSigUnit>::const_iterator stack_iter;

    ostringstream oss;
    oss << "offset patterns:" << offset_sig.show() << endl;
    oss << "offsets: " << printVector( logical_offset ) << endl;
    mlog(IDX_WARN, "%s", oss.str().c_str());

    int range_start = 0, range_end; //the range currently processing
    for (stack_iter = offset_sig.begin();
            stack_iter != offset_sig.end();
            stack_iter++ )
    {
        //cout << stack_iter->init << " " ;
        IdxSigEntry idx_entry;
        range_end = range_start 
                    + max( int(stack_iter->cnt * stack_iter->seq.size()), 1 );
        mlog(IDX_WARN, "range_start:%d, range_end:%d, logical_offset.size():%d",
                range_start, range_end, logical_offset.size());
        assert( range_end <= logical_offset.size() );

        idx_entry.original_chunk = proc;
        idx_entry.logical_offset = *stack_iter;
        
        vector<off_t> length_seg(length.begin()+range_start,
                                 length.begin()+range_end);
        idx_entry.length = generateComplexPatterns(length_seg);

        vector<off_t> physical_offset_seg( physical_offset.begin()+range_start,
                                           physical_offset.begin()+range_end);
        idx_entry.physical_offset = generateComplexPatterns(physical_offset_seg);
        
        idx_entry_list.push_back( idx_entry);

        range_start = range_end;
    }
    entrylist.append(idx_entry_list);
    //printIdxEntries(idx_entry_list);
    return entrylist;
}



//find out pattern of a number sequence 
void IdxSignature::discoverPattern(  vector<off_t> const &seq )
{
    vector<off_t>::const_iterator p_lookahead_win; // pointer(iterator) to the lookahead window
    PatternStack<PatternUnit> pattern_stack;

    p_lookahead_win = seq.begin();
    pattern_stack.clear();

    //cout << endl << "this is discoverPattern() :)" << endl;

    while ( p_lookahead_win != seq.end() ) {
        //lookahead window is not empty
        Tuple cur_tuple = searchNeighbor( seq, p_lookahead_win );
        //cur_tuple.show();
        if ( cur_tuple.isRepeatingNeighbor() ) {
            if ( pattern_stack.isPopSafe( cur_tuple.length ) ) {
                //safe
                pattern_stack.popElem( cur_tuple.length );

                vector<off_t>::const_iterator first, last;
                first = p_lookahead_win;
                last = p_lookahead_win + cur_tuple.length;

                PatternUnit pu;
                pu.seq.assign(first, last);
                pu.cnt = 2;

                pattern_stack.push( pu );
                p_lookahead_win += cur_tuple.length;
            } else {
                //unsafe
                PatternUnit pu = pattern_stack.top();

                if ( pu.seq.size() == cur_tuple.length ) {
                    //the subseq in lookahead window repeats
                    //the top pattern in stack
                    pu.cnt++;
                    pattern_stack.popPattern();
                    pattern_stack.push(pu);

                    p_lookahead_win += cur_tuple.length;
                } else {
                    //cannot pop out cur_tuple.length elems without
                    //totally breaking any pattern.
                    //So we simply add one elem to the stack
                    PatternUnit pu2;
                    pu2.seq.push_back( *p_lookahead_win );
                    pu2.cnt = 1;
                    pattern_stack.push(pu2);
                    p_lookahead_win++;
                }
            }


        } else {
            //(0,0,x)
            PatternUnit pu;
            pu.seq.push_back(cur_tuple.next_symbol);
            pu.cnt = 1;

            pattern_stack.push(pu); 
            p_lookahead_win++;
        }
        //pattern_stack.show();
    }

}

SigStack<IdxSigUnit> IdxSignature::generateComplexPatterns( vector<off_t> inits )
{
    mlog(IDX_WARN, "Entering %s. inits: %s", __FUNCTION__, printVector(inits).c_str());
    vector<off_t> deltas = buildDeltas(inits);
    SigStack<IdxSigUnit> pattern;

    pattern = findPattern( deltas );
    int pos = 0;
    vector<IdxSigUnit>::iterator it;
    for ( it = pattern.the_stack.begin() ;
          it != pattern.the_stack.end() ;
          it++ )
    {
        it->init = inits[pos];
        pos += it->seq.size() * it->cnt;
    }
    // handle the missing last
    if ( ! inits.empty() && ! pattern.the_stack.empty() ) {
        IdxSigUnit last = pattern.the_stack.back(); //for short
        if ( last.seq.size() == 1 && 
             last.init + last.seq[0]*last.cnt  == inits.back() ) 
        {
            pattern.the_stack.back().cnt++;
        } else {        
            IdxSigUnit punit;
            punit.init = inits.back();
            punit.cnt = 1;
            pattern.push(punit);
        }
    }    
    
    ostringstream oss;
    oss << pattern.show() << endl; 
    mlog(IDX_WARN, "Leaving %s:{%s}", __FUNCTION__, oss.str().c_str());

    return pattern;
}



SigStack<IdxSigUnit> IdxSignature::findPattern( vector<off_t> deltas )
{
    // pointer(iterator) to the lookahead window, bot should move together
    mlog(IDX_WARN, "Entering %s", __FUNCTION__);
    mlog(IDX_WARN, "deltas: %s", printVector(deltas).c_str());
    vector<off_t>::const_iterator p_lookahead_win;
    SigStack<IdxSigUnit> pattern_stack;

    p_lookahead_win = deltas.begin();
    pattern_stack.clear();

    while ( p_lookahead_win < deltas.end() ) {
        //lookahead window is not empty
        mlog(IDX_WARN, "window position:%d", p_lookahead_win - deltas.begin());
        Tuple cur_tuple = searchNeighbor( deltas, p_lookahead_win );
        mlog(IDX_WARN, "%s", cur_tuple.show().c_str());
        if ( cur_tuple.isRepeatingNeighbor() ) {
            if ( pattern_stack.isPopSafe( cur_tuple.length ) ) {
                mlog(IDX_WARN, "SAFE" );
                //safe
                pattern_stack.popElem( cur_tuple.length );

                vector<off_t>::const_iterator first, last;
                first = p_lookahead_win;
                last = p_lookahead_win + cur_tuple.length;

                IdxSigUnit pu;
                pu.seq.assign(first, last);
                pu.cnt = 2;

                pattern_stack.push( pu );
                p_lookahead_win += cur_tuple.length;
            } else {
                //unsafe
                mlog(IDX_WARN, "UNSAFE" );
                IdxSigUnit pu = pattern_stack.top();

                if ( cur_tuple.length % pu.seq.size() == 0
                     && cur_tuple.length <= pu.seq.size() * pu.cnt ) {
                    mlog(IDX_WARN, "-REPEATING LAST Pattern");
                    //the subseq in lookahead window repeats
                    //the top pattern in stack.
                    //initial remains the same.
                    pu.cnt += cur_tuple.length / pu.seq.size() ;
                    pattern_stack.popPattern();
                    pattern_stack.push(pu);
                    //tmp debug.

                    p_lookahead_win += cur_tuple.length;
                } else {
                    mlog(IDX_WARN, "-Just Pust it in");
                    //cannot pop out cur_tuple.length elems without
                    //totally breaking any pattern.
                    //So we simply add one elem to the stack
                    IdxSigUnit pu;
                    pu.seq.push_back( *p_lookahead_win );
                    pu.cnt = 1;
                    pattern_stack.push(pu);
                    p_lookahead_win++;
                }
            }
        } else {
            //(0,0,x)
            IdxSigUnit pu;
            pu.seq.push_back(cur_tuple.next_symbol);
            pu.cnt = 1;

            pattern_stack.push(pu); 
            p_lookahead_win++;
        }
        pattern_stack.the_stack.back().compressRepeats();
        mlog(IDX_WARN, "LOOP: %s", pattern_stack.show().c_str());
    }
   
    mlog(IDX_WARN, "Leaving %s:%s", __FUNCTION__, pattern_stack.show().c_str());
    return pattern_stack;

}

//find out pattern of a number sequence(deltas) with its
//original sequence
//if seq and orig have the same sizes.
//  the function returns pattern representing all orig numbers.
//else if orig has one more than seq (including seq.size()==0)
//  the function returns pattern representing all orig numbers, 
//  with the last orig num with seq.size()==0
//else
//  error
SigStack<IdxSigUnit> IdxSignature::discoverSigPattern( vector<off_t> const &seq,
        vector<off_t> const &orig )
{
    // pointer(iterator) to the lookahead window, bot should move together
    vector<off_t>::const_iterator p_lookahead_win, 
        p_lookahead_win_orig; 
    SigStack<IdxSigUnit> pattern_stack;

    p_lookahead_win = seq.begin();
    p_lookahead_win_orig = orig.begin();
    pattern_stack.clear();

    if (! (seq.size() == orig.size()
            || seq.size() + 1 == orig.size() ) )
    {
        ostringstream oss;
        oss << "seq.size():" << seq.size()
            << " orig.size():" << orig.size() << endl;
        mlog(IDX_ERR, "discoverSigPattern() needs to be used with "
                "seq.size==orig.size or seq.size+1==orig.size. \n %s", 
                oss.str().c_str() );
        exit(-1);
    }

    //cout << endl << "this is discoverPattern() :)" << endl;
   
    //TODO:
    //There's bug in handling the case of only one entry.
    //And whether 1,(3,4)^2 represnets five or four numbers.
    //Go back to handle this when protoc buffer is integrated.
    /*
    cout << "seq.size(): " << seq.size() << "orig.size():" << orig.size() << endl;
    assert(seq.size() == (orig.size()-1));

    //for the case there is only one entry
    if ( seq.size() == 0 ) {
        IdxSigUnit pu;
        pu.init = *p_lookahead_win_orig;
        pu.cnt = 0;
    }
    */    



    //TODO: WHY p_lookahead_win != seq.end() is a dead loop????
    while ( p_lookahead_win < seq.end() ) {
        //lookahead window is not empty
        Tuple cur_tuple = searchNeighbor( seq, p_lookahead_win );
        //cur_tuple.show();
        if ( cur_tuple.isRepeatingNeighbor() ) {
            if ( pattern_stack.isPopSafe( cur_tuple.length ) ) {
                //safe
                pattern_stack.popElem( cur_tuple.length );

                vector<off_t>::const_iterator first, last;
                first = p_lookahead_win;
                last = p_lookahead_win + cur_tuple.length;

                IdxSigUnit pu;
                pu.seq.assign(first, last);
                pu.cnt = 2;
                pu.init = *(p_lookahead_win_orig - cur_tuple.length);

                pattern_stack.push( pu );
                p_lookahead_win += cur_tuple.length;
                p_lookahead_win_orig += cur_tuple.length;
            } else {
                //unsafe
                IdxSigUnit pu = pattern_stack.top();

                if ( cur_tuple.length % pu.seq.size()
                     && cur_tuple.length <= pu.seq.size() * pu.cnt ) {
                    //the subseq in lookahead window repeats
                    //the top pattern in stack.
                    //initial remains the same.
                    pu.cnt += cur_tuple.length / pu.seq.size() ;
                    pattern_stack.popPattern();
                    pattern_stack.push(pu);
                    pu.init = *p_lookahead_win_orig; //should delete this. keep if only for 
                    //tmp debug.

                    p_lookahead_win += cur_tuple.length;
                    p_lookahead_win_orig += cur_tuple.length;
                } else {
                    //cannot pop out cur_tuple.length elems without
                    //totally breaking any pattern.
                    //So we simply add one elem to the stack
                    IdxSigUnit pu;
                    pu.seq.push_back( *p_lookahead_win );
                    pu.init = *p_lookahead_win_orig;
                    pu.cnt = 1;
                    pattern_stack.push(pu);
                    p_lookahead_win++;
                    p_lookahead_win_orig++;
                }
            }
        } else {
            //(0,0,x)
            IdxSigUnit pu;
            pu.seq.push_back(cur_tuple.next_symbol);
            pu.init = *p_lookahead_win_orig;
            pu.cnt = 1;

            pattern_stack.push(pu); 
            p_lookahead_win++;
            p_lookahead_win_orig++;
        }
    }
   
    if ( p_lookahead_win_orig < orig.end() ) {
        assert(p_lookahead_win_orig + 1 == orig.end());
        IdxSigUnit pu;
        pu.init = *p_lookahead_win_orig;
        pu.cnt = 0;

        pattern_stack.push(pu); 
    }
   
    SigStack<IdxSigUnit> pattern_stack_compressed;
    vector<IdxSigUnit>::iterator it;
    for ( it = pattern_stack.the_stack.begin();
          it != pattern_stack.the_stack.end();
          it++ )
    {
        it->compressRepeats();
        if (pattern_stack_compressed.the_stack.empty()) {
            mlog(IDX_WARN, "Empty");
            pattern_stack_compressed.the_stack.push_back(*it);
        } else {
            bool ret;
            ret = pattern_stack_compressed.the_stack.back().append(*it);
            if (ret == false) {
                pattern_stack_compressed.the_stack.push_back(*it);
            }
        }
        //ostringstream oss;
        //oss << pattern_stack_compressed.show();
        //mlog(IDX_WARN, "%s", oss.str().c_str());
    }
   
    

    if ( pattern_stack.size() != orig.size() ) {
        ostringstream oss;
        oss<< "pattern_stack.size() != orig.size() in"
               << __FUNCTION__ << pattern_stack.size() 
               << "," << orig.size() << endl;
        oss << "seq.size():" << seq.size() << endl;
        oss << pattern_stack.show() << endl;
        vector<off_t>::const_iterator it;
        for ( it = orig.begin();
              it != orig.end();
              it++ )
        {
            oss << *it << ",";
        }
        oss << endl;
        mlog(IDX_ERR, "%s", oss.str().c_str());
        exit(-1);
    }

    return pattern_stack_compressed;
}

Tuple IdxSignature::searchNeighbor( vector<off_t> const &seq,
        vector<off_t>::const_iterator p_lookahead_win ) 
{
    vector<off_t>::const_iterator i;     
    int j;
    //cout << "------------------- I am in searchNeighbor() " << endl;

    //i goes left util the begin or reaching window size
    i = p_lookahead_win;
    int remain = seq.end() - p_lookahead_win;
    while ( i != seq.begin() 
            && (p_lookahead_win - i) < win_size
            && (p_lookahead_win - i) < remain ) {
        i--;
    }
    //termination: i == seq.begin() or distance == win_size

    /*
    //print out search buffer and lookahead buffer
    //cout << "search buf: " ;
    vector<off_t>::const_iterator k;
    for ( k = i ; k != p_lookahead_win ; k++ ) {
    cout << *k << " ";
    }
    cout << endl;

    cout << "lookahead buf: " ;
    vector<off_t>::const_iterator p;
    p = p_lookahead_win;
    for ( p = p_lookahead_win ; 
    p != seq.end() && p - p_lookahead_win < win_size ; 
    p++ ) {
    cout << *p << " ";
    }
    cout << endl;
    */

    //i points to a element in search buffer where matching starts
    //j is the iterator from the start to the end of search buffer to compare
    for ( ; i != p_lookahead_win ; i++ ) {
        int search_length = p_lookahead_win - i;
        for ( j = 0 ; j < search_length ; j++ ) {
            if ( *(i+j) != *(p_lookahead_win + j) ) {
                break;
            }
        }
        if ( j == search_length ) {
            //found a repeating neighbor
            return Tuple(search_length, search_length, 
                    *(p_lookahead_win + search_length));
        }
    }

    //Cannot find a repeating neighbor
    return Tuple(0, 0, *(p_lookahead_win));
}

void IdxSigEntryList::append( vector<IdxSigEntry> &other ) 
{
    vector<IdxSigEntry>::iterator iter;
    for (iter = other.begin();
            iter != other.end();
            iter++ )
    {
        append(*iter, true);
    }

}

void IdxSigEntryList::append( IdxSigEntryList other ) 
{
    append(other.list);
}

void IdxSigEntryList::append( IdxSigEntry other, bool compress ) 
{
    if ( compress == false || list.empty() ) {
        list.push_back(other);
        return ;
    } else {
        if ( ! list.back().append(other) ) {
            list.push_back(other);
            return;
        }
    }
    return;
}

string 
IdxSigEntryList::show() 
{
    ostringstream showstr;
    showstr << printIdxEntries(list);

    vector<HostEntry>::const_iterator it;
    for ( it = messies.begin() ;
          it != messies.end() ;
          it++ ) 
    {
        double begin_timestamp = 0, end_timestamp = 0;
        begin_timestamp = (*it).begin_timestamp;
        end_timestamp  = (*it).end_timestamp;
        showstr  << setw(5)
            << (*it).id             << " w "
            << setw(16)
            << (*it).logical_offset << " "
            << setw(8) << (*it).length << " "
            << setw(16) << fixed << setprecision(16)
            << begin_timestamp << " "
            << setw(16) << fixed << setprecision(16)
            << end_timestamp   << " "
            << setw(16)
            << (*it).logical_tail() << " "
            << " [" << setw(10) << (*it).physical_offset << "]";
    }

    return showstr.str();
}

// [bodysize][type][body:[entry0][entry1]...]
void IdxSigEntryList::saveMessiesToFile(const int fd)
{
    header_t entrybodysize = sizeof(HostEntry)*messies.size();
    if ( entrybodysize == 0 ) {
        return;

    }
    char entrytype = 'M'; //means Messies

    Util::Writen(fd, &entrybodysize, sizeof(entrybodysize));
    Util::Writen(fd, &entrytype, sizeof(entrytype));

    if ( entrybodysize > 0 ) {
        Util::Writen(fd, &messies[0], entrybodysize);
    }
    return ;
}

// It appends to the buffer, to save time on copying
string &IdxSigEntryList::serializeMessies( string &buf )
{
    header_t entrybodysize = sizeof(HostEntry)*messies.size();
    char entrytype = 'M'; //means Messies

    appendToBuffer(buf, &entrybodysize, sizeof(entrybodysize));
    appendToBuffer(buf, &entrytype, sizeof(entrytype));

    if ( entrybodysize > 0 ) {
        appendToBuffer(buf, &messies[0], entrybodysize);
    }   
    return buf;
}


void IdxSigEntryList::deSerializeMessies( string &buf )
{
    header_t entrybodysize;
    char entrytype; //means Messies
    int cur = 0;

    readFromBuf(buf, &entrybodysize, cur, sizeof(entrybodysize));
    readFromBuf(buf, &entrytype, cur, sizeof(entrytype));
    messies.resize( entrybodysize/sizeof(HostEntry) );

    if ( entrybodysize > 0 ) {
        readFromBuf(buf, &messies[0], cur, entrybodysize);
    }
    return;
}


void IdxSigEntryList::saveListToFile(const int fd)
{
    string buf = serialize();
    if ( buf.size() > 0 ) {
        Util::Writen(fd, &buf[0], buf.size());
    }
}



void IdxSigEntryList::saveToFile(const int fd)
{
    saveListToFile(fd);
    saveMessiesToFile(fd);
}


void IdxSigEntryList::clear()
{
    list.clear();
    messies.clear();
}


void appendToBuffer( string &to, const void *from, const int size )
{
    if ( size > 0 ) { //make it safe
        to.append( (char *)from, size );
    }
}

//Note that this function will increase start
void readFromBuf( string &from, void *to, int &start, const int size )
{
    //'to' has to be treated as plain memory
    memcpy(to, &from[start], size);
    start += size;
}

//Serialiezd IdxSigUnit: [head:bodysize][body]
header_t IdxSigUnit::bodySize()
{
    header_t totalsize;
    totalsize = sizeof(init) //init
                + sizeof(cnt) //cnt
                + sizeof(header_t) //length of seq size header
                + seq.size()*sizeof(off_t);
    return totalsize;
}

// check if this follows other and merge
// has to satisfy two:
// 1. seq are exactly the same OR (repeating and the same, size can be diff)
//    (3,3,3)==(3,3,3)             (3,3,3)==(3), (3,3,3)==()
// 2. AND init1 + sum of deltas == init2
// return true if appended successfully
bool IdxSigUnit::append( IdxSigUnit &other )
{

    mlog(IDX_WARN, "in %s", __FUNCTION__);

    if ( this->isSeqRepeating() 
        && other.isSeqRepeating() )
    {
        if ( this->size() > 1 && other.size() > 1 ) {
            //case 1. both has size > 1
            if ( this->seq[0] == other.seq[0] 
                 && this->init + this->seq[0]*this->size() == other.init ) {
                int newsize = this->size() + other.size();
                this->seq.clear();
                this->seq.push_back(other.seq[0]);
                this->cnt = newsize;
                return true;
            } else {
                return false;
            }               
        } else if ( this->size() == 1 && other.size() == 1 ) {
            //case 2. both has size == 1
            //definitely follows
            this->seq.clear();
            this->seq.push_back(other.init - this->init);
            this->cnt = 2; //has two now
            return true;
        } else if ( this->size() == 1 && other.size() > 1 ) {
            if ( other.init - this->init == other.seq[0] ) {
                int newsize = this->size() + other.size();
                this->seq.clear();
                this->seq.push_back(other.seq[0]);
                this->cnt = newsize;
                return true;
            } else {
                return false;
            }
        } else if ( this->size() > 1 && other.size() == 1) {
            if ( this->init + this->seq[0]*this->size() == other.init ) {
                int newsize = this->size() + other.size();
                off_t tmp = this->seq[0];
                this->seq.clear();
                this->seq.push_back(tmp);
                this->cnt = newsize;
                return true;
            } else {
                return false;
            }
        }
    } else {
        return false;  //TODO:should handle this case
    }
}

// (3,3,3)^4 is repeating
// (0,0)^0 is also repeating
bool IdxSigUnit::isSeqRepeating()
{
    vector<off_t>::iterator it;
    bool allrepeat = true;
    for ( it = seq.begin();
          it != seq.end();
          it++ )
    {
        if ( it != seq.begin()
             && *it != *(it-1) ) {
            allrepeat = false;
            break;
        }
    }
    return allrepeat;
}

void IdxSigUnit::compressRepeats()
{
    if ( isSeqRepeating() && size() > 1 ) {
        cnt = size();
        off_t tmp = seq[0];
        seq.clear();
        seq.push_back(tmp);
    }
}


string 
IdxSigUnit::serialize()
{
    string buf; //let me put it in string and see if it works
    header_t seqbodysize;
    header_t totalsize;

    totalsize = bodySize(); 
    
    appendToBuffer(buf, &totalsize, sizeof(totalsize));
    appendToBuffer(buf, &init, sizeof(init));
    appendToBuffer(buf, &cnt, sizeof(cnt));
    seqbodysize = seq.size()*sizeof(off_t);
    appendToBuffer(buf, &(seqbodysize), sizeof(header_t));
    if (seqbodysize > 0 ) {
        appendToBuffer(buf, &seq[0], seqbodysize);
    }
    return buf;
}

//This input buf should be [data size of the followed data][data]
void 
IdxSigUnit::deSerialize(string buf)
{
    header_t totalsize;
    int cur_start = 0;
    header_t seqbodysize;

    readFromBuf(buf, &totalsize, cur_start, sizeof(totalsize));
    readFromBuf(buf, &init, cur_start, sizeof(init));
    readFromBuf(buf, &cnt, cur_start, sizeof(cnt));
    readFromBuf(buf, &seqbodysize, cur_start, sizeof(header_t));
    if ( seqbodysize > 0 ) {
        seq.resize(seqbodysize/sizeof(off_t));
        readFromBuf(buf, &seq[0], cur_start, seqbodysize); 
    }
}

//byte size is in [bodysize][data]
//it is the size of data
int IdxSigEntry::bodySize()
{
    int totalsize = 0;
    totalsize += sizeof(original_chunk);
    totalsize += sizeof(new_chunk_id);
    totalsize += sizeof(header_t) * 3; //the header size of the following 
    totalsize += logical_offset.bodySize();
    totalsize += length.bodySize();
    totalsize += physical_offset.bodySize();

    return totalsize;
}

string IdxSigEntry::serialize()
{
    header_t totalsize = 0;
    string buf, tmpbuf;
    header_t datasize;
    
    totalsize = bodySize();
    //cout << "IdxSigEntry totalsize put in: " << totalsize << endl;
    appendToBuffer(buf, &totalsize, sizeof(totalsize));
    appendToBuffer(buf, &original_chunk, sizeof(original_chunk));
    appendToBuffer(buf, &new_chunk_id, sizeof(new_chunk_id));
    //cout << "IdxSigEntry original_chunk put in: " << original_chunk << endl; 
    
    //this tmpbuf includes [data size][data]
    tmpbuf = logical_offset.serialize(); 
    appendToBuffer(buf, &tmpbuf[0], tmpbuf.size());
    
    tmpbuf = length.serialize();
    appendToBuffer(buf, &tmpbuf[0], tmpbuf.size());

    tmpbuf = physical_offset.serialize();
    appendToBuffer(buf, &tmpbuf[0], tmpbuf.size());

    return buf;
}



void IdxSigEntry::deSerialize(string buf)
{
    header_t totalsize = 0; 
    int cur_start = 0;
    header_t datasize = 0;
    string tmpbuf;


    readFromBuf(buf, &totalsize, cur_start, sizeof(totalsize));
    //cout << "IdxSigEntry totalsize read out: " << totalsize << endl;
    
    readFromBuf(buf, &original_chunk, cur_start, sizeof(original_chunk));
    //cout << "IdxSigEntry id read out: " << id << endl; 
    
    readFromBuf(buf, &new_chunk_id, cur_start, sizeof(new_chunk_id));
   
    tmpbuf.clear();
    readFromBuf(buf, &datasize, cur_start, sizeof(datasize));
    if ( datasize > 0 ) {
        int headanddatasize = sizeof(datasize) + datasize;
        tmpbuf.resize(headanddatasize);
        cur_start -= sizeof(datasize);
        readFromBuf(buf, &tmpbuf[0], cur_start, headanddatasize); 
    }
    logical_offset.deSerialize(tmpbuf);
    //cout << "deSerialized logical offset data size: " << datasize << endl;
    
    tmpbuf.clear();
    readFromBuf(buf, &datasize, cur_start, sizeof(datasize));
    if ( datasize > 0 ) {
        int headanddatasize = sizeof(datasize) + datasize;
        tmpbuf.resize(headanddatasize);
        cur_start -= sizeof(datasize);
        readFromBuf(buf, &tmpbuf[0], cur_start, headanddatasize); 
    }
    length.deSerialize(tmpbuf);

    tmpbuf.clear();
    readFromBuf(buf, &datasize, cur_start, sizeof(datasize));
    if ( datasize > 0 ) {
        int headanddatasize = sizeof(datasize) + datasize;
        tmpbuf.resize(headanddatasize);
        cur_start -= sizeof(datasize);
        readFromBuf(buf, &tmpbuf[0], cur_start, headanddatasize); 
    }
    physical_offset.deSerialize(tmpbuf);
}

string IdxSigEntryList::serialize()
{
    header_t bodysize, realbodysize = 0;
    char entrytype = 'P'; // means Pattern
    string buf;
    vector<IdxSigEntry>::iterator iter;
    
    bodysize = bodySize();

    appendToBuffer(buf, &bodysize, sizeof(bodysize));
    appendToBuffer(buf, &entrytype, sizeof(entrytype));
    
    //cout << "list body put in: " << bodysize << endl;

    for ( iter = list.begin() ;
          iter != list.end() ;
          iter++ )
    {
        string tmpbuf;
        tmpbuf = iter->serialize();
        if ( tmpbuf.size() > 0 ) {
            appendToBuffer(buf, &tmpbuf[0], tmpbuf.size());
        }
        realbodysize += tmpbuf.size();
    }
    assert(realbodysize == bodysize);
    //cout << realbodysize << "==" << bodysize << endl;

    return buf;
}

void IdxSigEntryList::deSerialize(string buf)
{
    header_t bodysize, bufsize;
    char entrytype;
    int cur_start = 0;

    list.clear();
    readFromBuf(buf, &bodysize, cur_start, sizeof(bodysize));
    readFromBuf(buf, &entrytype, cur_start, sizeof(entrytype));
    assert(entrytype == 'P'); 
    bufsize = buf.size();
    assert(bufsize == bodysize + sizeof(bodysize) + sizeof(entrytype));
    while ( cur_start < bufsize ) {
        header_t unitbodysize, sizeofheadandbody;
        string unitbuf;
        IdxSigEntry unitentry;

        readFromBuf(buf, &unitbodysize, cur_start, sizeof(unitbodysize));
        sizeofheadandbody = sizeof(unitbodysize) + unitbodysize; 
        unitbuf.resize(sizeofheadandbody);
        if ( unitbodysize > 0 ) {
            cur_start -= sizeof(unitbodysize);
            readFromBuf(buf, &unitbuf[0], cur_start, sizeofheadandbody);
        }
        unitentry.deSerialize(unitbuf);
        list.push_back(unitentry); //it is OK to push a empty entry
    }
    assert(cur_start==bufsize);
}

int IdxSigEntryList::bodySize()
{
    int bodysize = 0;
    vector<IdxSigEntry>::iterator iter;
    
    for ( iter = list.begin() ;
          iter != list.end() ;
          iter++ )
    {
        bodysize += iter->bodySize() + sizeof(header_t);
    }

    return bodysize;
}

void IdxSigEntryList::dumpMessies()
{
    vector<IdxSigEntry> nicepatterns;
    vector<IdxSigEntry>::iterator it;

    for ( it =  list.begin() ;
          it != list.end() ;
          it++ )
    {
        IdxSigEntry &patternentry = *it;
        if ( patternentry.logical_offset.size() <= 1 ) {
            //this is mess that had not been compressed
            HostEntry hentry;
            hentry.logical_offset = patternentry.logical_offset.init;
            assert( patternentry.length.the_stack.size() > 0 );
            hentry.length = patternentry.length.the_stack.front().init;
            assert( patternentry.physical_offset.the_stack.size() > 0 );
            hentry.physical_offset = patternentry.physical_offset.the_stack.front().init;
            hentry.begin_timestamp = 0; //TODO add timestamp to Sig Entry
            hentry.end_timestamp = 0;
            hentry.id = patternentry.original_chunk;

            messies.push_back(hentry);
            list.pop_back();
        } else {
            nicepatterns.push_back(patternentry);
        }
    }
    list = nicepatterns;
}


//return number of elements in total
int PatternUnit::size() const 
{
    if ( cnt == 0 ) {
        return 1; //not repetition
    } else {
        return seq.size()*cnt;
    }
}

string 
PatternUnit::show() const
{
    vector<off_t>::const_iterator iter;
    ostringstream showstr;
    showstr << "( " ;
    for (iter = seq.begin();
            iter != seq.end();
            iter++ )
    {
        showstr << *iter << " ";
    }
    showstr << ") ^" << cnt << endl;
    return showstr.str();
}

string
IdxSigUnit::show() const
{
    ostringstream showstr;
    showstr << init << " ... ";
    showstr << PatternUnit::show();
    return showstr.str();
}



string IdxSigEntry::show()
{
    ostringstream showstr;

    showstr << "[" << original_chunk << "]" 
         << "[" << new_chunk_id << "]" << endl;
    showstr << "----Logical Offset----" << endl;
    showstr << logical_offset.show();
    
    vector<IdxSigUnit>::const_iterator iter2;

    showstr << "----Length----" << endl;
    for (iter2 = length.begin();
            iter2 != length.end();
            iter2++ )
    {
        showstr << iter2->show(); 
    }

    showstr << "----Physical Offset----" << endl;
    for (iter2 = physical_offset.begin();
            iter2 != physical_offset.end();
            iter2++ )
    {
        showstr << iter2->show(); 
    }
    showstr << "-------------------------------------" << endl;

    return showstr.str();
}

// At this time, we only append when:
// For logical off, length and physical off, each of them 
// has only one SigUnit
bool IdxSigEntry::append(IdxSigEntry &other)
{
    IdxSigEntry tmpentry = *this;

    if ( this->length.the_stack.size() == 1
         && this->physical_offset.the_stack.size() == 1
         && other.length.the_stack.size() == 1
         && other.physical_offset.the_stack.size() == 1 
         && tmpentry.logical_offset.append( other.logical_offset )
         && tmpentry.length.the_stack[0].append( other.length.the_stack[0] )
         && tmpentry.physical_offset.the_stack[0].append( 
                                          other.physical_offset.the_stack[0] ) )
    {
        *this = tmpentry;
        return true;
    } else {
        return false;
    }
}




