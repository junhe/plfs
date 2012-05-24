#include "idxanalyzer.h"
#include "Util.h"
#include "plfs_private.h"

#include <algorithm>

//for debugging
void printIdxEntries( vector<IdxSigEntry> &idx_entry_list )
{
    vector<IdxSigEntry>::const_iterator iter;

    //cout << "this is printIdxEntries" << endl;

    for ( iter = idx_entry_list.begin();
            iter != idx_entry_list.end();
            iter++ )
    {
        cout << "[" << iter->original_chunk << "]" 
             << "[" << iter->new_chunk_id << "]" << endl;
        cout << "----Logical Offset----" << endl;
        iter->logical_offset.show();
        
        vector<IdxSigUnit>::const_iterator iter2;

        cout << "----Length----" << endl;
        for (iter2 = iter->length.begin();
                iter2 != iter->length.end();
                iter2++ )
        {
            iter2->show(); 
        }

        cout << "----Physical Offset----" << endl;
        for (iter2 = iter->physical_offset.begin();
                iter2 != iter->physical_offset.end();
                iter2++ )
        {
            iter2->show(); 
        }
        cout << "-------------------------------------" << endl;
    }
}

vector<off_t> buildDeltas( vector<off_t> seq ) 
{
    vector<off_t>::iterator it;
    vector<off_t> deltas;
    for ( it = seq.begin() ; it != seq.end() ; it++ )
    {
        if ( it > seq.begin() ) {
            deltas.push_back( *it - *(it-1) );
        }
    }
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
    
    logical_offset_delta = buildDeltas(logical_offset);
    length_delta = buildDeltas(length);
    physical_offset_delta = buildDeltas(physical_offset);

    SigStack<IdxSigUnit> offset_sig = 
        discoverSigPattern(logical_offset_delta, logical_offset);
    //offset_sig.show();
    //cout << "Just showed offset_sig" << endl;
    //Now, go through offset_sig one by one and build the IdxSigEntry s
    vector<IdxSigEntry>idx_entry_list;
    vector<IdxSigUnit>::const_iterator stack_iter;

    int range_start = 0, range_end; //the range currently processing
    for (stack_iter = offset_sig.begin();
            stack_iter != offset_sig.end();
            stack_iter++ )
    {
        //cout << stack_iter->init << " " ;
        IdxSigEntry idx_entry;
        range_end = range_start + stack_iter->size();

        vector<off_t>::iterator lstart, lend;
        lstart = length_delta.begin() + range_start;
        if ( length_delta.end() - (length_delta.begin() + range_end) > 0 ) {
            lend = length_delta.begin() + range_end;
        } else {
            lend = length_delta.end();
        }
        SigStack<IdxSigUnit> length_stack = 
            discoverSigPattern( 
                    vector<off_t> (lstart, lend),
                    vector<off_t> (length.begin()+range_start,
                        length.begin()+range_end) ); //this one pointed by length.begin()+range_end 
                                                     //won't be paseed in. The total size passed is
                                                     //stack_iter.size();
        //cout << "************End length" << endl;

        SigStack<IdxSigUnit> physical_offset_stack = 
            discoverSigPattern( 
                    vector<off_t>(lstart, lend),
                    vector<off_t> (physical_offset.begin()+range_start,
                        physical_offset.begin()+range_end) );

        idx_entry.original_chunk = proc;
        idx_entry.logical_offset = *stack_iter;
        idx_entry.length = length_stack;
        idx_entry.physical_offset = physical_offset_stack;
        
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
        fprintf(stderr, "discoverSigPattern() needs to be used with "
                "seq.size==orig.size or seq.size+1==orig.size");
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

                if ( pu.seq.size() == cur_tuple.length ) {
                    //the subseq in lookahead window repeats
                    //the top pattern in stack.
                    //initial remains the same.
                    pu.cnt++;
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

    return pattern_stack;
}

Tuple IdxSignature::searchNeighbor( vector<off_t> const &seq,
        vector<off_t>::const_iterator p_lookahead_win ) 
{
    vector<off_t>::const_iterator i;     
    int j;
    //cout << "------------------- I am in searchNeighbor() " << endl;

    //i goes left util the begin or reaching window size
    int distance = 0;
    i = p_lookahead_win;
    while ( i != seq.begin() && distance < win_size ) {
        i--;
        distance++;
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
        list.push_back(*iter);
    }

}

void IdxSigEntryList::append( IdxSigEntryList other ) 
{
    append(other.list);
}

void IdxSigEntryList::show()
{
    printIdxEntries(list);
}



void idxSigUnit2PBSigUnit( const IdxSigUnit &iunit, idxfile::SigUnit *pbunit )
{
    pbunit->set_init(iunit.init);
    vector<off_t>::const_iterator iter;
    for ( iter = iunit.seq.begin();
            iter != iunit.seq.end();
            iter++ )
    {
        pbunit->add_deltas(*iter);
    }
    pbunit->set_cnt(iunit.cnt);
}

void IdxSigEntryList::siglistToPblist(vector<IdxSigEntry> &slist,
        idxfile::EntryList &pblist)
{
    //read out every entry in slist and put it to pblist
    vector<IdxSigEntry>::iterator iter;

    for ( iter = slist.begin();
            iter != slist.end();
            iter++)
    {
        idxfile::Entry *fentry = pblist.add_entry();
        fentry->set_proc( (*iter).original_chunk );
        idxfile::SigUnit *su = fentry->mutable_logical_offset();
        idxSigUnit2PBSigUnit( (*iter).logical_offset, su );

        //length
        vector<IdxSigUnit>::const_iterator iter2;
        for ( iter2 = (*iter).length.begin();
                iter2 != (*iter).length.end();
                iter2++ )
        {
            idxfile::SigUnit *su = fentry->add_length();
            idxSigUnit2PBSigUnit( (*iter2), su);
        }

        //physical offset
        for ( iter2 = (*iter).length.begin();
                iter2 != (*iter).length.end();
                iter2++ )
        {
            idxfile::SigUnit *su = fentry->add_physical_offset();
            idxSigUnit2PBSigUnit( (*iter2), su);
        }
    }
}

void IdxSigEntryList::siglistToPblist()
{
    siglistToPblist(list, pb_list);
}


void IdxSigEntryList::saveToFile(const char *filename)
{
    siglistToPblist();
    fstream output(filename, ios::out | ios::trunc | ios::binary);
    if ( !pb_list.SerializeToOstream(&output) ) {
        cerr<<"failed to write to myfile."<<endl;
    } else {
        cout<<"Write to myfile: OK"<<endl;
    }
    output.close();
}

void IdxSigEntryList::saveToFile(const int fd)
{
    string buf = serialize();
    if ( buf.size() > 0 ) {
        Util::Writen(fd, &buf[0], buf.size());
    }
}

void IdxSigEntryList::readFromFile(char *filename)
{
    fstream input(filename, ios::in | ios::binary);
    if ( !input ) {
        cerr << "can not open my file.\n";
    } else if ( !pb_list.ParseFromIstream(&input) ) {
        cerr << "failed to parse from myfile\n";
    }
    input.close();
}

void IdxSigEntryList::clear()
{
    list.clear();
    pb_list.Clear();
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
    string buf;
    vector<IdxSigEntry>::iterator iter;
    
    bodysize = bodySize();

    appendToBuffer(buf, &bodysize, sizeof(bodysize));
    
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
    int cur_start = 0;

    list.clear();
    
    readFromBuf(buf, &bodysize, cur_start, sizeof(bodysize));
   
    bufsize = buf.size();
    assert(bufsize == bodysize + sizeof(bodysize));
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

//return number of elements in total
int PatternUnit::size() const 
{
    if ( cnt == 0 ) {
        return 1; //not repetition
    } else {
        return seq.size()*cnt;
    }
}

void PatternUnit::show() const
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

void IdxSigUnit::show() const
{
    cout << init << " ... ";
    PatternUnit::show();
}

off_t sumVector( vector<off_t> seq )
{
    vector<off_t>::const_iterator iiter;
    
    off_t sum = 0;
    for ( iiter = iter->seq.begin() ;
          iiter != iter->seq.end() ;
          iiter++ )
    {
        sum += *iiter;
    }
    
    return sum;
}

off_t IdxSigEntry::getLengthByPos( int pos ) const
{
    int cur_pos = 0; //it should always point to sigunit.init

    vector<IdxSigUnit>::const_iterator iter;
    vector<off_t>::const_iterator iiter;


    for ( iter = length.begin() ;
          iter != length.end() ;
          iter++ )
    {
        int range_start, range_end;
        int numoflen = iter->seq.size()*iter->cnt;
        
        if ( cur_pos <= pos && pos < cur_pos + numoflen ) {
            //in the range that pointed to by iter
            int rpos = pos - cur_pos;
            off_t delta_sum = sumVector(iter->seq);
            int remain = rpos % iter->seq.size();
            int factor = rpos / iter->seq.size();

        } else {
            //not in the range of iter
            cur_pos += numoflen; //keep track of current pos
        }

    }
}

bool IdxSigEntry::contain( off_t offset ) const
{
    vector<IdxSigUnit>::const_iterator iter;
    vector<off_t>::const_iterator iiter;
    vector<off_t>::const_iterator off_delta_iter;
    off_t &logical = offset;

    off_t off_init = logical_offset.init;
    off_t delta_sum = 0;
    for ( off_delta_iter = logical_offset.seq.begin();
          off_delta_iter != logical_offset.seq.end();
          off_delta_iter++ )
    {
        delta_sum += *off_delta_iter;
    }
    ostringstream oss;
    oss << delta_sum;
    mlog(IDX_WARN, "delta_sum:%s", oss.str().c_str());
    
    off_t logical_remain = (logical - off_init) % delta_sum;
    assert( logical_remain >= 0 );
    int factor = (logical - off_init - logical_remain) / delta_sum;

    if ( factor >= logical_offset.cnt ) {
        // logical is very large.
        // check the last offset
        //
        // note that there are totally cnt*seq.size() offsets
        // in this class
        // they are:
        // [init][init+d0][init+d0+d1]...[init+(d0+d1+..+dp)*cnt-dp]
        // [init+(d0+d1+..+dp)*cnt] is the 'last+1' offset
    } else {
        off_t sum = 0;
        for ( off_delta_iter = logical_offset.seq.begin();
              off_delta_iter != logical_offset.seq.end()
              && sum < logical_remain;
              off_delta_iter++ )
        {
            sum += *off_delta_iter;
        }
        assert(off_delta_iter != logical_offset.seq.end());

        int pos = off_delta_iter - logical_offset.seq.begin();
        int pre = pos - 1; //should check the one before pos
        // should check logical_offset.seq.size()*factor+pre+1: this
        // is the position of pre if offsets are expanded.
        // Assume the rank for the offset is like:
        // init:0, init+d0:1
    }

    for ( iter = length.begin() ;
          iter != length.end() ;
          iter++ )
    {
        for ( iiter = iter->seq.begin() ;
              iiter != iter->seq.end() ;
              iiter++ )
        {
            // iiter iterates through all deltas of length 
        }
    }
    return false;
}


