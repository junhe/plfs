#ifndef __Index_H__
#define __Index_H__

#include "COPYRIGHT.h"
#include <set>
#include <map>
#include <vector>
#include <list>
using namespace std;

#include "Util.h"
#include "Metadata.h"
#include "idxanalyzer.h"
#include "patternanalyzer.h"


enum IndexEntryType {
    SINGLEHOST = 0,  // HostEntry class
    SIMPLEFORMULA = 1, // SimpleFormulaEntry class
    COMPLEXPATTERN = 2, //IdxSigEntryList class
    MULTILEVEL = 3,
    UNKNOWNTYPE = 7, //IdxSigEntryList class 
};



// the LocalEntry (HostEntry) and the ContainerEntry should maybe be derived
// from one another. there are two types of index files
// on a write, every host has a host index
// on a read, all the host index files get aggregated into one container index

class IndexFileInfo
{
    public:
        IndexFileInfo();
        void *listToStream(vector<IndexFileInfo> &list,int *bytes);
        vector<IndexFileInfo> streamToList(void *addr);
        //bool operator<(IndexFileInfo d1);
        double timestamp;
        string hostname;
        IndexEntryType type;
        pid_t  id;
};

// this is the class that represents one record in the in-memory
// data structure that is
// the index for the entire container (the aggregation of the multiple host
// index files).
// this in-memory structure is used to answer read's by finding the appropriate
// requested logical offset within one of the physical host index files
class ContainerEntry : HostEntry
{
    public:
        bool mergable( const ContainerEntry& );
        bool abut( const ContainerEntry& );
        bool follows( const ContainerEntry& );
        bool preceeds( const ContainerEntry& );
        ContainerEntry split(off_t); //split in half, this is back, return front

    protected:
        pid_t original_chunk;   // we just need to track this so we can
        // rewrite the index appropriately for
        // things like truncate to the middle or
        // for the as-yet-unwritten index flattening

        friend ostream& operator <<(ostream&,const ContainerEntry&);

        friend class Index;
};

// this is a way to associate a integer with a local file
// so that the aggregated index can just have an int to point
// to a local file instead of putting the full string in there
typedef struct {
    string path;
    int fd;
} ChunkFile;

class IdxSigEntryList;
class IdxSignature;
class Index : public Metadata
{
    public:
        Index( string );
        Index( string path, int fd );
        ~Index();

        int readIndex( string hostindex );

        void setPath( string );

        bool ispopulated( );

        void addWrite( off_t offset, size_t bytes, pid_t, double, double );

        size_t memoryFootprintMBs();    // how much area the index is occupying

        void flushCompressedIndexBuf();
        int flushHostIndexBuf();
        int flush();

        off_t lastOffset( );

        void lock( const char *function );
        void unlock(  const char *function );

        int getFd() {
            return fd;
        }
        void resetFd( int fd ) {
            this->fd = fd;
        }

        int resetPhysicalOffsets();

        size_t totalBytes( );

        int getChunkFd( pid_t chunk_id );

        int setChunkFd( pid_t chunk_id, int fd );

        int globalMultiLevelLookup( int *fd, off_t *chunk_off, size_t *length,
                string& path, bool *hole, pid_t *chunk_id,
                off_t logical );
        int globalComplexLookup( int *fd, off_t *chunk_off, size_t *length,
                string& path, bool *hole, pid_t *chunk_id,
                off_t logical );
        int globalLookup( int *fd, off_t *chunk_off, size_t *length,
                string& path, bool *hole, pid_t *chunk_id,
                off_t logical );

        int insertGlobal( ContainerEntry * );
        int insertGlobal( IdxSigEntry *g_entry);
        void merge( Index *other);
        void truncate( off_t offset );
        int rewriteIndex( int fd );
        void truncateHostIndex( off_t offset );

        void compress();
        int debug_from_stream(void *addr);
        int global_to_file(int fd);
        int global_from_stream(void *addr);
        int global_to_stream(void **buffer,size_t *length);
        int global_to_stream( string &buf );
        friend ostream& operator <<(ostream&, Index&); //TODO: make all my show() const
        // Added to get chunk path on write
        string index_path;
        map<int, string> index_paths;
        void startBuffering();
        void stopBuffering();
        bool isBuffering();
        int getHostIndexSize();
        IndexEntryType type;  //TODO: figure out when to set this and when 
                              //      can use it. Probably use getType() setType
                              //      Shall I decide the type at the time creating
                              //      object
    private:
        void init( string );
        int chunkFound( int *, off_t *, size_t *, off_t,
                string&, pid_t *, ContainerEntry * );
        int chunkFound( int *fd, off_t *chunk_off, size_t *chunk_len,
                off_t shift, string& path, pid_t *chunk_id,
                IdxSigEntry *entry, int pos );
        int readComplexIndex( string hostindex );
        int readMultiLevelIndex( string hostindex );
        int cleanupReadIndex(int, void *, off_t, int, const char *,
                const char *);
        void *mapIndex( string, int *, off_t * );
        void resetFD( int fd, string indexpath );
        int handleOverlap( ContainerEntry& g_entry,
                pair< map<off_t,ContainerEntry>::iterator,
                bool > &insert_ret );
        map<off_t,ContainerEntry>::iterator insertGlobalEntryHint(
                ContainerEntry *g_entry ,map<off_t,ContainerEntry>::iterator hint);
        pair<map<off_t,ContainerEntry>::iterator,bool> insertGlobalEntry(
                ContainerEntry *g_entry);
        void insertGlobalEntry( IdxSigEntry *g_entry);
        size_t splitEntry(ContainerEntry *,set<off_t> &,
                multimap<off_t,ContainerEntry> &);
        void findSplits(ContainerEntry&,set<off_t> &);
        // where we buffer the host index (i.e. write)
        vector< HostEntry > hostIndex;
        
        //my own buffer. Complex pattern
        //analysis is on this buffer.
        IdxSigEntryList complexIndexBuf; //this is used to hold the
                                         //recognized complex patterns
        
        map< off_t, int > global_complex_index_map; // map from logical offset to its
                                                    // position in global_complex_index_list
        IdxSigEntryList global_complex_index_list;
        

        // For multi-level patterns
        MultiLevel::PatternCombo multilevelIndexbuf; // to hold recognized
                                                     // multi-level patterns
        MultiLevel::PatternCombo global_multilevel_index;  // use for reading.


        IdxSignature complexIndexUtil;    //a tool class to analyze pattern

        // this is a global index made by aggregating multiple locals
        map< off_t, ContainerEntry > global_index;

        // this is a way to associate a integer with a local file
        // so that the aggregated index can just have an int to point
        // to a local file instead of putting the full string in there
        vector< ChunkFile >       chunk_map;

        // keep index file descriptor here, one for each index entry type
        map< IndexEntryType, int > fds;

        // need to remember the current offset position within each chunk
        map<pid_t,off_t> physical_offsets;

        bool   populated;
        pid_t  mypid;
        string physical_path;
        int    chunk_id;
        off_t  last_offset;
        size_t total_bytes;
        int    fd;
        bool buffering;    // are we buffering the index on write?
        bool buffer_filled; // were we buffering but ran out of space?
        pthread_mutex_t    fd_mux;   // to allow thread safety

        bool compress_contiguous; // set true for performance. 0 for tracing.
        bool enable_hash_lookup;
};

inline
int Index::globalMultiLevelLookup( int *fd, off_t *chunk_off, size_t *chunk_len,
        string& path, bool *hole, pid_t *chunk_id,
        off_t logical )
{
    *hole = false;
    *chunk_id = (pid_t)-1;
}


#define MAP_ITR map<off_t,ContainerEntry>::iterator
#define COMPLEXMAP_ITR map<off_t,int>::iterator
inline
int Index::globalComplexLookup( int *fd, off_t *chunk_off, size_t *chunk_len,
        string& path, bool *hole, pid_t *chunk_id,
        off_t logical )
{
    *hole = false;
    *chunk_id = (pid_t)-1;

    /////////////////
    // For debug
    // Return an offset instantly
    /*
    if ( global_complex_index_list.list.size() == 0 
         || logical >= 1048576) {
        *fd = -1;
        *chunk_len = 0;
        return 0;
    }
    
    IdxSigEntry entry = global_complex_index_list.list[0];
    
    *chunk_off = 0;
    *chunk_len = 1;
    *chunk_id = 0;
    *fd = chunk_map[0].fd;
    path = chunk_map[0].path;
    return 0;
    */

    //////////////////////////////////////////
    /////////// without hashtable  /////////////////////////
    if ( enable_hash_lookup == false ) {
        vector<IdxSigEntry>::iterator itr;

        if ( global_complex_index_list.list.size() == 0 ) {
            *fd = -1;
            *chunk_len = 0;
            return 0;
        }

        for ( itr =  global_complex_index_list.list.begin();
              itr != global_complex_index_list.list.end();
              itr++ )
        {
            int pos = -1;
            if ( itr->contains( logical, pos ) ) {
                //oss << " !!!!!!!!Found at pos:" << pos << endl; 
                //mlog(IDX_WARN, "%s", oss.str().c_str());
                return chunkFound( fd, chunk_off, chunk_len,
                        logical - itr->logical_offset.getValByPos(pos), path,
                        chunk_id, &(*itr), pos );
            } //else {
                //oss << "Not Found!" << endl;
                //mlog(IDX_WARN, "%s", oss.str().c_str());
            //}

        }

        //TODO: handle it...
        mlog(IDX_WARN, "%s in a hole. or off the end of the file", __FUNCTION__);
        
        *fd = -1;
        *chunk_len = 0;
        return 0;
    } else {
        //////////////////////////////////////////
        /////////// with hashtable  /////////////////////////
        COMPLEXMAP_ITR itr;
        COMPLEXMAP_ITR prev = (COMPLEXMAP_ITR)NULL;

        itr = global_complex_index_map.lower_bound(logical);
        
        if ( global_complex_index_list.list.size() == 0 ) {
            *fd = -1;
            *chunk_len = 0;
            return 0;
        }
        
        if ( itr == global_complex_index_map.end() ) {
            itr--;
        }

        if ( itr != global_complex_index_map.begin() ) {
            prev = itr;
            prev--;
        }
      
        // Now we have:
        // itr->first >= logical
        // AND pre-first < logical

        if ( itr->first > logical ) {
            int previous_pos = -1;
            if ( prev != (COMPLEXMAP_ITR)NULL ) {
                IdxSigEntry &previous 
                    = global_complex_index_list.list[prev->second];
                if ( previous.contains( logical, previous_pos ) ) {
                    return chunkFound( fd, chunk_off, chunk_len,
                            logical - previous.logical_offset.getValByPos(previous_pos), 
                            path, chunk_id, &previous, previous_pos );
                }
            }
        } else {
            // Now we have:
            // itr->first <= logical
            IdxSigEntry &entry = global_complex_index_list.list[itr->second];

            int entry_pos = -1;
            if ( entry.contains( logical, entry_pos ) ) {
                return chunkFound( fd, chunk_off, chunk_len,
                        logical - entry.logical_offset.getValByPos(entry_pos), path,
                        chunk_id, &entry, entry_pos );
            }
        }

        //TODO: handle it...
        mlog(IDX_WARN, "%s in a hole. or off the end of the file", __FUNCTION__);
        
        *fd = -1;
        *chunk_len = 0;
        return 0;
    }
}



#endif
