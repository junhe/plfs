#define _FILE_OFFSET_BITS 64 

#include <string.h>
#include <stdlib.h>
#include "plfs_tool_common.h"
#include "plfs.h"
#include "Container.h"
#include "idxanalyzer.h"


void show_usage(char* app_name) {
    fprintf(stderr, "Usage: %s [-version] <filename>\n", app_name);
}

int main (int argc, char **argv) {

    /*
    IdxSignature mysig;
    ifstream idx_file;
    vector<off_t> off_deltas;
    idxfile::EntryList pb_entrylist;
    IdxSigEntryList sig_entrylist;
    IdxSigEntry myentry;


    IdxSigUnit myunit;
    myunit.init = 2;
    myunit.cnt = 3;
    myunit.seq.push_back(4);
    myunit.seq.push_back(5);
    
    myentry.logical_offset = myunit;

    myunit.seq.clear();
    myunit.seq.push_back(0);
    myunit.seq.push_back(0);
    myunit.init = 1;

     
    SigStack<IdxSigUnit> stack;
    stack.push(myunit);
    
    myentry.length = stack;
    myentry.physical_offset = stack;

    vector<IdxSigEntry> mylist;
    mylist.push_back(myentry);

    printIdxEntries(mylist);

    int pos;
    for ( int i = 0 ; i < 50 ; i++ ) {
        cout << i << ":" << myentry.contains( i, pos ) << endl;
    }

    return 0;

    stack.deSerialize( stack.serialize() );

    myentry.original_chunk = 8848;
    myentry.logical_offset = myunit;
    myentry.length.push(myunit);
    myentry.length.push(myunit);
    myentry.physical_offset.push(myunit);
    myentry.physical_offset.push(myunit);
    myentry.physical_offset.push(myunit);

    vector<IdxSigEntry> alist;
    alist.push_back(myentry);
    printIdxEntries(alist);
    alist[0].deSerialize(alist[0].serialize());
    printIdxEntries(alist);







    return -1;
    //////////////////////
    */
    int i;
    char *target;
    bool found_target = false;
    for (i = 1; i < argc; i++) {
        plfs_handle_version_arg(argc, argv[i]);
        if (strcmp(argv[i], "-nc") == 0) {
            // silently ignore deprecated argument
        } else if (!found_target) {
            target = argv[i];
            found_target = true;
        } else {
            // Found more than one target. This is an error.
            show_usage(argv[0]);
            exit(1);
        }
    }
    if (!found_target) {
        show_usage(argv[0]);
        exit(1);
    }

    int ret = plfs_dump_index(stderr,target,0);
    if ( ret != 0 ) {
        fprintf(stderr, "Error: %s is not in a PLFS mountpoint"
                " configured with 'workload n-1'\n", target);
    }
    exit( ret );
}
