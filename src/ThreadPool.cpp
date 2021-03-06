#include <pthread.h>
#include <errno.h>
#include <string.h>
#include "ThreadPool.h"
#include "Util.h"
#include "mlogfacs.h"

ThreadPool::ThreadPool( size_t size, void *(*func) (void *), void *args )
{
    pthread_t *threads = new pthread_t[size];
    thread_error = 0;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    mlog(INT_DAPI, "THREAD_POOL: Creating %lu threads", (unsigned long)size );
    for( size_t t = 0; t < size; t++ ) {
        if ( 0 != pthread_create(&threads[t], &attr, func, args) ) {
            thread_error = errno;
            mlog(INT_DRARE, "THREAD_POOL: create error %s", strerror(errno) );
            break;
        }
    }
    if ( thread_error == 0 ) {
        for(size_t t=0; t < size; t++) {
            void *status;
            if ( 0 != pthread_join(threads[t], &status) ) {
                thread_error = errno;
                mlog(INT_DRARE, "THREAD_POOL: join error %s", strerror(errno));
                break;
            } else {
                stati.push_back(status);
            }
        }
    }
    // clean up the thread stuff
    pthread_attr_destroy(&attr);
    delete []threads;
}

ThreadPool::~ThreadPool()
{
}

int ThreadPool::threadError()
{
    return thread_error;
}

vector<void *> * ThreadPool::getStati()
{
    return &stati;
}
