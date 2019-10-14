#ifndef _THREADS_H_
#define _THREADS_H_

extern int threads_ncpu;
extern int threads_nthreads;

void threads_init(void);

void threads_shutdown(void);

int add_converge_request(const MAP *map, int id);

int get_converge_result(MAP *map, int *id);

int queue_size(void);

#endif
