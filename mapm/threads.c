#ifdef THREAD

#include <sys/sysinfo.h>

#include <stdlib.h>
#include <pthread.h>

#include "mapm.h"

int threads_ncpu = 0;
int threads_nthreads = 0;


typedef struct tpool_work {
    void (*routine)(void *, int);

    void *arg;
    struct tpool_work *next;
} tpool_work_t;


typedef struct tpool {
    /* pool characteristics */
    int active_threads;
    int allowed_threads;
    int num_threads;
    int max_queue_size;

    int do_not_block_when_full;
    /* pool state */
    pthread_t *threads;
    int index;
    pthread_mutex_t index_lock;
    int cur_queue_size;
    tpool_work_t *queue_head;
    tpool_work_t *queue_tail;
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_not_empty;
    pthread_cond_t queue_not_full;
    pthread_cond_t queue_empty;

    int cur_results_size;
    tpool_work_t *results_head;
    tpool_work_t *results_tail;
    pthread_mutex_t results_lock;
    pthread_cond_t results_not_empty;

    int queue_closed;
    int shutdown;
} *tpool_t;


typedef struct {
    MAP *map;
    int id;
    int converged;
} converge_request_t;


static time_t threads_last_throttle = 0;
static tpool_t converge_pool;


static void tpool_init(tpool_t *tpoolp, int num_worker_threads, int max_queue_size, int do_not_block_when_full);

static void *tpool_thread(void *arg);

static int tpool_add_work(tpool_t tpool, void (*routine)(void *, int), void *arg);

static int tpool_get_result(tpool_t tpool, void **argp);

static void tpool_destroy(tpool_t tpool, int finish);

static void thread_converge(void *arg, int);

static int get_ncpu(void);

static void tpool_throttle(void);


int get_ncpu(void) {
    int ncpu = get_nprocs();
    return ncpu;
}


void threads_init(void) {
    threads_ncpu = get_ncpu();
    threads_nthreads = threads_ncpu;

    tpool_init(&converge_pool, threads_nthreads, threads_nthreads * 5, 0);
}


int add_converge_request(const MAP *map, int id) {
    converge_request_t *p = malloc(sizeof(converge_request_t));
    p->map = allocate_map(map->num_loci);
    mapcpy(p->map, map, FALSE);
    p->id = id;

    tpool_throttle();

    return tpool_add_work(converge_pool, thread_converge, p);
}


int get_converge_result(MAP *map, int *id) {
    converge_request_t *p = NULL;
    int rv;

    tpool_get_result(converge_pool, (void **) &p);

    mapcpy(map, p->map, FALSE);
    *id = p->id;
    rv = p->converged;
    free_map(p->map);
    free(p);

    return rv;
}


int queue_size(void) {
    return converge_pool->max_queue_size;
}


void thread_converge(void *arg, int index) {
    converge_request_t *p = (converge_request_t *) arg;
    init_rec_fracs(p->map);
    p->converged = converge_to_map_mt(p->map, index);
}


void tpool_init(tpool_t *tpoolp, int num_worker_threads, int max_queue_size, int do_not_block_when_full) {
    int i, rtn;
    tpool_t tpool;

    /* allocate a pool data structure */
    if ((tpool = malloc(sizeof(struct tpool))) == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    /* initialize the fields */
    tpool->active_threads = 0;
    tpool->allowed_threads = num_worker_threads;
    tpool->num_threads = num_worker_threads;
    tpool->max_queue_size = max_queue_size;
    tpool->do_not_block_when_full = do_not_block_when_full;
    if ((tpool->threads = malloc(sizeof(pthread_t) * num_worker_threads)) == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    tpool->cur_queue_size = 0;
    tpool->queue_head = NULL;
    tpool->queue_tail = NULL;
    tpool->queue_closed = 0;
    tpool->cur_results_size = 0;
    tpool->results_head = NULL;
    tpool->results_tail = NULL;
    tpool->shutdown = 0;
    tpool->index = 1;
    if ((rtn = pthread_mutex_init(&tpool->index_lock, NULL)) != 0) {
        perror("pthread_mutex_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_mutex_init(&tpool->queue_lock, NULL)) != 0) {
        perror("pthread_mutex_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_cond_init(&tpool->queue_not_empty, NULL)) != 0) {
        perror("pthread_cond_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_cond_init(&tpool->queue_not_full, NULL)) != 0) {
        perror("pthread_cond_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_cond_init(&tpool->queue_empty, NULL)) != 0) {
        perror("pthread_cond_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_mutex_init(&tpool->results_lock, NULL)) != 0) {
        perror("pthread_mutex_init");
        exit(EXIT_FAILURE);
    }
    if ((rtn = pthread_cond_init(&tpool->results_not_empty, NULL)) != 0) {
        perror("pthread_cond_init");
        exit(EXIT_FAILURE);
    }

    /* create threads */
    for (i = 0; i < num_worker_threads; i++) {
        if ((rtn = pthread_create(&tpool->threads[i], NULL, tpool_thread, (void *) tpool)) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    *tpoolp = tpool;
}


void threads_shutdown(void) {
    tpool_destroy(converge_pool, 1);
}


void *tpool_thread(void *arg) {
    tpool_t tpool = (tpool_t) arg;
    tpool_work_t *my_workp;
    int index;
    int decrement = 0;

    pthread_mutex_lock(&tpool->index_lock);
    index = tpool->index++;
    pthread_mutex_unlock(&tpool->index_lock);

    while (1) {
        pthread_mutex_lock(&tpool->queue_lock);
        if (decrement) tpool->active_threads--;
        while ((tpool->cur_queue_size == 0 ||
                tpool->active_threads >= tpool->allowed_threads) &&
               !tpool->shutdown) {
            pthread_cond_wait(&tpool->queue_not_empty, &tpool->queue_lock);
        }

        if (tpool->shutdown) {
            pthread_mutex_unlock(&tpool->queue_lock);
            pthread_exit(NULL);
        }

        tpool->active_threads++;

        my_workp = tpool->queue_head;
        tpool->cur_queue_size--;
        if (tpool->cur_queue_size == 0) {
            tpool->queue_head = tpool->queue_tail = NULL;
        } else {
            tpool->queue_head = my_workp->next;
        }

        if (!tpool->do_not_block_when_full &&
            tpool->cur_queue_size == tpool->max_queue_size - 1) {
            pthread_cond_signal(&tpool->queue_not_full);
        }

        if (tpool->cur_queue_size == 0) {
            pthread_cond_signal(&tpool->queue_empty);
        }

        pthread_mutex_unlock(&tpool->queue_lock);

        (*(my_workp->routine))(my_workp->arg, index);

        //free(my_workp);

        pthread_mutex_lock(&tpool->results_lock);
        my_workp->next = NULL;
        if (tpool->cur_results_size == 0) {
            tpool->results_tail = tpool->results_head = my_workp;
            pthread_cond_signal(&tpool->results_not_empty);
        } else {
            tpool->results_tail->next = my_workp;
            tpool->results_tail = my_workp;
        }
        tpool->cur_results_size++;
        pthread_mutex_unlock(&tpool->results_lock);
        decrement = 1;
    }
}


int tpool_add_work(tpool_t tpool, void (*routine)(void *, int), void *arg) {
    tpool_work_t *workp;

    pthread_mutex_lock(&tpool->queue_lock);

    if (tpool->cur_queue_size == tpool->max_queue_size &&
        tpool->do_not_block_when_full) {
        pthread_mutex_unlock(&tpool->queue_lock);
        return 0;
    }

    while (tpool->cur_queue_size == tpool->max_queue_size &&
           !(tpool->shutdown || tpool->queue_closed)) {
        pthread_cond_wait(&tpool->queue_not_full, &tpool->queue_lock);
    }

    if (tpool->shutdown || tpool->queue_closed) {
        pthread_mutex_unlock(&tpool->queue_lock);
        return -1;
    }

    /* allocate work structure */
    workp = malloc(sizeof(tpool_work_t));
    workp->routine = routine;
    workp->arg = arg;
    workp->next = NULL;
    if (tpool->cur_queue_size == 0) {
        tpool->queue_tail = tpool->queue_head = workp;
        pthread_cond_signal(&tpool->queue_not_empty);
    } else {
        tpool->queue_tail->next = workp;
        tpool->queue_tail = workp;
        pthread_cond_signal(&tpool->queue_not_empty);
    }
    tpool->cur_queue_size++;

    pthread_mutex_unlock(&tpool->queue_lock);

    return 1;
}


int tpool_get_result(tpool_t tpool, void **argp) {
    tpool_work_t *workp;

    pthread_mutex_lock(&tpool->results_lock);

    while (tpool->cur_results_size == 0 && !tpool->shutdown) {
        pthread_cond_wait(&tpool->results_not_empty, &tpool->results_lock);
    }

    if (tpool->shutdown) {
        pthread_mutex_unlock(&tpool->results_lock);
        return -1;
    }

    workp = tpool->results_head;
    tpool->cur_results_size--;
    if (tpool->cur_results_size == 0) {
        tpool->results_head = tpool->results_tail = NULL;
    } else {
        tpool->results_head = workp->next;
    }

    pthread_mutex_unlock(&tpool->results_lock);

    *argp = workp->arg;
    free(workp);

    return 1;
}


void tpool_destroy(tpool_t tpool, int finish) {
    int i, rtn;
    tpool_work_t *cur_nodep;

    if ((rtn = pthread_mutex_lock(&tpool->queue_lock)) != 0) {
        perror("pthread_mutex_lock");
        exit(EXIT_FAILURE);
    }

    /* Is a shutdown already in progress? */
    if (tpool->queue_closed || tpool->shutdown) {
        if ((rtn = pthread_mutex_unlock(&tpool->queue_lock)) != 0) {
            perror("pthread_mutex_unlock");
            exit(EXIT_FAILURE);
        }
        return;
    }

    tpool->queue_closed = 1;

    /* If finish flag is set, wait for workers to drain queue */
    if (finish) {
        while (tpool->cur_queue_size != 0) {
            if ((rtn = pthread_cond_wait(&tpool->queue_empty,
                                         &tpool->queue_lock)) != 0) {
                perror("pthread_cond_wait");
                exit(EXIT_FAILURE);
            }
        }
    }

    tpool->shutdown = 1;

    if ((rtn = pthread_mutex_unlock(&tpool->queue_lock)) != 0) {
        perror("pthread_mutex_unlock");
        exit(EXIT_FAILURE);
    }

    /* Wake up any workers so they recheck shutdown flag */
    if ((rtn = pthread_cond_broadcast(&tpool->queue_not_empty)) != 0) {
        perror("pthread_cond_broadcast");
        exit(EXIT_FAILURE);
    }

    if ((rtn = pthread_cond_broadcast(&tpool->queue_not_full)) != 0) {
        perror("pthread_cond_broadcast");
        exit(EXIT_FAILURE);
    }

    /* Wait for workers to exit */
    for (i = 0; i < tpool->num_threads; i++) {
        if ((rtn = pthread_join(tpool->threads[i], NULL)) != 0) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }

    /* Now free pool structure */
    free(tpool->threads);
    while (tpool->queue_head != NULL) {
        cur_nodep = tpool->queue_head->next;
        tpool->queue_head = tpool->queue_head->next;
        free(cur_nodep);
    }
    while (tpool->results_head != NULL) {
        cur_nodep = tpool->results_head->next;
        tpool->results_head = tpool->results_head->next;
        free(cur_nodep);
    }
    free(tpool);
}


void tpool_throttle(void) {
    time_t cur = time(NULL);

    if (labs(threads_last_throttle - cur) > 60) {
        int nthr;
        double avg[3];
        threads_last_throttle = cur;

        if (getloadavg(avg, 1) == -1) {
            perror("getloadavg");
        }

        pthread_mutex_lock(&converge_pool->queue_lock);

        nthr = (int) (threads_ncpu - (avg[0] - converge_pool->allowed_threads) + .5);
        if (nthr < 1) nthr = 1;
        if (nthr > threads_ncpu) nthr = threads_ncpu;
        printf("1min load avg: %f\tNCPU:%d\tNthreads:%d\tSetting num. of allowed threads to %d\n",
               avg[0], threads_ncpu, threads_nthreads, nthr);
        converge_pool->allowed_threads = nthr;

        pthread_mutex_unlock(&converge_pool->queue_lock);
    }
}

#endif
