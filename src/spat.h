#ifndef __SPAT_H_
#define __SPAT_H_

typedef struct{
    int seq;
    int pos;
    int next;
    int prev;
}sOffset_t;

typedef struct{
    char *string;
    int length;
    int support;
    sOffset_t *offset;
}sPat_t;



#endif  /*__SPAT_H_*/
