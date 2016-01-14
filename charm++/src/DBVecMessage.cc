#include "TabaSCo.decl.h"
#include "DBVecMessage.h"

DBVecMessage::DBVecMessage(std::vector<double> vec)
{
    this->length = vec.size();
    this->data = new double[this->length];
    std::memcpy(this->data, vec.data(), sizeof(double) * this->length);
}

DBVecMessage::DBVecMessage(double * vec, int dim)
{
    this->length = dim;
    this->data = new double[dim];
    std::memcpy(this->data, vec, sizeof(double) * dim);
}

void * DBVecMessage::pack(DBVecMessage* m)
{
    int msize = 0;
    msize += sizeof(int);//length
    msize += sizeof(double) * m->length;//data

    char * p = (char*) CkAllocBuffer(m, msize);
    char * buf = p;

    std::memcpy(buf, &m->length, sizeof(int));
    buf += sizeof(int);
    double * dbuf = (double *) buf;
    for(int i = 0; i < m->length; i++)
    {
        dbuf[i] = m->data[i];
    }
    //std::memcpy(buf, &m->data, sizeof(double) * m->length); //What am I stupid about? Why doesn't this equal the above?
    CkFreeMsg(m);
    return (void*) p;
}

DBVecMessage * DBVecMessage::unpack(void *inbuf)
{
    char * buf = (char *) inbuf;
    int len = *(int *) buf;
    buf += sizeof(int);

    DBVecMessage * pmsg =(DBVecMessage*)CkAllocBuffer(inbuf, sizeof(DBVecMessage));

    pmsg = new ((void*)pmsg) DBVecMessage((double *)buf, len);

    CkFreeMsg(inbuf);
    return pmsg;
}


