/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/29 17:38
 */

#ifndef IOE_SORW_CIRCULARQUEUE_HPP
#define IOE_SORW_CIRCULARQUEUE_HPP

#include <iostream>

#define MAX_QUEUE_SIZE 65536
typedef WalkDataType ElemType;

typedef struct QNode {
    ElemType data;
    QNode *next;

} QNode, *QueuePtr;      //节点

//循环队列
typedef struct {
    ElemType *base;
    int front;
    int rear;
} SqQueue;

class CircularQueue {
public:
    void InitQueue();            //初始化队列
    void DestroyQueue();        //销毁队列
    void ClearQueue();            //清空队列
    bool QueueEmpty() const;            //队列是否为空
    bool QueueFull() const;

    int QueueLength() const;            //队列长度
    void Enqueue(ElemType val);    //在队尾插入数据
    void DeQueue(ElemType &val);    //删除队头

private:
    SqQueue q;
};

//初始化队列
void CircularQueue::InitQueue() {
    q.base = (ElemType *) malloc(sizeof(ElemType) * MAX_QUEUE_SIZE);
    if (!q.base) {
        //如果分配失败
        assert(false);
    }
    q.front = q.rear = 0;
}

//销毁队列
void CircularQueue::DestroyQueue() {
    free(q.base);
    q.front = q.rear = 0;
}

//在队尾插入数据
void CircularQueue::Enqueue(ElemType val) {
    if ((q.rear + 1) % MAX_QUEUE_SIZE == q.front) {
        assert(false);
    }
    q.base[q.rear] = val;
    q.rear = (q.rear + 1) % MAX_QUEUE_SIZE;
}

//删除队头，并返回当前队头的值
void CircularQueue::DeQueue(ElemType &val) {
    if (q.front == q.rear) {
        assert(false);
        return;
    }
    val = q.base[q.front];
    q.front = (q.front + 1) % MAX_QUEUE_SIZE;
}

//清空队列
void CircularQueue::ClearQueue() {
    DestroyQueue();
    InitQueue();
}


//队列是否为空
bool CircularQueue::QueueEmpty() const {
    if (q.front == q.rear)
        return true;
    else
        return false;
}

//队列长度
int CircularQueue::QueueLength() const {
    return (q.rear - q.front + MAX_QUEUE_SIZE) % MAX_QUEUE_SIZE;
}

bool CircularQueue::QueueFull() const {
    return (q.rear + 1) % MAX_QUEUE_SIZE == q.front;
}

#endif //IOE_SORW_CIRCULARQUEUE_HPP
