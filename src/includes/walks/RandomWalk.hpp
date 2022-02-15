//
// Created by Ethan on 2021/1/9.
//

#ifndef IOE_SORW_RANDOMWALK_HPP
#define IOE_SORW_RANDOMWALK_HPP
class RandomWalk{
public:

    virtual void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m);

    virtual void
    updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
               WalkManager &walkManager);

    virtual bool allWalksFinished(WalkManager &walkManager);

    virtual bid_t getBlock(vid_t v);

    virtual void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex);

};

void RandomWalk::updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
                            WalkManager &walkManager) {

}

bool RandomWalk::allWalksFinished(WalkManager &walkManager) {
    return walkManager.currentNWalks == 0;
}

void RandomWalk::setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) {
}

void RandomWalk::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {

}

bid_t RandomWalk::getBlock(vid_t v) {
    return 0;
}

#endif //IOE_SORW_RANDOMWALK_HPP
