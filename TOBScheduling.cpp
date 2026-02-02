#include "TOBscheduling.hpp"
#include "strategy.hpp"

bool scheduling::verifyCellsCritical()
{
    mobility.resize(db->maxMobility + 1);
    for (auto i : hCells)
    {
        if (db->TopoSorts[i.ID] != db->TopoSortsReverse[i.ID])
        { // The mobility of cell is smaller, the priority is higher.
            mobility[db->TopoSortsReverse[i.ID] - db->TopoSorts[i.ID]].insert(i.ID);
            cellMob[i.ID] = db->TopoSortsReverse[i.ID] - db->TopoSorts[i.ID]; // must be >=1
        }
    }
    return true;
}

bool scheduling::init()
{
    bool ret = true;

    // db = strategyInst->getDBInst();
    segNum = db->segmentNum;

    topoPartNum.resize(db->topoOrderMAX_ + 1);
    topoAve.resize(db->topoOrderMAX_ + 1);
    topoTotalNum.resize(db->topoOrderMAX_ + 1);
    topoVar.resize(db->topoOrderMAX_ + 1);
    topoOrder.resize(hCells.size(), -1);
    cellMob.resize(hCells.size(), -1); // the mobility of critical nodes equal to -1
    // lock.resize(hCells.size(), false);
    calLock.resize(hCells.size(), false);
    trash.resize(hCells.size(), false);

    for (auto &t : topoPartNum)
    {
        t.resize(segNum, 0);
    }
    for (auto c : hCells)
    {
        ++topoPartNum[db->TopoSorts[c.ID]][c.segment];
        ++topoTotalNum[db->TopoSorts[c.ID]];
        topoOrder[c.ID] = db->TopoSorts[c.ID];
    }
    for (int t = 0; t < topoTotalNum.size(); ++t)
        topoAve[t] = (double)topoTotalNum[t] / (double)segNum;
    double SD = 0;
    for (int i = 0; i < topoPartNum.size(); ++i)
    {
        for (int j = 0; j < topoPartNum[i].size(); ++j)
        {
            topoVar[i] += (((double)topoPartNum[i][j] - topoAve[i]) * ((double)topoPartNum[i][j] - topoAve[i])) / (double)segNum;
        }
        SD += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    }
    LASTSD2 = SD;
    verifyCellsCritical();
    int cutSize = 0;
    for (auto &net : hNets)
    {
        if (!net.vaild)
            continue;
        std::vector<int> segflag;
        segflag.resize(segNum, -1);
        for (auto s = net.first; s != net.end; ++s)
        {
            int k = hCells[netPins[s]].segment;
            segflag[k] = 1;
        }
        int n = 0;
        for (auto &c : segflag)
            if (c != -1)
                n++;
        n = (n > 0) ? n : 1;
        cutSize += (n - 1) * net.wt;
    }
    printf("-    Before Scheduling Variance Metric =  %.8f \n", SD);
    printf("-    Before Scheduling cut =  %d \n", cutSize);
    
    return ret;
}
double scheduling::calOriginTopoGain(int asapTopo, int segment)
{
    double tmpGain = 0;
    double oldPartNum = (double)topoPartNum[asapTopo][segment];
    double newPartNum = oldPartNum - 1;
    double oldNum = (double)topoTotalNum[asapTopo];
    double newNum = oldNum - 1;
    double oldAve = topoAve[asapTopo];
    double newAve = oldAve - (double)(1.0 / (double)segNum);

    tmpGain = (((double)(1.0) / (double)segNum) * ((double)(1.0) / (double)segNum) +
               (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
               (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum));
    tmpGain = tmpGain / newNum;
    tmpGain = tmpGain / newNum;
    // tmpGain /= (newNum * newNum);
    double fixTerm1 = (double)(topoVar[asapTopo] / oldNum / oldNum) -
                      (double)(topoVar[asapTopo] / newNum / newNum);

    fixTerm1 = fixTerm1 - tmpGain;

    /*for debug*/
    // {
    //     double ss = 0;
    //     for (int i = 0; i <= db->topoOrderMAX_; ++i)
    //     {
    //         ss += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    //     }

    //     double t = ((double)(1.0) / (double)(segNum)) * ((double)(1.0) / (double)(segNum)) +
    //                (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
    //                (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum);
    //     // topoVar[asapTopo] += t;
    //     // topoTotalNum[asapTopo]--; // maintain the origin topo order info.
    //     // topoAve[asapTopo] -= (double)(1.0 / (double)segNum);
    //     // topoPartNum[asapTopo][segment]--;
    //     SD2 = 0;
    //     for (int i = 0; i <= db->topoOrderMAX_; ++i)
    //     {
    //         if (i == asapTopo)
    //         {
    //             double a = topoVar[asapTopo] + t;
    //             double b = (double)topoTotalNum[i] - 1;
    //             SD2 += (double)(a) / (double)(newNum) / (double)(newNum);
    //         }
    //         else
    //             SD2 += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    //     }

    //     if (fixTerm1 != ss - SD2)
    //         printf(" cal_Actual Variance :%.16f, Variance Metric =  %.16f \n\n", fixTerm1, ss - SD2);
    //     int dd = 0;
    // }
    return fixTerm1;
}
double scheduling::calOtherTopoGain(int topo, int segment)
{
    double tmpGain = 0;
    double oldPartNum = (double)topoPartNum[topo][segment];
    double newPartNum = oldPartNum + 1;
    double oldNum = (double)topoTotalNum[topo];
    double newNum = oldNum + 1;
    double oldAve = topoAve[topo];
    double newAve = topoAve[topo] + (double)(1.0 / (double)segNum);

    tmpGain = ((double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
               (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
               (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum));
    tmpGain /= (newNum * newNum);
    double fixTerm2 = (double)(topoVar[topo] / oldNum / oldNum) -
                      (double)(topoVar[topo] / newNum / newNum);
    fixTerm2 = fixTerm2 - tmpGain;
    /* for debug*/
    // {
    //     double ss = 0;
    //     for (int i = 0; i <= db->topoOrderMAX_; ++i)
    //     {
    //         ss += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    //     }
    //     double t = (double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
    //                (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
    //                (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum);
    //     // topoVar[asapTopo] += t;
    //     // topoTotalNum[asapTopo]--; // maintain the origin topo order info.
    //     // topoAve[asapTopo] -= (double)(1.0 / (double)segNum);
    //     // topoPartNum[asapTopo][segment]--;
    //     SD2 = 0;
    //     for (int i = 0; i < topoPartNum.size(); ++i)
    //     {
    //         if (i == topo)
    //         {
    //             double a = (double)topoVar[i] + t;
    //             double b = (double)topoTotalNum[i] + 1;
    //             SD2 += a / b / b;
    //         }
    //         else
    //             SD2 += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    //     }
    //     if (fixTerm2 != ss - SD2)
    //         printf(" cal_Actual Variance :%.16f, Variance Metric =  %.16f \n\n", fixTerm2, ss - SD2);
    //     int dd = 0;
    // }
    return fixTerm2;
}
int scheduling::obtainGain(uint32_t c)
{
    std::vector<Gain<double>> best;

    int asapTopo = topoOrder[c];
    int segment = hCells[c].segment;
    int bestTopo = -1;
    double bestGain = std::numeric_limits<double>::min();
    std::vector<double> gain;
    gain.resize(db->TopoSortsReverse[c] + 1, 0);

    /* base gain: the decrease var of current topological order db->TopoSorts[c]*/
    double tmpGain1 = calOriginTopoGain(asapTopo, segment);

    /* extention gain: the decrease var of potential topological order db->TopoSorts[c + i]*/
    // for (int i = topoOrder[c] + 1; i <= db->TopoSortsReverse[c]; i++)
    int i = topoOrder[c] + 1;
    for (int t = 0; t < cellMob[c]; ++t)
    {
        double tmpGain2 = calOtherTopoGain(i, segment);
        gain[i] = tmpGain1 + tmpGain2;
        if (gain[i] >= bestGain)
        {
            bestGain = gain[i];
            bestTopo = i;
        }
        i++;
    }
    best.emplace_back(Gain<double>(c, asapTopo, bestTopo, bestGain));
    // best.cellID = c;
    // best.cellGain = bestGain;
    // best.targetSegment = bestTopo;
    // best.sourceSegment = asapTopo;
    if (bestTopo > 0)
        trace.emplace_back(std::pair<double, std::vector<Gain<double>>>(bestGain, best));

    return bestTopo;
}

// bool scheduling::updateNeighbor(uint32_t c)
// {
//     for (auto &netID : db->incidentNets[c])
//     {
//         auto &net = hNets[netID];
//         uint32_t sourceCell = netPins[net.first];
//         if (sourceCell == c) // c is the source node
//         {
//             for (auto s = net.first; s != net.end; ++s)
//             {
//                 int cid = netPins[s];
//                 if (cid == c)
//                     continue;
//                 auto &cell = hCells[cid];
//                 if (topoOrder[cell.ID] == topoOrder[c])
//                 {
//                     if (cellMob[cell.ID] > 0)
//                     {
//                         auto &it = mobility[cellMob[cell.ID]];
//                         it.erase(cell.ID);
//                         cellMob[cell.ID]--;
//                         mobility[cellMob[cell.ID]].insert(cell.ID);
//                         topoOrder[cell.ID] = topoOrder[c];
//                         updateNeighbor(cell.ID);
//                     }
//                 }
//             }
//         }
//         else // c is the drain node
//         {
//             cellMob[sourceCell]--;
//         }
//     }
//     return true;
// }

bool scheduling::cancelOnePropagation()
{
    // int endIndex = trace.size() - 1;
    // for (auto t : trace[endIndex].second)
    // {
    //     uint32_t c = t.cellID;
    //     if (cellMob[c] >= 0 && mobility[cellMob[c]].find(c) != mobility[cellMob[c]].end())
    //         mobility[cellMob[c]].erase(c);
    //     int beforTopo = t.sourceSegment;
    //     int afterTopo = t.targetSegment;
    //     if (db->TopoSortsReverse[c] - t.sourceSegment)
    //         cellMob[c] = db->TopoSortsReverse[c] - t.sourceSegment;
    //     if (!calLock[c])
    //         mobility[cellMob[c]].insert(c);
    //     topoOrder[c] = beforTopo;
    //     maintainRev(c, beforTopo, afterTopo);
    // }
    // for (int i = topoReverseTrace.size() - 1; i >= 0; --i)
    // {
    //     uint32_t c = topoReverseTrace[i].first;
    //     db->TopoSortsReverse[c] = topoReverseTrace[i].second;
    //     bool t = false;
    //     if (mobility[cellMob[c]].find(c) != mobility[cellMob[c]].end())
    //         t = true;
    //     mobility[cellMob[c]].erase(c);
    //     cellMob[c] = db->TopoSortsReverse[c] - topoOrder[c];
    //     if (t)
    //         mobility[cellMob[c]].insert(c);
    // }
    for (int i = pt.size() - 1; i >= 0; --i)
    {
        uint32_t c = pt[i].cell;
        if (pt[i].backward)
        {
            // bool t = false;
            // if (mobility[cellMob[c]].find(c) != mobility[cellMob[c]].end())
            //     t = true;
            mobility[cellMob[c]].erase(c);
            int afterTopo = topoOrder[c];
            topoOrder[c] = pt[i].beforeTopo;
            cellMob[c] = pt[i].cellmob;
            if (!trash[c])
                // if (t)
                mobility[cellMob[c]].insert(c);
            maintainRev(c, topoOrder[c], afterTopo);
        }
        else
        {
            db->TopoSortsReverse[c] = pt[i].beforeTopo;
            // bool t = false;
            // if (mobility[cellMob[c]].find(c) != mobility[cellMob[c]].end())
            //     t = true;
            mobility[cellMob[c]].erase(c);
            cellMob[c] = pt[i].cellmob;
            if (!trash[c])
                // if (t)
                mobility[cellMob[c]].insert(c);
        }
    }
    trace.pop_back();
    return true;
}

bool scheduling::propagarting(uint32_t c, int bfTopo, int cmob)
{
    std::queue<uint32_t> q;
    q.push(c);
    std::vector<bool> lock;
    lock.resize(hCells.size(), false);
    pt.clear();
    pt.emplace_back(progaTrace{true, c, bfTopo, cmob});
    while (!q.empty())
    {
        uint32_t i = q.front();
        q.pop();
        if (lock[i])
            continue;
        lock[i] = true;
        for (auto &netID : db->incidentNets[i])
        {
            auto &net = hNets[netID];
            uint32_t sourceCell = netPins[net.first];
            if (sourceCell == i) // i is source node
            {
                for (auto s = net.first; s != net.end; ++s)
                {
                    int cid = netPins[s];
                    if (cid == i && lock[cid])
                        continue;
                    auto &cell = hCells[cid];
                    if (topoOrder[cell.ID] <= topoOrder[i])
                    {
                        if (cellMob[cell.ID] > 0)
                        {
                            if (db->TopoSortsReverse[cell.ID] - topoOrder[i] - 1 < 0)
                            {
                                cancelOnePropagation();
                                return false;
                            }
                            else
                            {
                                // double ss = 0;
                                // for (int i = 0; i < topoPartNum.size(); ++i)
                                // {
                                //     ss += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
                                // }
                                // checknet.insert(net.ID);
                                int beforeTopo = topoOrder[cell.ID];
                                pt.emplace_back(progaTrace{true, cell.ID, beforeTopo, cellMob[cell.ID]});
                                // bool t = false;
                                // if (mobility[cellMob[cell.ID]].find(cell.ID) != mobility[cellMob[cell.ID]].end())
                                //     t = true;
                                mobility[cellMob[cell.ID]].erase(cell.ID);
                                cellMob[cell.ID] = db->TopoSortsReverse[cell.ID] - topoOrder[i] - 1;
                                if (!trash[cell.ID])
                                    // if (t)
                                    mobility[cellMob[cell.ID]].insert(cell.ID);

                                double tg1 = calOriginTopoGain(beforeTopo, cell.segment);
                                topoOrder[cell.ID] = topoOrder[i] + 1;
                                double tg2 = calOtherTopoGain(topoOrder[cell.ID], cell.segment);
                                trace.back().first += (tg2 + tg1);
                                trace.back().second.emplace_back(Gain<double>(cell.ID, beforeTopo, topoOrder[cell.ID], tg2 + tg1));
                                maintain(cell.ID, beforeTopo, topoOrder[cell.ID]);
                                q.push(cell.ID);
                                // SD2 = 0;
                                // for (int i = 0; i < topoPartNum.size(); ++i)
                                // {
                                //     SD2 += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
                                // }
                                // if (tg1 + tg2 + SD2 != ss)
                                //     printf("popa——Actual gain :%.16f, cal gain =  %.16f \n\n", tg1 + tg2, ss - SD2);
                            } // the drain nodes of i should be propagarated since their topological order has changed
                        }
                        else
                        {
                            cancelOnePropagation();
                            return false;
                        }
                    }
                }
            }
            else // i is drain node
            {
                if (topoOrder[i] <= cellMob[sourceCell] + topoOrder[sourceCell])
                {
                    if (cellMob[sourceCell] > 0)
                    {
                        if (topoOrder[i] - 1 - topoOrder[sourceCell] < 0)
                        {
                            cancelOnePropagation();
                            return false;
                        }
                        else
                        {
                            pt.emplace_back(progaTrace{false, sourceCell, db->TopoSortsReverse[sourceCell], cellMob[sourceCell]});
                            db->TopoSortsReverse[sourceCell] = topoOrder[i] - 1;
                            // bool t = false;
                            // if (mobility[cellMob[sourceCell]].find(sourceCell) != mobility[cellMob[sourceCell]].end())
                            //     t = true;
                            auto &it = mobility[cellMob[sourceCell]];
                            it.erase(sourceCell);
                            cellMob[sourceCell] = topoOrder[i] - 1 - topoOrder[sourceCell];
                            if (!trash[sourceCell])
                                // if (t)
                                mobility[cellMob[sourceCell]].insert(sourceCell); // the source nodes of i only update their mobility
                        }
                    }
                    else
                    {
                        cancelOnePropagation();
                        return false;
                    }
                }
            }
        }
    }
    return true;
}
bool scheduling::maintain(uint32_t c, int beforeTopo, int afterTopo)
{
    double tmpGain = 0;
    double oldPartNum = (double)topoPartNum[beforeTopo][hCells[c].segment];
    double newPartNum = oldPartNum - 1;
    double oldNum = (double)topoTotalNum[beforeTopo];
    double newNum = oldNum - 1;
    double oldAve = topoAve[beforeTopo];
    double newAve = topoAve[beforeTopo] - (double)(1.0 / (double)segNum);

    tmpGain = (double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
              (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
              (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum);
    topoVar[beforeTopo] += tmpGain;
    topoTotalNum[beforeTopo]--; // maintain the origin topo order info.
    topoAve[beforeTopo] -= (double)(1.0 / (double)segNum);
    topoPartNum[beforeTopo][hCells[c].segment]--;

    oldPartNum = (double)topoPartNum[afterTopo][hCells[c].segment];
    newPartNum = oldPartNum + 1;
    oldNum = (double)topoTotalNum[afterTopo];
    newNum = oldNum + 1;
    oldAve = topoAve[afterTopo];
    newAve = topoAve[afterTopo] + (double)(1.0 / (double)segNum);
    tmpGain = 0;
    tmpGain = (double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
              (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
              (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum);
    topoVar[afterTopo] += tmpGain;
    topoTotalNum[afterTopo]++; // maintain the changed topo order info.
    topoAve[afterTopo] += (double)(1.0 / (double)segNum);
    topoPartNum[afterTopo][hCells[c].segment]++;
    return true;
}
bool scheduling::maintainRev(uint32_t c, int beforeTopo, int afterTopo)
{
    double tmpGain = 0;
    double oldPartNum = (double)topoPartNum[beforeTopo][hCells[c].segment];
    double newPartNum = oldPartNum + 1;
    double oldNum = (double)topoTotalNum[beforeTopo];
    double newNum = oldNum + 1;
    double oldAve = topoAve[beforeTopo];
    double newAve = topoAve[beforeTopo] + (double)(1.0 / (double)segNum);

    tmpGain = (double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
              (newPartNum - newAve) * (newPartNum - newAve) / (double)(segNum) -
              (oldPartNum - newAve) * (oldPartNum - newAve) / (double)(segNum);
    topoVar[beforeTopo] += tmpGain;

    topoTotalNum[beforeTopo]++; // maintain the origin topo order info.
    topoAve[beforeTopo] += (double)(1.0 / (double)segNum);
    topoPartNum[beforeTopo][hCells[c].segment]++;

    oldPartNum = (double)topoPartNum[afterTopo][hCells[c].segment];
    newPartNum = oldPartNum - 1;
    oldNum = (double)topoTotalNum[afterTopo];
    newNum = oldNum - 1;
    oldAve = topoAve[afterTopo];
    newAve = topoAve[afterTopo] - (double)(1.0 / (double)segNum);

    tmpGain = (double)(1.0 / (double)segNum) * (double)(1.0 / (double)segNum) +
              (newPartNum - newAve) * (double)(newPartNum - newAve) / (double)(segNum) -
              (oldPartNum - newAve) * (double)(oldPartNum - newAve) / (double)(segNum);
    topoVar[afterTopo] += tmpGain;
    topoTotalNum[afterTopo]--; // maintain the changed topo order info.
    topoAve[afterTopo] -= (double)(1.0 / (double)segNum);
    topoPartNum[afterTopo][hCells[c].segment]--;
    return true;
}

bool scheduling::traceBack()
{
    double maxGain = std::numeric_limits<double>::min();
    int maxIndex = -1;
    double gainSum = 0;
    for (int i = 0; i < trace.size(); ++i)
    {
        gainSum += trace[i].first;
        if (gainSum >= maxGain)
        {
            maxGain = gainSum;
            maxIndex = i;
        }
    }
    for (int i = maxIndex + 1; i < trace.size(); i++)
    {
        for (int j = trace[i].second.size() - 1; j >= 0; --j)
        {
            uint32_t c = trace[i].second[j].cellID;
            int beforTopo = trace[i].second[j].sourceSegment;
            int afterTopo = trace[i].second[j].targetSegment;
            topoOrder[c] = beforTopo;
            maintainRev(c, beforTopo, afterTopo);
        }
        // legalization();
    }
    printf("-    Scheduling Gain is: %f\n", maxGain);
    return true;
}

bool scheduling::legalization()
{
    printf("-    Prepare to proceeding the Legalization...\n");
    for (auto &n : hNets)
    {
        uint32_t s = netPins[n.first];
        for (auto i = n.first; i != n.end; ++i)
        {
            uint32_t cid = netPins[i];
            if (cid == s)
                continue;
            // if (topoOrder[cid] <= topoOrder[s])
            //     printf("In net %d, source node %d's topo is: %d, drain node %d's topo is:%d.\n",
            //            n.ID, s, topoOrder[s], cid, topoOrder[cid]);
            // return false;
        }
    }
    printf("-    Legalization successful!\n");
    return true;
}
bool scheduling::schedulingFlow()
{
    bool ret = true;

    init();
    for (int i = 1; i < mobility.size(); ++i)
    {
        if (mobility[i].size() == 0)
            continue;
        int sum = mobility[i].size();
        while (!mobility[i].empty())
        {
            int c = (*mobility[i].begin());
            if (calLock[c])
                continue;
            int bestTopo = -1;
            if (mobility[i].size() < sum / 50)
            {
                // printf(" Mobility %d : sum process :%d, actual process =  %ld \n", i, sum, mobility[i].size());
                sum = mobility[i].size();
            }
            bestTopo = obtainGain(c);
            if (bestTopo == -1) // calculate the TOB gain of different topological order, select the best topological order
            {
                calLock[c] = true;
                mobility[i].erase(c);
                trash[c] = true;
                continue;
            }
            // double ss = 0;
            // for (int j = 0; j < topoPartNum.size(); ++j)
            // {
            //     ss += (double)(topoVar[j] / (double)(topoTotalNum[j]) / (double)(topoTotalNum[j]));
            // }
            int bfTopo = topoOrder[c];
            int cmob = cellMob[c];
            maintain(c, bfTopo, bestTopo);
            topoOrder[c] = bestTopo; // assign the best topo to the node

            mobility[cellMob[c]].erase(c);
            cellMob[c] = db->TopoSortsReverse[c] - topoOrder[c];
            mobility[cellMob[c]].insert(c);
            if (cellMob[c] <= 0)
                calLock[c] = true; // the node can not be visited

            // double SD1 = 0;
            // for (int j = 0; j < topoPartNum.size(); ++j)
            // {
            //     SD1 += (double)(topoVar[j] / (double)(topoTotalNum[j]) / (double)(topoTotalNum[j]));
            // }
            // if (fabs(trace.back().first + SD1 - ss) > EPS)
            //     printf(" Actual gain :%.16f, cal gain =  %.16f \n\n", trace.back().first, ss - SD1);
            bool r = propagarting(c, bfTopo, cmob);
            if (!r)
            {
                mobility[cellMob[c]].erase(c);
                trash[c] = true;
            }
            // legalization();
            //  if (r)
            //  {
            //      double SD2 = 0;
            //      for (int j = 0; j < topoPartNum.size(); ++j)
            //      {
            //          SD2 += (double)(topoVar[j] / (double)(topoTotalNum[j]) / (double)(topoTotalNum[j]));
            //      }
            //      // if (trace.back().first + SD2 != LASTSD2)
            //      if (fabs(trace.back().first + SD2 - ss) > EPS)
            //          printf("AFTER Actual gain :%.16f, cal gain =  %.16f \n\n", trace.back().first, SD2 - ss);
            //  }
            //  LASTSD2 = SD2;
        }
    }
    // legalization();
    traceBack();
    // double SD2 = 0;
    // for (int j = 0; j < topoPartNum.size(); ++j)
    // {
    //     SD2 += (double)(topoVar[j] / (double)(topoTotalNum[j]) / (double)(topoTotalNum[j]));
    // }
    // printf("Final1 Variance Metric =  %.8f \n\n", SD2);

    topoPartNum.clear();
    topoAve.clear();
    topoTotalNum.clear();
    topoVar.clear();
    topoPartNum.resize(db->topoOrderMAX_ + 1);
    topoAve.resize(db->topoOrderMAX_ + 1);
    topoTotalNum.resize(db->topoOrderMAX_ + 1);
    topoVar.resize(db->topoOrderMAX_ + 1);

    for (auto &t : topoPartNum)
    {
        t.resize(segNum, 0);
    }
    for (auto c : hCells)
    {
        ++topoPartNum[topoOrder[c.ID]][c.segment];
        ++topoTotalNum[topoOrder[c.ID]];
    }
    for (int t = 0; t < topoTotalNum.size(); ++t)
        topoAve[t] = (double)topoTotalNum[t] / (double)segNum;
    double SD = 0;
    for (int i = 0; i < topoPartNum.size(); ++i)
    {
        for (int j = 0; j < topoPartNum[i].size(); ++j)
        {
            topoVar[i] += (((double)topoPartNum[i][j] - topoAve[i]) * ((double)topoPartNum[i][j] - topoAve[i])) / (double)segNum;
        }
        SD += (double)(topoVar[i] / (double)(topoTotalNum[i]) / (double)(topoTotalNum[i]));
    }
    printf("-    After Scheduling Variance Metric =  %.8f \n\n", SD);
    legalization();
    // calculateTimeStep();
    return ret;
}
