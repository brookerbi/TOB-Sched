#ifndef TOBSCHEDULING_HPP
#define TOBSCHEDULING_HPP

#include "global.hpp"

class strategy;

class scheduling
{
public:
    // scheduling(strategy *pStrategy)
    // {
    //     strategyInst = pStrategy;
    // }
    scheduling(std::vector<hCell> &c, std::vector<hNet> &n, std::vector<uint32_t> &p) : hCells(c),
                                                                                        hNets(n), netPins(p)
    {
        db = DB::getInst();
    }
    ~scheduling() {}

    // strategy *strategyInst;
    DB *db;

    std::vector<hNet> &hNets;
    std::vector<hCell> &hCells;
    std::vector<uint32_t> &netPins; // index starts from hNet[.].f

    int segNum = 0;
    std::vector<std::set<int>> mobility; // Maybe not the prefect datastructure
    std::vector<int> cellMob;
    std::vector<std::vector<int>> topoPartNum;
    std::vector<double> topoAve;
    std::vector<double> topoVar;
    std::vector<uint32_t> topoTotalNum;

    std::vector<int> topoOrder;
    std::vector<std::pair<double, std::vector<Gain<double>>>> trace;
    std::vector<bool> trash;
    // std::vector<Gain<double>> trace;

    std::vector<bool> calLock;
    struct progaTrace
    {
        bool backward = true;
        uint32_t cell;
        int beforeTopo; // backward is false:db->toposortsrervese; backward is true:beforetopo
        int cellmob;
    };
    std::vector<progaTrace> pt;

    std::vector<bool> tsCalFlag;

    // SparseMap mob;
    const double EPS = 1e-10;
    double SD2 = 0;
    double LASTSD2 = 0;

    bool verifyCellsCritical();

    bool init();

    double calOriginTopoGain(int asapTopo, int segment);

    double calOtherTopoGain(int topo, int segment);

    int obtainGain(uint32_t c);

    bool maintain(uint32_t c, int beforeTopo, int afterTopo);

    bool maintainRev(uint32_t c, int beforeTopo, int afterTopo);

    bool cancelOnePropagation();

    bool propagarting(uint32_t c, int beforeTopo, int cellmob);

    bool updateNeighbor(uint32_t c);

    bool traceBack();

    bool legalization();

    bool schedulingFlow();
};

#endif
