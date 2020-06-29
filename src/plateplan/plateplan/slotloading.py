
import gurobipy as grb
import numpy as np

from .utils import not_none, flatten


def _make_stock_layout(allowed, targets, nstages, priority=None):
    """Solve the fixed stage Mantis scheduling problem.
    
    Inputs
    ------
    allowed: bool np.ndarray (nslots x nwares)
        If `allowed[i,j]` is True, the labware j is allowed in slot i.
    targets: int np.ndarray, len=nwares
        Number of each labware type that must be scheduled.
    nstages: int
        Number of times the Mantis can be loaded with stocks. The total number
        of available slots is nstages * nslots.
    priority: optional, numeric np.ndarray (nwares x nslots)
        When scheduling, try to maximum the priority by placing labware in
        slots with the highest priority. Not every labware or slot needs a
        priority score. Slots with priority zero will still be filled.
        If None, all priorities are set to zero.
        
    Returns
    -------
    int np.ndarray (nslots, nstages)
    
    An integer layout matrix. If `layout[i,j] = k`, then slot i during stage j
    should be loaded with the labware with index k. Any slot that should not
    be loaded will be -1.
    """
    model = grb.Model()
    model.setParam("OutputFlag", 0)
    varname = lambda stage, slot, ware: "Stage{}Slot{}Ware{}".format(stage, slot, ware)
    nslots, nwares = allowed.shape
    
    if priority is None:
        priority = np.zeros((nwares, nslots))
    
    inds = [[[] for _ in range(nslots)] for __ in range(nstages)]
    for stage in range(nstages):
        for slot in range(nslots):
            inds[stage][slot] = [None for _ in range(nwares)]
            for ware in range(nwares):
                if allowed[slot,ware]:
                    inds[stage][slot][ware] = model.addVar(
                            vtype=grb.GRB.BINARY,
                            obj=priority[ware,slot],
                            name=varname(stage, slot, ware))
            # only one ware can be placed in each slot per stage
            terms = not_none(inds[stage][slot])
            if terms:
                linex = grb.LinExpr()
                linex.addTerms([1 for _ in terms], terms)
                model.addLConstr(linex, sense=grb.GRB.LESS_EQUAL, rhs=1)
    # sum_{ware i} ind[i] = target[i]
    for ware in range(nwares):
        ware_inds = not_none([inds[stage][slot][ware] 
                                for stage in range(nstages)
                                for slot in range(nslots)])
        linex = grb.LinExpr()
        linex.addTerms([1 for _ in ware_inds], ware_inds)
        model.addLConstr(linex, sense=grb.GRB.EQUAL, rhs=targets[ware])
    
    # maximize the sum of priorities for scheduled slots
    model.setAttr(grb.GRB.Attr.ModelSense, grb.GRB.MAXIMIZE)
    
    model.update()
    model.optimize()
    
    if model.getAttr("Status") != grb.GRB.OPTIMAL:
        # infeasible; probably because more stages are required
        return None
    
    layout = -np.ones((nslots, nstages), dtype=int)
    for stage in range(nstages):
        for slot in range(nslots):
            for ware in range(nwares):
                ind = inds[stage][slot][ware]
                if inds[stage][slot][ware]:
                    if ind.getAttr('x') > 0.5:
                        layout[slot,stage] = ware
                    
    return layout


def assign_mantis_stocks(stocks, config):
    """Group stocks into batches for the Mantis.
    
    Inputs
    ------
    stocks: List[Stock]
        List of Stock objects that should be scheduled in batches for the
        Mantis.
    config: dict[str: dict]
        Mantis configuration settings, i.e. where labware can and should be
        placed on the Mantis. The `config` dict contains two sub-dicts:
            'allowed': dict[str: List[str]]
                Each key is a str ID for a slot in the Mantis. The values are
                a list of str IDs for labware that are allowed in this slot.
                The order of the labware IDs in this list does not matter.
            'priority': dict (str : List[List[str]]), optional
                Each key is a labware ID. The values are a list of lists 
                identifying the preferred slots for the labware. Not all the
                slots in 'allowed' need to be assigned a priority. Any slot
                that is not given a priority will be assigned a priority of 
                zero -- it will still be used, but only if all slots with
                higher priority have been filled.
        Example config:
            { 
              'allowed': {
                    'slot1': ['Conical15ml', 'Conical50ml'],
                    'slot2': ['Conical15ml'],
                    'slot3': ['Conical50ml'],
                    'slot4': ['Conical15ml']},
              'priority': {
                    'Conical15ml': [['slot2','slot4'], ['slot1']],
                    'Conical50ml': [['slot3']]}
            }
        In this example, a 15ml tube can be placed in slots 1, 2, or 4, and a
        50 ml tube can be placed in slots 1 or 3. The 15 ml tubes will be 
        placed in slots 2 and 4 first (if possible), followed by slot 1. Any
        50 ml conicals will placed in slot 3 before slot 1 (if possible).
        
    Returns
    -------
    List[dict[str: Stock]]
    
    Stocks will be batched into `n` stages where each stage requires reloading
    the Mantis. The returned list has `n` items, one per stage. Each item is
    a dict with str keys identifying the slots used during the stage mapped
    to the Stock to be loaded.
    """
    nslots = len(config['allowed'])
    slots = list(config['allowed'].keys())
    allowed_wares = set(flatten(config['allowed'].values()))
    stock_wares = set([s.labware.id for s in stocks])
    
    # Ensure every labware in the list of stocks can be used with the Mantis.
    assert stock_wares <= allowed_wares
    # only schedule wares in the stock list
    scheduled_wares = list(stock_wares)
    nwares = len(scheduled_wares)
    
    # By default, all slots are given the same (0) priority. If there is no
    # 'priority' entry in the config, the assignment problem needs to only
    # find a feasible, not optimal, solution.
    priority = np.zeros((nwares,nslots))
    if 'priority' in config:
        for ware in range(nwares):
            priorities = config['priority'].get(scheduled_wares[ware])
            if priorities:
                for i, group in enumerate(priorities):
                    # if a ware's priority list has k levels, the priorities
                    # are k, k-1, ..., 1; and 0 for all other slots
                    obj = len(priorities) - i
                    for member in group:
                        priority[ware,slots.index(member)] = obj
    
    allowed = np.zeros((nslots, nwares), bool)
    for i, slot in enumerate(config['allowed']):
        for ware in config['allowed'][slot]:
            if ware in scheduled_wares:
                allowed[i,scheduled_wares.index(ware)] = True
    
    # find how many of each ware we need to schedule
    countof = lambda w: len([s for s in stocks if s.labware.id == w])
    targets = [countof(ware) for ware in scheduled_wares]
    
    # Our solution uses a sliding horizon approach. To find a scheduling
    # strategy with the minimum number of stages, we increment the number
    # of stages until a feasible solution is found.
    
    # To find an upper bound on the number of stages, we schedule each ware
    # separately using all allowed slots for that ware.
    max_stages = np.ceil(targets / allowed.sum(0)).sum().astype(int)
    # The minimum number of stages would only be feasible if every slot 
    # allowed every ware; however, it's a good starting point.
    nstages = len(stocks) // nslots
    while nstages <= max_stages:
        layout_matrix = _make_stock_layout(allowed, targets, nstages, priority)
        if layout_matrix is not None:
            # we found a feasible solution
            break
        else:
            # no feasible solution; add a stage and try again
            nstages += 1
    else:
        # something is wrong; we should always find a feasible solution
        return None
    
    # For each ware, compile of list of (stage, slot) tuples for final
    # assignment.
    ware_slots = {ware: [] for ware in scheduled_wares}
    for stage in range(nstages):
        for i, slot in enumerate(slots):
            if layout_matrix[i,stage] > -1:
                ware = scheduled_wares[layout_matrix[i,stage]]
                ware_slots[ware].append((stage, slot))
    
    # For each stock, grab an available (stage, slot) location.
    layout = [{} for _ in range(nstages)]
    for stock in stocks:
        stage, slot = ware_slots[stock.labware.id].pop(0)
        layout[stage][slot] = stock
        
    return layout
    

if __name__ == '__main__':
    from pprint import pprint
    import constants
    #pprint(group_mantis_stocks(44, 13, 16, 4, 8))
    #pprint(assign_mantis_stocks({}, constants.mantis_lc3_config))
    stocks = list(CDM_stocks.values())
    stocks[-1].labware = constants.tubes['Bottle250ml']
    pprint(assign_mantis_stocks(stocks, mantis_lc3_config))

        
