import heapq
import warnings
import numpy as np

class PriorityQueue:
    def __init__(self):
        self.elements = []

    def empty(self):
        return len(self.elements) == 0

    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))

    def get(self):
        return heapq.heappop(self.elements)[1]

    def get_cost_item(self):
        return heapq.heappop(self.elements)

    def heap_querry(self):
        return self.elements[0]

    def get_max(self):
        return heapq.nlargest(1, self.elements)[0][0]

from routing import dijkstra_search_bheap_all

def network_strong_connectivity(network, demands_factors, numClasses):
    if len(network['non_passing_nodes']) == 0:
        ori_nodes = [list(network['forward_star'])[-1]]
    else:
        ori_nodes = network['non_passing_nodes']

    for node in ori_nodes:
        forward_came_from, _ = dijkstra_search_bheap_all(start=node,
                                                         pred_dict=network['forward_star'],
                                                         non_passing_nodes=network['non_passing_nodes'])
        backward_came_from, _ = dijkstra_search_bheap_all(start=node,
                                                          pred_dict=network['backward_star'],
                                                          non_passing_nodes=network['non_passing_nodes'])

        forward_not_reachable_nodes = list(set(network['nodes'])
                                           - set(list(forward_came_from)))
        backward_not_reachable_nodes = list(set(network['nodes'])
                                            - set(list(forward_came_from)))

        #print('forward_not_reachable_nodes', node, forward_not_reachable_nodes)
        #print('backward_not_reachable_nodes', node, backward_not_reachable_nodes)

        if len(forward_not_reachable_nodes) > 0 or len(backward_not_reachable_nodes) > 0:
            for to_node in forward_not_reachable_nodes:
                new_link = (node, to_node)
                network['edges'][new_link] = {'Capacity': 1,
                                              'Length': 9999,
                                              'ff_tt': 9999,
                                              'alpha': 0,
                                              'beta': 1,
                                              'speedLimit': 0,
                                              'toll': 0,
                                              'linkType': 1,
                                              'flow': 0,
                                              'cost': 9999,
                                              'classToll': {claasID: 9999 for claasID in range(numClasses)},
                                              'classFlow': {claasID: 0 for claasID in range(numClasses)},
                                              'classDistance': {claasID: 9999 for claasID in range(numClasses)},
                                              'classCost': {claasID: 9999 for claasID in range(numClasses)},
                                              'virtual': True
                                             }

            for to_node in backward_not_reachable_nodes:
                new_link = (node, to_node)
                network['edges'][new_link] = {'Capacity': 1,
                                              'Length': 9999,
                                              'ff_tt': 9999,
                                              'alpha': 0,
                                              'beta': 1,
                                              'speedLimit': 0,
                                              'toll': 0,
                                              'linkType': 1,
                                              'flow': 0,
                                              'cost': 9999,
                                              'classToll': {claasID: 9999 for claasID in range(numClasses)},
                                              'classFlow': {claasID: 0 for claasID in range(numClasses)},
                                              'classDistance': {claasID: 9999 for claasID in range(numClasses)},
                                              'classCost': {claasID: 9999 for claasID in range(numClasses)},
                                              'virtual': True
                                             }

            network['edges'], network['nodes'], network['forward_star'], network['backward_star'] = \
                network_import_postprocessing(network['edges'], demands_factors)

    return network

def network_import_postprocessing(edges, demands_factors):
    forward_star = {}
    backward_star = {}
    nodes = set()
    for link in edges:
        forward_star.setdefault(link[0], []).append(link[1]) # All outgoing nodes (links) from node_0
        backward_star.setdefault(link[1], []).append(link[0]) # All incoming nodes (links) from node_1
        nodes.add(link[0])
        nodes.add(link[1])
        for classID in demands_factors:
            edges[link]['classToll'][classID] = edges[link]['classToll'][classID] \
                                                * demands_factors[classID]['tollFactor']
            edges[link]['classFlow'][classID] = 0
            edges[link]['classDistance'][classID] = edges[link]['Length'] * demands_factors[classID]['distanceFactor']
            # Initialize the class cost
            # The first term is the flow-dependent term (initialized with travel time)
            # The second term is the fixed cost with respect to the each class (could be scaled by class 'tollFactor')
            edges[link]['classCost'][classID] = edges[link]['cost'] + edges[link]['classToll'][classID]

    return edges, list(nodes), forward_star, backward_star


def topological_ordering(ori, network, pred):
    indegree = {node: 1 for node in network['nodes']}
    indegree_bush_nodes = {node: len(pred[node]) for node in pred if pred[node][0] != None}
    indegree.update(indegree_bush_nodes)
    indegree[ori] = 0

    order_list = [node for node in indegree if indegree[node] == 0] # actually could be simplified as origin only

    bush_topology = [] # topology of the bush
    while len(order_list) > 0:
        current = order_list.pop(0) # obtain the first element in the order_list
        bush_topology.append(current)
        if current not in network['forward_star']:  # only incoming links
            continue
        for next in network['forward_star'][current]: # check all outgoing links
            if next in pred and current in pred[next]:
                #print(current, next, pred[next], indegree[next])
                indegree[next] -= 1
                if indegree[next] == 0: # in case more than one incoming link, the topology of each node is processed once
                    order_list.append(next)

    #print('indegree', [(node, indegree[node]) for node in indegree if indegree[node] > 0])
    if len(bush_topology) < len(set(network['nodes']) - set(network['non_passing_nodes'])):
        raise AssertionError('Bush given for topology ordering contains a cycle {} {} {}'.format(len(bush_topology),
                                                                                              len(set(set(network['nodes']))),
                                                                                              len(network['non_passing_nodes'])))
    #print(ori, set(bush_topology) == set(network['nodes']), len(bush_topology), len(set(network['nodes'])), len(network['non_passing_nodes']))
    return bush_topology

def link_travel_time(t0, alpha, beta, flow, cap):
    return t0 * (1 + alpha * np.power((flow/max(cap, 1e-8)), beta))

def objective_function_value(network, flow_loading=None):
    if not flow_loading:
        return sum([(edge['ff_tt']*edge['flow'])
                     * (1 + ((edge['alpha']*np.power((edge['flow']/max(edge['Capacity'], 1e-8)), edge['beta'])/(1+edge['beta']))))
                     for edge in network['edges'].values()])
    else:
        return sum([(network['edges'][link]['ff_tt']*flow_loading[link])
                     * (1 + ((network['edges'][link]['alpha']
                              * np.power((flow_loading[link]/max(network['edges'][link]['Capacity'], 1e-8)), network['edges'][link]['beta']))
                             /(1+network['edges'][link]['beta'])))
                     for link in network['edges']])

def network_cost_update(network, updated_links=None):
    if not updated_links:
        updated_links = network['edges']


    for link in updated_links:
        network['edges'][link]['cost'] = link_travel_time(t0=network['edges'][link]['ff_tt'],
                                                          alpha=network['edges'][link]['alpha'],
                                                          beta=network['edges'][link]['beta'],
                                                          flow=network['edges'][link]['flow'],
                                                          cap=network['edges'][link]['Capacity'])
        for classID in list(network['edges'][link]['classToll']):
            network['edges'][link]['classCost'][classID] = network['edges'][link]['cost'] \
                                                           + network['edges'][link]['classToll'][classID]

    return network


def bush_pred_to_forward_star(pred):
    forward_star = {}

    for node in pred:
        for pred_node in pred[node]:
            if pred_node != None:
                forward_star.setdefault(pred_node, []).append(node)

    return forward_star


def bush_scan(start, bush, edges, flow_ignore_threshold, classCost=False, LPRule='LONGEST_USED_BUSH_PATH'):
    bush['LPcost'] = {start: 0}
    bush['SPcost'] = {start: 0}
    bush['merges']['shortest_pred'] = {}
    bush['merges']['longest_pred'] = {}
    assert LPRule in ['LONGEST_BUSH_PATH', 'LONGEST_USED_BUSH_PATH', 'LONGEST_USED_BUSH_PATH_OR_SP'], \
        'Unknown longest route rule!'

    for node in bush['bush_topology']: # traverse the bush by topological order
        if node not in bush['LPcost']:
            bush['LPcost'][node] = -1 * np.inf
        if node not in bush['SPcost']:
            bush['SPcost'][node] = np.inf

        for pred_node in bush['pred'][node]: # check all approaches (incoming links)
            added = False
            l_added = False
            if pred_node != None:
                if classCost == False:
                    shortest_new_cost = bush['SPcost'][pred_node] + edges[(pred_node, node)]['cost']
                    longest_new_cost = bush['LPcost'][pred_node] + edges[(pred_node, node)]['cost']
                else:
                    shortest_new_cost = bush['SPcost'][pred_node] + edges[(pred_node, node)]['classCost'][classCost]
                    longest_new_cost = bush['LPcost'][pred_node] + edges[(pred_node, node)]['classCost'][classCost]

                #print('LSP', pred_node, bush['SPcost'][pred_node], bush['LPcost'][pred_node])

                try:
                    if node not in bush['SPcost'] or (node in bush['SPcost'] and shortest_new_cost < bush['SPcost'][node]):
                        bush['SPcost'][node] = shortest_new_cost
                        bush['merges']['shortest_pred'][node] = pred_node
                        added = True
                except Exception:
                    raise Exception("{} {} {} {}".format(shortest_new_cost, bush['SPcost'][node], bush['LPcost'][pred_node], edges[(pred_node, node)]['classCost'][classCost]))

                # Default use longest route in the bush, regardless if it carries flows
                if bush['LPcost'][node] == -1 * np.inf or (node in bush['LPcost'] and longest_new_cost > bush['LPcost'][node]):
                    if bush['LPcost'][node] == -1 * np.inf:
                        bush['LPcost'][node] = longest_new_cost
                        bush['merges']['longest_pred'][node] = pred_node
                        l_added = True
                    elif LPRule == 'LONGEST_BUSH_PATH':
                        bush['LPcost'][node] = longest_new_cost
                        bush['merges']['longest_pred'][node] = pred_node
                    # if LONGEST_USED_BUSH_PATH, then only approaches (incoming links) with flows are considered
                    elif LPRule == 'LONGEST_USED_BUSH_PATH' \
                            and bush['merges']['approach_flows'][node][pred_node] > flow_ignore_threshold:
                        bush['LPcost'][node] = longest_new_cost
                        bush['merges']['longest_pred'][node] = pred_node
                        l_added = True
                    elif LPRule == 'LONGEST_USED_BUSH_PATH_OR_SP' and \
                            (bush['merges']['approach_flows'][node][pred_node] > 0
                             or bush['merges']['shortest_pred'][node] == pred_node):
                        bush['LPcost'][node] = longest_new_cost
                        bush['merges']['longest_pred'][node] = pred_node

    return bush


def BPR_derivative(link):
    #print('In Derivative {} {} {}'.format(link['flow'], link['Capacity'], link['beta']))
    if link['flow'] == 0 and (link['beta']-1) < 0:
        return 0.0
    else:
        return link['ff_tt'] * link['alpha'] * link['beta'] * np.power((link['flow']/max(link['Capacity'], 1e-8)), (link['beta']-1)) / max(link['Capacity'], 1e-8)


def compute_shifted_flow(bush, network, merge_node, inner_iteration):
    # Apply simple projected Newton
    # The network costs are computed once for each outer iteration, therefore the SPcost and LPcost are directly used

    newton_step = 1
    min_derivative = 1e-20

    divergent_node = bush['merges']['shiftable_links'][merge_node]['SP_links'][0][0]

    g = (bush['LPcost'][merge_node] - bush['LPcost'][divergent_node]) \
        - (bush['SPcost'][merge_node] - bush['SPcost'][divergent_node])

    h_links = set(bush['merges']['shiftable_links'][merge_node]['SP_links']) \
        .union(set(bush['merges']['shiftable_links'][merge_node]['LP_links']))

    h = sum([BPR_derivative(network['edges'][link]) for link in h_links])

    if h == 0:
        h = min_derivative

    dx = max(0, (newton_step/1) * g / h)

    if dx < 0:
        warnings.warn('Something wrong with negative projected newton')

    #print('\n', dx, bush['SPcost'][merge_node], bush['LPcost'][merge_node])
    #print('LP', [(link, bush['bushFlow'][link]) for link in
    #                      bush['merges']['shiftable_links'][merge_node]['LP_links']])

    if len(bush['merges']['shiftable_links'][merge_node]['LP_links']) > 0:
        dx = min(dx, min([bush['bushFlow'][link] for link in
                          bush['merges']['shiftable_links'][merge_node]['LP_links']]))

    #print('dx after LP', dx)




    return dx
