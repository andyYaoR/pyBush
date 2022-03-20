import warnings
import hashlib
import pickle
from network_import import dataset_import
from routing import dijkstra_search_bheap_all, reconstruct_path
from util import topological_ordering, objective_function_value, network_cost_update, bush_scan, \
    bush_pred_to_forward_star, compute_shifted_flow
import time

flow_ignore_threshold = 0  # for links with near-zero flows, these flows will be ignored (removed)
thresholdGap = 1e-12
thresholdAEC = 1e-12
LPRule = 'LONGEST_USED_BUSH_PATH'

def bush_flow_push(ori, network, demand, bush, network_identifier=None, rectify=False, flow_shift=False,
                   classID=None, inner_iteration=None):
    bush['nodeFlow'] = {node: 0 for node in network['nodes']}
    demand_node_flows = {des: demand[ori][des] for des in demand[ori]}
    bush['nodeFlow'].update(demand_node_flows)

    if not flow_shift:
        bush['bushFlow'] = {edge: 0 for edge in network['edges'] if (edge[0] not in network['non_passing_nodes']
                                                                     or edge[0] == ori)}

    backtrack_bush_topology = bush['bush_topology'].copy()
    backtrack_bush_topology.reverse()
    updated_links = []

    for node in backtrack_bush_topology:
        if len(bush['pred'][node]) == 1:
            pred_node = bush['pred'][node][0]
            if pred_node:
                bush['bushFlow'][(pred_node, node)] = bush['nodeFlow'][node]
                bush['nodeFlow'][pred_node] += bush['nodeFlow'][node]  # pushing flows from destination to origin
                if node not in bush['merges']['approach_flows']:
                    bush['merges']['approach_flows'][node] = {}
                bush['merges']['approach_flows'][node][pred_node] = bush['nodeFlow'][node]
                bush['merges']['shortest_pred'][node] = pred_node
                bush['merges']['longest_pred'][node] = pred_node
        else: # this is a merge node
            if rectify:
                total_flows = sum([bush['merges']['approach_flows'][node][pred_node]
                                   for pred_node in bush['merges'][node]])
                if total_flows > 0:
                    for pred_node in bush['merges']['approach_flows'][node]:
                        # In case of changes demands, rescale it by the same approach proportion
                        bush['merges']['approach_flows'][node][pred_node] *= bush['nodeFlow'][node] / total_flows
                else:  # In case of removed demands, push all existing flows to shortest path
                    for pred_node in bush['merges']['approach_flows'][node]:
                        bush['merges']['approach_flows'][node][pred_node] = 0
                    bush['merges']['approach_flows'][node][bush['merges']['shortest_pred']] = bush['nodeFlow'][node]

            if flow_shift and node in bush['merges']['shiftable_links'] \
                    and len(bush['merges']['shiftable_links'][node]['LP_links']) > 0:

                if bush['merges']['shiftable_links'][node]['SP_links'][0][0] != bush['merges']['shiftable_links'][node]['LP_links'][0][0] \
                        or bush['merges']['shiftable_links'][node]['SP_links'][-1][1] != bush['merges']['shiftable_links'][node]['LP_links'][-1][1]:
                    raise Exception('Backtracking failed! {} {}'.format( bush['merges']['shiftable_links'][node]['SP_links'], bush['merges']['shiftable_links'][node]['LP_links']))

                dx = compute_shifted_flow(bush=bush,
                                          network=network,
                                          merge_node=node,
                                          inner_iteration=inner_iteration)


                # applying the shifts
                for link in bush['merges']['shiftable_links'][node]['SP_links']:
                    network['edges'][link]['flow'] += dx
                    network['edges'][link]['classFlow'][classID] += dx
                    updated_links.append(link)
                    bush['bushFlow'][link] += dx
                    if link[1] in bush['merges']['approach_flows']:
                        bush['merges']['approach_flows'][link[1]][link[0]] += dx

                for link in bush['merges']['shiftable_links'][node]['LP_links']:
                    tmp = network['edges'][link]['flow']
                    network['edges'][link]['flow'] -= dx
                    network['edges'][link]['classFlow'][classID] -= dx
                    if network['edges'][link]['flow'] < 0 and abs(network['edges'][link]['flow']) < 1e-10:
                        network['edges'][link]['flow'] = 0
                    if network['edges'][link]['classFlow'][classID] < 0 \
                            and abs(network['edges'][link]['classFlow'][classID]) < 1e-10:
                        network['edges'][link]['classFlow'][classID] = 0
                    updated_links.append(link)
                    bush['bushFlow'][link] -= dx

                    if link[1] in bush['merges']['approach_flows']:
                        bush['merges']['approach_flows'][link[1]][link[0]] -= dx

                        if bush['merges']['approach_flows'][link[1]][link[0]] < 0:
                            raise ValueError('Negative approach flow on link {} {} {}'
                                             .format(link, bush['merges']['approach_flows'][link[1]][link[0]], dx))

                    if network['edges'][link]['flow'] < 0 or network['edges'][link]['classFlow'][classID] < 0 \
                            or bush['bushFlow'][link] < 0:
                        raise ValueError('Negative edge flow on link {} {} {} {} {} {}'
                                         .format(link, network['edges'][link]['flow'],
                                                 network['edges'][link]['classFlow'][classID],
                                                 bush['bushFlow'][link], dx, tmp))


            for pred_node in bush['pred'][node]:
                # the approach flows are computed after update bushes
                approach_flow = bush['merges']['approach_flows'][node][pred_node]
                #if flow_shift:
                #    print('Previous {}'.format((pred_node, node)), bush['bushFlow'][(pred_node, node)], 'After', approach_flow)
                bush['bushFlow'][(pred_node, node)] = approach_flow
                # pushing flows from destination to origin
                bush['nodeFlow'][pred_node] += approach_flow


    if network_identifier:
        bush['network_identifier'] = network_identifier

    if not flow_shift:
        return bush
    else:
        return bush, network, updated_links


def bush_init_worker(network, demands, bushID, network_identifier):
    pred, cost_so_far = dijkstra_search_bheap_all(start=bushID[1],
                                                  pred_dict=network['forward_star'],
                                                  edges=network['edges'],
                                                  non_passing_nodes=network['non_passing_nodes'],
                                                  classCost=bushID[0])
    if len(set(network['nodes']) - set(list(pred)) - set(network['non_passing_nodes'])) != 0:
        raise ValueError('Cannot find initial connected bush!')

    # At initialization, the shortest path tree is acyclic, therefore, each node has only one predecessor
    pred = {node: [pred[node]] for node in pred}

    bush_topology = topological_ordering(bushID[1], network, pred)

    bush = {
            'network': network,
            'LPcost': {},
            'SPcost': {},
            'bushFlow': {},
            'nodeFlow': {},
            'pred': pred,
            'bush_topology': bush_topology,
            'merges': {'shortest_pred': {},
                       'longest_pred': {},
                       'approach_flows': {},
                       'shiftable_links': {}
                       },
            'network_identifier': network_identifier
            }
    bush = bush_flow_push(ori=bushID[1],
                          network=network,
                          demand=demands[bushID[0]], #class-specific demand
                          bush=bush,
                          network_identifier=network_identifier)
    return bush


def bush_initialization(network, demands, previous_bushes=None):

    if not network['Initialized']:
        print('Re-initialization for imported network ...')
        for link in network['edges']:
            network['edges'][link]['flow'] = 0
            network['edges'][link]['cost'] = network['edges'][link]['ff_tt']
            for classID in list(demands):
                network['edges'][link]['classFlow'][classID] = 0
                network['edges'][link]['classCost'][classID] = network['edges'][link]['cost'] \
                                                               + network['edges'][link]['classToll'][classID]
        network['Initialized'] = True

    # For identify if rectify is needed
    network_md5 = hashlib.md5(pickle.dumps(network)).hexdigest()
    demands_md5 = hashlib.md5(pickle.dumps(demands)).hexdigest()

    # Initialize bushes using shortest path trees
    bushes = {}
    rectify = False # only when the bushes is reloaded, or demand has been changed
    if previous_bushes:
        # To do: implement warm-up by loading existing bushes
        #
        #
        raise NotImplementedError('Warm-up not implemented yet!')
    else:
        # Default: create bushes from scratch for each physical origin node and for each class
        bush_index = [(classID, ori) for classID in demands for ori in list(demands[classID])]

        print('\n Creating initial bushes and pushing bush flows ...')
        for bushID in bush_index:
            bushes[bushID] = bush_init_worker(network, demands, bushID, network_md5)

        bushes['md5'] = (network_md5, demands_md5)
        print('-> Created initial bushes and pushed bush flows')

    # If changes in demand or network, rectify the bush flows after reloading the bushes
    if (network_md5, demands_md5) != bushes['md5']:
        bushes = {bushID: bush_flow_push(ori=bushID[1],
                                         network=network,
                                         demand=demands[bushID[0]],
                                         bush=bushes[bushID],
                                         network_identifier=network_md5,
                                         rectify=True
                                         ) for bushID in bushes if bushID != 'md5'}
        bushes['md5'] = (network_md5, demands_md5)

    #print('\n Updating initial network link flows and costs...')
    # After creating the bushes, update the link flows in the network
    for link in network['edges']:
        network['edges'][link]['flow'] = sum([bushes[bushID]['bushFlow'][link] for bushID in bushes
                                              if bushID != 'md5'
                                              if link in bushes[bushID]['bushFlow']])



        for classID in demands:
            network['edges'][link]['classFlow'][classID] = sum([bushes[bushID]['bushFlow'][link] for bushID in bushes
                                                                if bushID != 'md5'
                                                                if bushID[0] == classID
                                                                if link in bushes[bushID]['bushFlow']])

    network = network_cost_update(network=network)

    #print('-> Updated initial network link flows and costs')

    #print('\n Updating gradients...')
    # Could be implemented with flow shifting? Or computes the closed-form part of the gradient here?
    #print('-> Updated gradients')

    return bushes, network

def bush_structure_update_worker(network, demands, bushID, bush, network_identifier):
    # if the network in the bush is outdated, update it with latest network
    if bush['network_identifier'] != network_identifier:
        bush['network'] = network
        bush['network_identifier'] = network_identifier

    # Update LPcost (longest_pred) and SPcost (shortest_pred) in the bush
    bush = bush_scan(start=bushID[1],
                     bush=bush,
                     edges=network['edges'],
                     flow_ignore_threshold=flow_ignore_threshold,
                     classCost=bushID[0],
                     LPRule=LPRule)

    new_links = []
    P1_links = []
    P2_links = []
    for link in network['edges']:
        # For numerical stability, remove near-zero link flows
        if link in bush['bushFlow'] and bush['bushFlow'][link] > 0 and bush['bushFlow'][link] < flow_ignore_threshold:
            bush['bushFlow'][link] = 0

        if (link in bush['bushFlow'] and bush['bushFlow'][link] > 0) or link[1] == bushID[1] \
                or (link[0] in network['non_passing_nodes'] and link[0] != bushID[1]):
            continue  # this link has already been included in the bush

        # otherwise, check if this link provides a shortcut, based on a stricter condition (P1 and P2 in Nie 2010)
        if link[0] not in bush['LPcost'] or link[1] not in bush['LPcost']:
            warnings.warn('Not reachable node {} found in bush of origin {}!'.format(link, bushID[1]))
            continue  # should we raise an exception? since every node should be reachable from the origin

        if bush['SPcost'][link[0]] + network['edges'][link]['classCost'][bushID[0]] < bush['SPcost'][link[1]]:
            P1_links.append(link)

        #if bush['SPcost'][link[0]] + network['edges'][link]['classCost'][bushID[0]] == bush['SPcost'][link[1]] \
        #        and link in bush['bushFlow'] and bush['bushFlow'][link] == 0:
        #    #new_links.append(link)
        #    pass

        if bush['LPcost'][link[0]] + network['edges'][link]['classCost'][bushID[0]] < bush['LPcost'][link[1]]:
            P2_links.append(link)  # Nie (2010) P2, stricter than Bar-Gera (2002), i.e., potentially less links



    # In case no new link is found, and P1 and P2 are not empty, P2 links are used to avoid breakdown
    new_links += list(set(P1_links).intersection(set(P2_links)))
    new_links = set(new_links)


    if LPRule == 'LONGEST_BUSH_PATH' and len(new_links) == 0 and len(P1_links) > 0 and len(P2_links) > 0:
        new_links = P2_links
        warnings.warn('P2 links are directly added for bush of origin {} to avoid breakdown!'.format(bushID[1]))

    # Keep all links in the bush with positive flows
    for link in bush['bushFlow']:
        if bush['bushFlow'][link] > 0:
            new_links.add(link)

    if len(new_links) > 0:
        # Update pred merges
        bush['pred'] = {bushID[1]: [None]}  # recompute the pred
        for link in new_links:
            if link not in bush['bushFlow']: # Make sure existing flow is not override
                bush['bushFlow'][link] = 0

            # Initialize the next_node in bush['pred'] for all nodes
            if link[1] not in bush['pred']:
                bush['pred'][link[1]] = []
            bush['pred'][link[1]].append(link[0])

            # Initialize approach flow for the merge node, in case a new link
            if link[1] not in bush['merges']['approach_flows']:
                bush['merges']['approach_flows'][link[1]] = {}

            # Initialize approach flow for the approach, in case a new link
            if link[0] not in bush['merges']['approach_flows'][link[1]]:
                bush['merges']['approach_flows'][link[1]][link[0]] = 0

        shortest_links = []
        for node in set(network['nodes']) - set(list(bush['pred'])):
            shortest_path = reconstruct_path(bush['merges']['shortest_pred'], bushID[1], node)
            shortest_links += list(zip(shortest_path[:-1], shortest_path[1:]))
        shortest_links = set(shortest_links)

        for link in shortest_links:
            # Initialize the next_node in bush['pred'] for all nodes
            if link[1] not in bush['pred']:
                bush['pred'][link[1]] = [link[0]]
                if link not in bush['bushFlow']:  # Make sure existing flow is not override
                    bush['bushFlow'][link] = 0

            # Initialize approach flow for the merge node, in case a new link
            if link[1] not in bush['merges']['approach_flows']:
                bush['merges']['approach_flows'][link[1]] = {link[0]: 0}


        #tmp = {link[1] for link in new_links}
        #print(bushID[1], len(bush['pred']), len(set(network['nodes'])), set(network['nodes']) - set(list(bush['pred'])), len(tmp), len(bush['merges']['shortest_pred']))
        # Update bush topology
        #fl = {(pred_node, next_node) : bush['bushFlow'][(pred_node, next_node)] for next_node in bush['pred'] for pred_node in bush['pred'][next_node] if pred_node}
        #print('Empty links', {link: fl[link] for link in fl if fl[link] == 0})
        bush['bush_topology'] = topological_ordering(bushID[1], network, bush['pred'])

        # If any update to the bush structure, re-initialize the node and bush flows
        bush = bush_flow_push(ori=bushID[1],
                              network=network,
                              demand=demands[bushID[0]],  # class-specific demand
                              bush=bush,
                              network_identifier=network_identifier)

    return bush

def update_bush_structures(bushes, network, demands):
    bush_index = [(classID, ori) for classID in demands for ori in list(demands[classID])]

    network_identifier = hashlib.md5(pickle.dumps(network)).hexdigest()
    if network_identifier != bushes['md5'][0]:
        prev_network_identifier = network_identifier
        network = network_cost_update(network=network) # for consistency, update network costs with current flows
        network_identifier = hashlib.md5(pickle.dumps(network)).hexdigest()
        #print('\n Updated network ---', prev_network_identifier, bushes['md5'][0], network_identifier)

    #print('\n Updating bush structures ...')
    for bushID in bush_index:
        bushes[bushID] = bush_structure_update_worker(network=network,
                                                      demands=demands,
                                                      bushID=bushID,
                                                      bush=bushes[bushID],
                                                      network_identifier=network_identifier)

    # Update the network has been updated, update the identifier for all the bushes
    bushes['md5'] = (network_identifier, bushes['md5'][1])
    #print('-> Updated bush structures')

    return bushes


def find_shiftable_links(ori, bush):  # only flows between a divergent and merge nodes need to be shifted
    bush['merges']['shiftable_links'] = {}
    backtrack_bush_topology = bush['bush_topology'].copy()
    backtrack_bush_topology.reverse()

    for i in range(len(backtrack_bush_topology)):
        merge_node = backtrack_bush_topology[i]
        if len(bush['pred'][merge_node]) > 1: # this is an actual merge node
            if merge_node not in bush['merges']['longest_pred']:
                continue  # when there is no longest used (or shortest) path incidents the merge_node, continue
            if merge_node not in bush['merges']['shortest_pred']:
                warnings.warn('Merge node in pred list, but not (used) on the shortest path!')

            shortest_path = reconstruct_path(bush['merges']['shortest_pred'], ori, merge_node)
            longest_path = reconstruct_path(bush['merges']['longest_pred'], ori, merge_node)

            common_divergent_nodes = set(shortest_path).intersection(set(longest_path)) - {merge_node}

            SP_links = shortest_path[max([shortest_path.index(node) for node in common_divergent_nodes]):]
            SP_links = list(zip(SP_links[:-1], SP_links[1:]))


            LP_links = longest_path[max([longest_path.index(node) for node in common_divergent_nodes]):]
            LP_links = list(zip(LP_links[:-1], LP_links[1:]))
            #print('Merge node', ori, merge_node, SP_links, LP_links)

            bush['merges']['shiftable_links'][merge_node] = {'SP_links': SP_links, 'LP_links': LP_links}

    return bush


def update_bush_flows(bushes, network, demands, inner_iteration):
    """
    This function should be performed in a sequential manner for better convergence performances
    :param bushes:
    :param network:
    :param demands:
    :return:
    """
    bush_index = [(classID, ori) for classID in demands for ori in list(demands[classID])]
    #print('\n Updating bush flows ...')
    updated = []
    inner_gap = []
    time_updated = False
    for bushID in bush_index:
        # Update LPcost (longest_pred) and SPcost (shortest_pred) in the bush
        #a = time.time()
        bushes[bushID] = bush_scan(start=bushID[1],
                                   bush=bushes[bushID],
                                   edges=network['edges'],
                                   flow_ignore_threshold=flow_ignore_threshold,
                                   classCost=bushID[0],
                                   LPRule=LPRule)
        #print('bush_scan', time.time() - a)

        #a = time.time()
        bushes[bushID] = find_shiftable_links(ori=bushID[1], bush=bushes[bushID])
        #print('find_shiftable_links', time.time() - a)

        # Heuristic to stop processing current bush if the gap is relatively small, and move to next bush
        #a = time.time()
        max_bush_gap = []
        bushSPTT = 0
        bushExcess = 0
        for des in demands[bushID[0]][bushID[1]]:
            try:
                bushSPTT += demands[bushID[0]][bushID[1]][des] * bushes[bushID]['SPcost'][des]
            except Exception:
                raise Exception("{} {} {}".format(bushSPTT, des in demands[bushID[0]][bushID[1]], des in bushes[bushID]['SPcost']))

        nodes = list(bushes[bushID]['SPcost'])
        max_bush_gap += [abs(bushes[bushID]['LPcost'][node] - bushes[bushID]['SPcost'][node])
                         for node in nodes if bushes[bushID]['LPcost'][node] != -999999999]
        bushExcess += sum([bushes[bushID]['bushFlow'][link]
                           * (bushes[bushID]['SPcost'][link[0]]
                              + network['edges'][link]['classCost'][bushID[0]]
                              - bushes[bushID]['SPcost'][link[1]]) for link in bushes[bushID]['bushFlow']])

        #print(max(max_bush_gap), sum(max_bush_gap), (bushExcess/bushSPTT), nodes[max_bush_gap.index(max(max_bush_gap))])
        #print('stop check', time.time() - a)
        # If either condition is satisfied, moved to next bush
        if max(max_bush_gap) < thresholdGap or (bushExcess/bushSPTT) < thresholdAEC:
            continue

        #print(max(max_bush_gap), bushExcess / bushSPTT)
        # Additional flag for stopping the inner iterations
        updated.append(bushID)
        inner_gap.append(max(max_bush_gap))


        # Update bush flow
        #a = time.time()
        bushes[bushID], network, updated_links = bush_flow_push(ori=bushID[1],
                                                               network=network,
                                                               demand=demands[bushID[0]],
                                                               bush=bushes[bushID],
                                                               flow_shift=True,
                                                               classID=bushID[0],
                                                               inner_iteration=inner_iteration)
        #print('bush_flow_push', time.time() - a)

        #a = time.time()

        #print('network_cost_update', time.time() - a)

        if inner_iteration % 1 == 0:
            network = network_cost_update(network=network, updated_links=updated_links)
            time_updated = True

    if not time_updated:
        network = network_cost_update(network=network)

    return bushes, network, updated


def algorithmBush(network_filepath, demand_filepaths, network_type='TNTP', numClasses=1,
                  maxIterations=1000, convergence_gap=1e-4, reload=False):

    start_time = time.time()
    best_known_obj = 0
    reloaded = False
    if not reload:
        print('\n** Importing network **')
        network, demands = dataset_import(network_filepath, network_type, demand_filepaths, numClasses)
        print('-> Imported network')

        print('\n** Start initialization **')
        bushes, network = bush_initialization(network, demands)
        print('\n-> Initialization done!')

        with open('network/Barcelona_flow.txt', 'r') as solution:
            start_reading = False
            flow_loading = {}
            for line in solution:
                if len(line.strip()):
                    if 'From' in line:
                        start_reading = True
                        continue
                    if start_reading:
                        data = line.strip().split('\t')
                        flow_loading[(int(data[0]), int(data[1]))] = float(data[2])

            tmp = {link: 0 for link in network['edges'] if link not in flow_loading}
            flow_loading.update(tmp)

        best_known_obj = objective_function_value(network, flow_loading)

        obj = [objective_function_value(network)]
        print('\nInitial Beckmann objective function value:', obj[-1])

        iteration = 0
        gap = 9999999
        start_inner_iteration = 0

        print(sum([link['flow'] for link in network['edges'].values()]))
    else:
        with open('checkpoint.pkl', 'rb') as handle:
            data = pickle.load(handle)
        bushes = data['bushes']
        network = data['network']
        demands = data['demands']
        iteration = data['iteration']
        gap = data['gap']
        start_inner_iteration = data['inner_iteration'] % 10
        obj = data['obj']
        best_known_obj = data['best_known_obj']

    print('\n** Start iterations **')
    maxIterations = 20
    while iteration < maxIterations and gap > convergence_gap:
        iteration += 1
        print('\n<<< Iteration {} >>>'.format(iteration))
        #gap = 0 # to be accumulated?

        if not reload or (reload and reloaded):
            bushes = update_bush_structures(bushes=bushes,
                                            network=network,
                                            demands=demands)


        # Solve the RMP, with flows sequentially updated, and link costs are only updated once per iteration
        # Consequently, only the gradient is affected, but the links to be updated have been defined
        for inner_iteration in range(start_inner_iteration, 10):
            inner_iteration += 1


            bushes, network, updated = update_bush_flows(bushes=bushes,
                                                         network=network,
                                                         demands=demands,
                                                         inner_iteration=inner_iteration)



            obj.append(objective_function_value(network))
            #print('\nBeckmann objective function value:', len(updated), best_known_obj, obj)
            #print(sum([link['flow'] for link in network['edges'].values()]))

            print('Finished Inner iteration', inner_iteration, 'Benchmark', best_known_obj,
                  'Found', obj[-1], 'Runtime', round(time.time()-start_time, 2))

            if len(updated) == 0:
                break


            checkpoint = {'bushes': bushes,
                          'network': network,
                          'demands': demands,
                          'iteration': iteration,
                          'gap': gap,
                          'inner_iteration': inner_iteration,
                          'obj': obj,
                          'best_known_obj': best_known_obj}
            with open('checkpoint.pkl', 'wb') as handle:
                pickle.dump(checkpoint, handle, protocol=pickle.HIGHEST_PROTOCOL)

            reloaded = True

        print('Finished in {} inner iterations. Benchmark {} : Found {}'
              .format(inner_iteration, best_known_obj, obj[-1]))
        start_inner_iteration = 0





if __name__ == "__main__":
    network_filepath = 'network/Barcelona_net.txt'
    demand_filepaths = {0: 'network/Barcelona_trips.txt'}
    algorithmBush(network_filepath=network_filepath,
                  demand_filepaths=demand_filepaths,
                  numClasses=1,
                  )