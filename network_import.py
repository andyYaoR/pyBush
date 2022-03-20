from util import network_strong_connectivity, network_import_postprocessing, link_travel_time

def TNTP_import(network_filepath, demand_filepaths, defaultDemandMultiplier, numClasses, classToll=None):
    if not classToll:
        classToll = {}
    numZones = None
    numNodes = None
    numLinks = None
    defaultTollFactor = 1
    defaultDistanceFactor = 1
    non_passing_nodes = []
    edges = {}

    demands = {}
    demands_factors = {}
    # Read network file
    with open(network_filepath, 'r') as net:
        read_mode = False
        meta_finished = False
        for line in net:
            if len(line.strip()):
                line = line.replace('\t', ' ')
                if 'NUMBER OF ZONES' in line:
                    numZones = int(line.strip().split(' ')[-1].strip())
                if 'NUMBER OF NODES' in line:
                    numNodes = int(line.strip().split(' ')[-1].strip())
                if 'FIRST THRU NODE' in line:
                    non_passing_nodes += list(range(1, int(line.strip().split(' ')[-1].strip())))
                if 'NUMBER OF LINKS' in line:
                    numLinks = int(line.strip().split(' ')[-1].strip())
                if 'DISTANCE FACTOR' in line:
                    defaultDistanceFactor = float(line.strip().split(' ')[-1].strip())
                if 'TOLL FACTOR' in line:
                    defaultDistanceFactor = float(line.strip().split(' ')[-1].strip())
                if '<END OF METADATA>' in line:
                    meta_finished = True
                if meta_finished and '~' in line and ('Init node' in line or 'init_node' in line):
                    read_mode = True
                    continue

                if read_mode:
                    line_data = line.strip().split(' ')[:-1]
                    # edges[(start node, end node)] = {...}
                    edges[(int(line_data[0]), int(line_data[1]))] = {'Capacity': float(line_data[2]),
                                                                     'Length': float(line_data[3]),
                                                                     'ff_tt': float(line_data[4]),
                                                                     'alpha': float(line_data[5]),
                                                                     'beta': float(line_data[6]),
                                                                     'speedLimit': float(line_data[7]),
                                                                     'toll': float(line_data[8]),
                                                                     'linkType': int(line_data[9]),
                                                                     'flow': 0,
                                                                     'cost': link_travel_time(t0=float(line_data[4]),
                                                                                              alpha=float(line_data[5]),
                                                                                              beta=float(line_data[6]),
                                                                                              flow=0,
                                                                                              cap=float(line_data[2])),
                                                                     'classToll': {},
                                                                     'classFlow': {},
                                                                     'classDistance': {},
                                                                     'classCost': {},
                                                                     }
                    # Define the toll for each class, default using the same toll, otherwise should be provided
                    for classID in range(numClasses):
                        if classID not in classToll: # If not provided, use the default value from the network file
                            edges[(int(line_data[0]), int(line_data[1]))]['classToll'][classID] = float(line_data[8])
                        elif type(classToll[classID]) is dict: # If provided for each class and each link
                            if (int(line_data[0]), int(line_data[1])) not in classToll[classID]:
                                raise AttributeError('Link not provided in class toll')
                            edges[(int(line_data[0]), int(line_data[1]))]['classToll'][classID] = \
                                classToll[classID][(int(line_data[0]), int(line_data[1]))]
                        else: # If provided for each class as a uniform value
                            edges[(int(line_data[0]), int(line_data[1]))]['classToll'][classID] = classToll[classID]

        if numZones == None or numNodes == None or numLinks == None:
            raise AttributeError('Metadata for TNTP file not correct! Please check your network file ...')

    # Read demand file
    for classID in range(numClasses): # read class-specific demand files
        demands[classID] = {}
        demands_factors[classID] = {'tollFactor': defaultTollFactor,
                                    'distanceFactor': defaultDistanceFactor,
                                    'demandMultiplier': defaultDemandMultiplier
                                    }
        with open(demand_filepaths[classID], 'r') as demand:
            origin = -1
            for line in demand:
                if len(line.strip()):
                    if 'NUMBER OF ZONES' in line and int(line.strip().split(' ')[-1].strip()) != numZones:
                        raise ValueError('Demand zones not match between demand file and network file')
                    if 'TOTAL OD FLOW' in line:
                        class_total_flows = float(line.strip().split(' ')[-1].strip())
                    if 'DEMAND MULTIPLIER' in line:
                        demands_factors[classID]['demandMultiplier'] = float(line.strip().split(' ')[-1].strip())
                    if 'DISTANCE FACTOR' in line:
                        demands_factors[classID]['distanceFactor'] = float(line.strip().split(' ')[-1].strip())
                    if 'TOLL FACTOR' in line:
                        demands_factors[classID]['tollFactor'] = float(line.strip().split(' ')[-1].strip())
                    if 'Origin' in line:
                        try:
                            origin = int(line.strip().split('\t')[-1].strip())
                        except ValueError:
                            origin = int(line.strip().split(' ')[-1].strip())

                        demands[classID][origin] = {}
                        continue

                    if origin > 0:
                        to_demands = line.strip().split(';')
                        for to in to_demands:
                            if len(to):
                                demands[classID][origin][int(to.split(':')[0].strip())] = \
                                    float(to.split(':')[1].strip()) * demands_factors[classID]['demandMultiplier']

    edges, nodes, forward_star, backward_star = network_import_postprocessing(edges, demands_factors)

    network = {'edges': edges,
               'nodes': nodes,
               'forward_star': forward_star,
               'backward_star': backward_star,
               'non_passing_nodes': non_passing_nodes,
               'Initialized': True #If the network has ever been modified (any attribute) marked as False
               }

    # Check strong connectivity, i.e., there is at least 1 route from any node to all other nodes
    network = network_strong_connectivity(network, demands_factors, numClasses)

    return network, demands




def dataset_import(network_filepath, network_type, demand_filepaths, numClasses, defaultDemandMultiplier=1):
    if type(demand_filepaths) is not dict:
        demand_filepaths = {0: demand_filepaths}

    if numClasses != len(demand_filepaths):
        raise IndexError('Class-specific demand file not match with number of classes!')

    if network_type == 'TNTP':
        return TNTP_import(network_filepath, demand_filepaths, defaultDemandMultiplier, numClasses)

