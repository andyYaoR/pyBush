from util import PriorityQueue

def reconstruct_path(came_from, start, goal):
    try:
        current = goal
        path = []
        while current != start:
            path.append(current)
            if current not in came_from:
                print('SPT error from {} to {} via {}'.format(start, goal, current))
            current = came_from[current]
        path.append(start) # optional
        path.reverse() # optional
        return path
    except Exception:
        print(came_from)
        print('SPT error from {} to {}'.format(start, goal))
        raise Exception('SPT error from {} to {}'.format(start, goal))



""" One-to-ALL Shortest Paths using Dijkstra Algorithm Implemented with Binary Heaps"""
def dijkstra_search_bheap_all(start, pred_dict, edges=None, non_passing_nodes=[], classCost=False):
    frontier = PriorityQueue()
    frontier.put(start, 0)
    came_from = {}
    cost_so_far = {}
    came_from[start] = None
    cost_so_far[start] = 0
    init_nodes = list(pred_dict)

    def passable(nodeID):
        return nodeID not in non_passing_nodes

    while not frontier.empty():
        current = frontier.get()
        neighbors = pred_dict[current]
        #neighbors = filter(passable, neighbors)
        for next in neighbors:
            if not next:
                continue
            if edges:
                if classCost == False:
                    new_cost = cost_so_far[current] + edges[(current, next)]['cost']
                else:
                    new_cost = cost_so_far[current] + edges[(current, next)]['classCost'][classCost]
            else:
                new_cost = cost_so_far[current] + 1

            if next not in cost_so_far or new_cost < cost_so_far[next]:
                cost_so_far[next] = new_cost
                priority = new_cost
                #if next in init_nodes:
                if next not in non_passing_nodes and next in init_nodes:
                    frontier.put(next, priority)
                came_from[next] = current

    return came_from, cost_so_far
