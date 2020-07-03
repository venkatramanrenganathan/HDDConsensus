import numpy as np
import matplotlib.pyplot as plt
import networkx as netx
from matplotlib.animation import FuncAnimation
import copy

def animate_nodes(G, node_colors, edge_colors, *args, **kwargs):

    pos = netx.drawing.layout.circular_layout(G)

    # draw graph
    Gnodes = netx.draw_networkx_nodes(G, pos, *args, **kwargs)
    Gedges = netx.draw_networkx_edges(G, pos, *args, **kwargs)
    Glabels = netx.draw_networkx_labels(G, pos)
    # display(Gedges.__dict__)
    plt.axis('off')

    def update(i):
        Gnodes.set_color(node_colors[i])
        Gedges.set_color(edge_colors[i])
        return Gnodes, Gedges,

    fig = plt.gcf()
    animation = FuncAnimation(fig, update, interval=50, frames=len(node_colors), blit=True)
    return animation


n_nodes = 8

"""
Hex color codes for the default nodes and edges, and for malicious nodes and outgoing edges
"""

node_base_color = '#329EDC'
edge_base_color = '#0E0E0E'
node_malicious_color = '#E55151'
edge_malicious_color = '#E1BC45'


A = np.zeros([n_nodes,n_nodes])
A[0,4]=A[0,5]=A[0,6] = 1
A[1,3]=A[1,4]=A[1,5]=A[1,7] = 1
A[2,3]=A[2,4]=A[2,6]=A[2,7] = 1
A[3,1]=A[3,2]=A[3,7] = 1
A[4,0]=A[4,1]=A[4,2]=A[4,5]=A[4,6]=A[4,7] = 1
A[5,0]=A[5,1]=A[5,4]=A[5,6]=A[5,7] = 1
A[6,0]=A[6,2]=A[6,4]=A[6,5]=A[6,7] = 1
A[7,1]=A[7,2]=A[7,3]=A[7,4]=A[7,5]=A[7,6] = 1
G = netx.convert_matrix.from_numpy_matrix(A)
nodelabels = {0:1,1:2,2:3,3:4,4:5,5:6,6:7,7:8}
G = netx.relabel_nodes(G,nodelabels)

n_edges = G.number_of_edges()
node_list = list(G.nodes)
edge_list = list(G.edges)

"""
Number of timesteps for the simulation
Minimum number of malicious nodes set at 0
Maximum number and position of malicious nodes bounded by constraint w.r.t degree of the network and each node
"""

Time_steps = 50
min_attackers = 0
max_attackers = 1
attacker_choice_set = copy.deepcopy(node_list)
#
# display(node_list)
# display(edge_list)

node_colors=[]
edge_colors=[]
#
for t in range(Time_steps):
    node_colors_t = [node_base_color]*n_nodes
    edge_colors_t = [edge_base_color]*n_edges

    """
    Following lines of code choose a random number of attackers between min_attackers and max_attackers from the attacker_choice_set at each time step.
    Can implement update functions on the nodes from here to use malicious nodes over time data from an implemented algorithm
    Each attacker node and connected edges are highlighted.
    Both video and gif formats are generated. Both output formats can be adjusted for speed and fps
    """

    n_malicious = np.random.randint(min_attackers,max_attackers+1)
    if n_malicious > 0:
        malicious_list = np.random.choice(attacker_choice_set, size = n_malicious)
        for m_node in malicious_list:
            node_colors_t[m_node-1] = node_malicious_color
            edge_connect = list(G.edges(m_node))
            for m_edge in edge_connect:
                try:
                    edge_id = edge_list.index(m_edge)
                except:
                    edge_id = edge_list.index(m_edge[::-1])
                edge_colors_t[edge_id] = edge_malicious_color
    node_colors.append(node_colors_t)
    edge_colors.append(edge_colors_t)

animation = animate_nodes(G, node_colors, edge_colors)
animation.save('basic_animation.gif', writer='imagemagick', savefig_kwargs={'facecolor':'white'}, fps=2)
animation.save('basic_video.mp4', fps=2, extra_args=['-vcodec', 'libx264'])
