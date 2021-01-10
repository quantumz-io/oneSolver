import qubogen
import networkx as nx
import dimod
import sys

__copyright__ = "Piotr Gawron"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"

def bqm_to_qubo(bqm):
    max_nodes = len(bqm.linear.keys())
    num_nodes = max(bqm.linear.keys()) + 1
    num_links = len(bqm.quadratic.keys())

    s = "c qubo Target MaxNodes NumNodes NumLinks\n"
    s += f"p qubo 0 {max_nodes} {num_nodes} {num_links}\n"
    s += bqm.to_coo()
    return s

def generate_mvc(graph_nodes):
    g = qubogen.Graph.from_networkx(nx.random_graphs.complete_graph(graph_nodes)) 
    mvc_problem = dimod.AdjArrayBQM(qubogen.qubo_mvc(g), dimod.BINARY)
    return mvc_problem

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: qubo_mvp_generator.py graph_nodes")
        print("see: https://github.com/dwavesystems/dwavebinarycsp")
        exit(-1)
    graph_nodes = int(sys.argv[1])
    
    qubo = generate_mvc(graph_nodes)
    s = bqm_to_qubo(qubo)
    print(s)

