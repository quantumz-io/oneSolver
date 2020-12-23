import dwavebinarycsp
import dimod
import sys

__copyright__ = "hyperQ – Ewa Hendzel"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"

def bqm_to_qubo(bqm):
    max_nodes = len(bqm.linear.keys())
    num_nodes = max(bqm.linear.keys())
    num_links = len(bqm.quadratic.keys())

    s = "c qubo Target MaxNodes NumNodes NumLinks\n"
    s += f"p qubo 0 {max_nodes} {num_nodes} {num_links}\n"
    s += bqm.to_coo()
    return s

def generate_bqm(csp_vars, csp_clauses):
    csp = dwavebinarycsp.factories.random_2in4sat(csp_vars, csp_clauses) 
    bqm = dwavebinarycsp.stitch(csp)
    return bqm

 
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: qubo_csp_generator.py csp_vars csp_clauses")
        print("see: https://github.com/dwavesystems/dwavebinarycsp")
        exit(-1)
    csp_vars = int(sys.argv[1])
    csp_clauses = int(sys.argv[2])
    bqm = generate_bqm(csp_vars, csp_clauses)
    s = bqm_to_qubo(bqm)
    print(s)