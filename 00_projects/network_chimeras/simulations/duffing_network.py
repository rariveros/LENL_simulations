import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np


def create_fixed_degree_network(num_nodes_x, num_nodes_y, degree_x, degree_y, seed=None):
    # Set the seed for reproducibility
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Define degree sequence
    degree_sequence = [degree_x] * num_nodes_x + [degree_y] * num_nodes_y

    # Check if the sum of degrees is even
    if sum(degree_sequence) % 2 != 0:
        raise ValueError("The sum of degrees must be even to form a valid graph.")

    # Generate a random graph using the configuration model
    G = nx.configuration_model(degree_sequence)

    # Convert the multigraph to a simple graph (remove parallel edges and self-loops)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))

    return G


def create_fixed_degree_network_v2(num_nodes_x, num_nodes_y, degree_x, degree_y, seed=None):
    """
    Creates a fixed degree network with the specified parameters, ensuring the graph is a single connected component.

    Parameters:
    - num_nodes_x (int): Number of nodes of type X.
    - num_nodes_y (int): Number of nodes of type Y.
    - degree_x (int): Degree of each node of type X.
    - degree_y (int): Degree of each node of type Y.
    - seed (int, optional): Random seed for reproducibility.

    Returns:
    - G (networkx.Graph): A simple, connected graph with the specified degree sequence.
    """
    # Set the seed for reproducibility
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Define degree sequence
    degree_sequence = [degree_x] * num_nodes_x + [degree_y] * num_nodes_y

    # Check if the sum of degrees is even
    if sum(degree_sequence) % 2 != 0:
        raise ValueError("The sum of degrees must be even to form a valid graph.")

    while True:
        # Generate a random graph using the configuration model
        G = nx.configuration_model(degree_sequence)

        # Convert the multigraph to a simple graph (remove parallel edges and self-loops)
        G = nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))

        # Check if the graph is connected
        if nx.is_connected(G):
            return G

# Generate a random graph using the Erdős-Rényi model
# random_graph = nx.erdos_renyi_graph(n=num_nodes, p=probability)

eed = 142

N_nodes = 20

degree1 = 3
degree2 = 5

node1_prop = 0.5
node2_prop = 0.5

N_node1 = int(node1_prop*N_nodes)
N_node2 = int(node2_prop*N_nodes)

graph = create_fixed_degree_network(N_node1, N_node2, degree1, degree2, seed=seed)

# Draw the graph
plt.figure(figsize=(8, 6))
nx.draw(graph, with_labels=False, node_color='skyblue', edge_color='gray', node_size=80, font_size=10)
plt.title("Graph")
plt.show()

print("connected components: "+str(nx.number_connected_components(graph)))

# Get the adjacency matrix
adj_matrix = nx.adjacency_matrix(graph).toarray()

# Compute the Laplacian matrix
laplacian_matrix = nx.laplacian_matrix(graph).toarray()

# Visualize the adjacency matrix as a lattice
plt.figure(figsize=(6, 6))
plt.imshow(adj_matrix, cmap='Greys', interpolation='none')
plt.colorbar(label="Edge Weight")
plt.title("Adjacency Matrix Visualization")
plt.xlabel("Node Index")
plt.ylabel("Node Index")
plt.show()

def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate): #implementa rouge-kutta
    t_grid = grids[0]
    x_grid = grids[1]
    y_grid = grids[2]
    fields_history = []
    time_grid = []
    for i in range(Nt - 1):
        old_fields = fields
        k_1 = equations_FD(eq, old_fields, t_grid[i], x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t_grid[i], x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t_grid[i], x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid[i], x_grid, y_grid, parameters, operators)
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields
        if i % t_rate == 0:
            fields_history.append(fields)
            time_grid.append(t_grid[i])
    return fields, fields_history, time_grid

def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators): #ecuaciones
    if eq == 'duffing':
        U = field_slices[0]
        V = field_slices[1]

        alpha = parameters[0]
        mu = parameters[1]
        gamma = parameters[2]
        k = parameters[3]
        w = parameters[4]
        DD = operators[0]

        ddU = Der(DD, U)

        F = V
        G = - U + alpha * U ** 3 - U ** 5 - mu * V + gamma * np.cos(w * t_i) + k * ddU

        fields = np.array([F, G])
    return fields

def Der(D, f): #función de diferenciación
    d_f = D @ f
    return d_f