import networkx as nx  # version == 2.3
import pandas as pd
import requests  # HTTP Client for Python
import json  # Standard JSON library
from py2cytoscape import util as cy
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np
from colour import Color

#   loading external data

#   load disease genes associated in experiments
with open('./addData/CADgens_DisGenNET.txt', 'r') as file:
    CADgens = file.readline().strip().split(' ')
file.close()

#   load interactome from String and convert to nx Graph
interactome = pd.read_csv('./baseData/String_interactome.csv')
G = nx.from_pandas_edgelist(interactome, 'Source', 'Target')


def target_list_creator(miR_name, miR_names, miR_dict):
    """
    miRTarBase содержит разные формы miRNA, так например, miR-21- может соответствовать
    miR-21-3p и miR-21-5p. Функция конкатенирует таргеты всех форм микроРНК и удалаяет дубли,
    которые возникают из-за того, что одной и тойже мишени может соответствовать несколько строк из-за
    разные методов её подтверждения.
    """

    res = set()
    mir_name_app = []

    for name in miR_names:
        if miR_name in name:
            res.update(miR_dict[name])
            mir_name_app.append(name)

    if not mir_name_app:
        print('miRNA', '"{}"'.format(miR_name), 'not found, use another name')
        return

    print('I found a miRNA with name:', *mir_name_app)
    print('and ', len(res), 'unique targets')

    return res


def tissue_selector():
    ans = str(input('"Human Protein Atlas"(0) or "GTEx"(1) ? '))
    if ans == '1':
        dt = pd.read_csv(
            '/Users/german/Dropbox/CardioCenter/Papers/2020/375_РФФИ/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
            sep='\t')  # загрузка med(TPM) from GTEx

        print('Gene universe is...')

        labels = sorted(dt.columns)
        index = [i for i in range(0, len(dt.columns))]
        for i in index:
            print(index[i], '-----', labels[i])

        print("give me int, or input 'all'", sep='\n')
        tissue_id = str(input('Your choice: '))

        if tissue_id != 'all':
            tissue = labels[int(tissue_id)]
            tissue_genes = set(dt['Description'][dt[tissue] > 0])
            print('your tissue is ', tissue, ' number of genes: ', len(tissue_genes))
            
            return tissue_genes
        else:
            return 'all'

    elif ans == '0':
        dt = pd.read_csv(
            './addData/Hum_Atas_normal_tissue.tsv',
            sep='\t')
        print('Gene universe is...')

        labels = sorted(list(map(str, set(dt['Tissue']))))
        index = [i for i in range(0, len(labels))]
        for i in index:
            print(index[i], '-----', labels[i])

        print("give me int, or input 'all'", sep='\n')
        tissue_id = str(input('Your choice: '))

        if tissue_id != 'all':
            tissue = labels[int(tissue_id)]
            tissue_genes = set(dt[(dt['Tissue'] == tissue) & (dt['Level'] != 'Not detected')]['Gene name'])
            print('your tissue is ', tissue, ' number of genes: ', len(tissue_genes))
            return tissue_genes
        else:
            return 'all'
    else:
        return tissue_selector()


def get_LCC(G):
    """
    :param G: networkX Graph
    :return: the Largest Connected Component, as NetrworkX object
    """
    CC_G = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    return CC_G[0]


def inflection_finder(card_LCC, n_CC, sigma):
    """
    :param card_LCC: cardinality of the LCC
    :param n_CC: number of connected components in the Network
    :return: the index of the last key node, after the removal of which the network stops rapidly falling apart
    """

    y = gaussian_filter1d(card_LCC, sigma=sigma)
    dy = np.diff(y)  # first derivative
    idx_max_dy = np.argmax(dy)
    if card_LCC[idx_max_dy] > n_CC[idx_max_dy]:
        return inflection_finder(card_LCC, n_CC, sigma + 0.2)
    else:
        return idx_max_dy


def node_centrality_calc(G):
    """
    :param G: a Graph
    :return: sorted dict of node centrality
    """
    centrality_node = nx.betweenness_centrality(G)
    degree_centrality = nx.degree_centrality(G)
    for k, v in centrality_node.items():
        centrality_node[k] = centrality_node[k] + degree_centrality[k]
    centrality_node = {k: v for k, v in sorted(centrality_node.items(), key=lambda item: item[1], reverse=True)}
    return centrality_node


def edge_centrality_calc(G):
    """
    :param G: a Graph
    :return: dict of edge centrality
    """
    edge_centrality = nx.edge_betweenness_centrality(G)
    return edge_centrality


def key_nodes_extractor(G, node_centrality):
    """

    :param G: NetworkX graph
    :param node_centrality: sorted dict of node centrality
    :return: list of key nodes and dict of some graph characteristics
    """
    G_copy = nx.Graph(G)
    graph_char = {'card_LCC': [len(G_copy.nodes())],
                  'n_CC': [len(list(nx.connected_component_subgraphs(G_copy)))],
                  'transitivity': [nx.transitivity(G_copy)],
                  'sh_path': [nx.average_shortest_path_length(G_copy)/len(G_copy.nodes())]}

    for k, v in node_centrality.items():

        if v == 0:
            break
        G_copy.remove_node(k)
        if len(G_copy.nodes()) == 0:
            break
        LCC_curent = get_LCC(G_copy)
        graph_char['card_LCC'].append(len(LCC_curent.nodes()))
        graph_char['n_CC'].append(len(list(nx.connected_component_subgraphs(G_copy))))
        graph_char['transitivity'].append(nx.transitivity(LCC_curent))
        graph_char['sh_path'].append(nx.average_shortest_path_length(LCC_curent)/len(LCC_curent.nodes()))


    # find inflection point of function

    idx_max_dy = inflection_finder(card_LCC=graph_char['card_LCC'],
                                   n_CC=graph_char['n_CC'],
                                   sigma=0.0001)

    key_nodes = dict()

    # idx_max_dy = np.argmax(np.array(graph_char['n_CC']) > np.array(graph_char['card_LCC']))

    graph_char['cutoff_point'] = idx_max_dy

    for i in range(0, idx_max_dy + 1):
        nods = list(node_centrality.keys())[i]
        key_nodes[nods] = node_centrality[nods]

    return key_nodes, graph_char


class Targets:
    """
    The class extracts and stores the targets of one microRNA and its name
    """

    miR_dict = {}
    with open('./baseData/hsa_miRTarBase.csv') as interact:  # import targets from miRTarBase
        for line in interact:
            (key, val) = line.strip().split(';')
            if key in miR_dict:
                miR_dict[key] += [val]
            else:
                miR_dict[key] = [val]

    def __init__(self, miR_name):
        self.miR_name = miR_name
        self.miR_targets = target_list_creator(miR_name=self.miR_name, miR_names=self.miR_dict.keys(),
                                               miR_dict=self.miR_dict)


class MainMirNetwork:

    def __init__(self, G, targets, CADgens, gene_set='all'):
        """

        :param G: NetworkX graph
        :param targets: targets of miRNA
        :param CADgens: disease gene
        :param gene_set: tissue or another gene set
        """
        self.CADgenes = CADgens
        self.tis_gene = gene_set
        self.LCC_G = get_LCC(G)
        if self.tis_gene != 'all':  # collapse an interactome to the tissue-specific interactome
            self.LCC_G = self.LCC_G.subgraph(self.tis_gene)
        self.LCC_miR_G = nx.Graph(get_LCC(self.LCC_G.subgraph(targets)))
        self.centrality_node = node_centrality_calc(G=self.LCC_miR_G)  # sorted dict
        self.centrality_edge = edge_centrality_calc(G=self.LCC_miR_G)

        for nods, dct in self.LCC_miR_G.nodes(data=True):
            if nods in CADgens:
                dct['CADgene'] = 1
            else:
                dct['CADgene'] = 0
            dct['BtweenCentrl'] = self.centrality_node[nods]

        for nods1, nods2, dct in self.LCC_miR_G.edges(data=True):
            dct['BtweenCentrl'] = self.centrality_edge[(nods1, nods2)]

        self.key_nodes, self.graph_char = \
            key_nodes_extractor(G=self.LCC_miR_G, node_centrality=self.centrality_node)

        for nods, dct in self.LCC_miR_G.nodes(data=True):
            if nods in self.key_nodes.keys():
                dct['key_nodes_flag'] = 1
            else:
                dct['key_nodes_flag'] = 0


#   visualisation tools

def draw_graph_to_cytoscape(miR_G, centrality_node):
    miR_G = nx.Graph(miR_G)  # unfreezing of the graph
    rem_CC = 0

    for CC in list(nx.connected_components(miR_G)):
        if len(CC) < 3:
            miR_G.remove_nodes_from(CC)
            rem_CC += 1

    print(rem_CC, ' network components with less than two nodes have been removed', end='\n')

    PORT_NUMBER = 1234
    IP = 'localhost'
    BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'

    requests.delete(BASE + 'session')  # Delete all networks in current session

    cytoscape_network = cy.from_networkx(miR_G)
    cytoscape_network['data']['name'] = 'miR_Net'
    res1 = requests.post(BASE + 'networks', data=json.dumps(cytoscape_network))
    res1_dict = res1.json()
    new_suid = res1_dict['networkSUID']
    requests.get(BASE + 'apply/layouts/force-directed/' + str(new_suid))

    # load and apply style

    res = requests.get(BASE + 'styles/miR_Net_Styles')
    if res.status_code != 200:

        with open('./options/cytoscape_styles/miR_Net_Styles.json') as json_file:
            miR_Net_Styles = json.load(json_file)

        for mapings in range(0, len(miR_Net_Styles['mappings'])):
            if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_LABEL_FONT_SIZE':
                miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())
            if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_SIZE':
                miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())
            if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_FILL_COLOR':
                miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())

        # Create new Visual Style
        res = requests.post(BASE + "styles", data=json.dumps(miR_Net_Styles))

    # Apply it to current network

    requests.get(BASE + 'apply/styles/' + 'miR_Net_Styles' + '/' + str(new_suid))  # !Это говно почему-то не работает


def draw_central_distr(miR_G, key_nodes, mir_name):
    """

    :param miR_G: a Graph of miRNA
    :param key_nodes: key_nodes from MainMiRNetwork.key_nodes
    :param mir_name: miR name for saving file
    :return: visualisation hist of centrality distribution
    """

    fig = plt.figure()
    ax = fig.add_subplot()
    N, bins, patches = ax.hist([miR_G.nodes[node]['BtweenCentrl'] for node in miR_G.nodes])
    ax.set(xlabel='Centrality',
           ylabel='Count of nodes')
    ax.yaxis.label.set_size(30)
    ax.xaxis.label.set_size(30)
    ax.tick_params(labelsize=20)

    ax.axvline(key_nodes[list(key_nodes.keys())[-1]], color='red', lw='4')
    
    right_side = ax.spines["right"]
    right_side.set_visible(False)
    top_side = ax.spines["top"]
    top_side.set_visible(False)

    # create gradient (grey_to_red hist path
    grey = Color('#cccccc')
    colors = list(grey.range_to(Color("red"), len(bins) - 1))
    for tmp_color, tmp_patch in zip(colors, patches):
        color = str(tmp_color)
        if len(color) < 7 and color[0] == '#':
            color = color + (7 - len(color)) * color[len(color) - 1]
        tmp_patch.set_facecolor(color)

    fig.set_figwidth(8.5)
    fig.set_figheight(8.5)

    plt.tight_layout()

    plt.savefig('./result/' + mir_name + '_centrality_distr.png', dpi=250)


def draw_key_nodes_extractor(card_LCC, n_CC, idx_max_dy, mir_name):
    """

    :param card_LCC:
    :param n_CC:
    :param idx_max_dy:
    :param mir_name:
    :return: visualisation plot of key nodes selection
    """

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(card_LCC, linewidth=4, label='Cardinality of the LCC')
    ax.plot(n_CC, linewidth=4, label='Count of CC', color='tab:green', linestyle='dashed')
    # ax.plot(idx_max_dy, card_LCC[idx_max_dy], marker='o', markersize=20, color="red")
    ax.axvline(idx_max_dy, color='red', lw='4')
    ax.minorticks_on()
    ax.grid(which='major',
            color='w',
            linewidth=1.3)
    ax.grid(which='minor',
            color='w',
            linestyle=':')
    ax.set(xlabel='Number of top nodes removed',
           ylabel='LCC cardinality / Count of CC')
    ax.legend()
    
    right_side = ax.spines["right"]
    right_side.set_visible(False)
    top_side = ax.spines["top"]
    top_side.set_visible(False)
    
#    plt.rc('font', size=10)  # controls default text sizes
    plt.rc('axes', labelsize=30)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=20)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=20)  # fontsize of the tick labels
    plt.rc('legend', fontsize=20)  # legend fontsize
    fig.set_figwidth(8)
    fig.set_figheight(8)
    plt.tight_layout()

    plt.savefig('./result/' + mir_name + 'key_nodes_selection.png', dpi=300)
