import mdtraj as md
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def build_neighbors(top):
    neigh = defaultdict(list)
    for a1, a2 in top.bonds:
        neigh[a1.index].append(a2.index)
        neigh[a2.index].append(a1.index)
    return neigh

def dfs_CA(start_CA, neigh, backbone_CA, backbone_N, backbone_C):
    all_chains = []

    def dfs(current_CA, used, chain):
        chain.append(current_CA)
        used.add(current_CA)

        N_candidates = [i for i in neigh[current_CA] if i in backbone_N and i not in used]

        if not N_candidates:
            all_chains.append(chain.copy())
            chain.pop()
            used.remove(current_CA)
            return

        for N in N_candidates:
            used.add(N)
            chain.append(N)

            C_candidates = [i for i in neigh[N] if i in backbone_C and i not in used]

            if not C_candidates:
                all_chains.append(chain.copy())
                chain.pop()
                used.remove(N)
                continue

            for C in C_candidates:
                used.add(C)
                chain.append(C)

                next_CA_candidates = [i for i in neigh[C] if i in backbone_CA and i not in used]

                if not next_CA_candidates:
                    all_chains.append(chain.copy())
                else:
                    for next_CA in next_CA_candidates:
                        dfs(next_CA, used, chain)

                chain.pop()
                used.remove(C)

            chain.pop()
            used.remove(N)

        chain.pop()
        used.remove(current_CA)

    dfs(start_CA, set(), [])
    return all_chains


def find_backbone_atoms(top, protect=True):
    neigh = build_neighbors(top)
    atoms = list(top.atoms)

    backbone_N = []
    backbone_CA = []
    backbone_C = []

    for atom in atoms:
        symbol = atom.element.symbol
        nb = [top.atom(i) for i in neigh[atom.index]]

        if symbol == "N":
            c_neighbors = [a for a in nb if a.element.symbol == "C"]
            if len(c_neighbors) > 1:
                backbone_N.append(atom.index)

        elif symbol == "C":
            o_neighbors = [a for a in nb if a.element.symbol == "O"]
            c_neighbors = [a for a in nb if a.element.symbol == "C"]

            if len(o_neighbors) >= 1 and any(a.element.symbol == "C" for a in nb):
                backbone_C.append(atom.index)
                continue

            #if any(a.element.symbol == "N" for a in nb) and any(a.element.symbol == "C" for a in nb):
            if any(a.element.symbol == "N" for a in nb):
                backbone_CA.append(atom.index)
    
    if protect:
        #to find [CA-N- C-CA-N- C-CA-N -C-CA-N- C] from N-terminus (protected by acetyl on N and methyl amine on C)
        
        end_CAs = []
        for ca in backbone_CA:
            connected_Cs = [i for i in neigh[ca] if i in backbone_C]
            if not connected_Cs:
                end_CAs.append(ca)
        '''
        chains = []
        for end_CA in end_CAs:
            chain = []
            used = set()
            current_CA = end_CA
            while True:
                chain.append(current_CA)

                N_candidates = [i for i in neigh[current_CA] if i in backbone_N and i not in used]
                if not N_candidates:
                    break
                N = N_candidates[0]
                chain.append(N)

                C_candidates = [i for i in neigh[N] if i in backbone_C and i not in used]
                if not C_candidates:
                    break
                C = C_candidates[0]

                used.update([current_CA, N, C])

                CA_candidates = [i for i in neigh[C] if i in backbone_CA and i not in used]
                if not CA_candidates:
                    break
                current_CA = CA_candidates[0]

            chains.append(chain)
        
        main_chain = [c for c in chains if len(c) == 12][0]
        '''

        all_chains = []
        for end_CA in end_CAs:
            chains_from_start = dfs_CA(end_CA, neigh, backbone_CA, backbone_N, backbone_C)
            all_chains.extend(chains_from_start)

        main_chain = [c for c in all_chains if len(c) == 12][0]

        # backbone is [CA-N- C-CA-N- C-CA-N -C-CA-N- C]
        return main_chain

import numpy as np

def make_phi_psi_indices(backbone_list):
    # backbon_list and phi_psi angle: [CA-N- C-CA-Ni+1- Ci-CAi-Ni -Ci-1-CA-N- C]

    phi_idx = []
    psi_idx = []

    N_i = backbone_list[7]
    CA_i = backbone_list[6]
    C_i = backbone_list[5]

    C_prev = backbone_list[8]
    N_next = backbone_list[4]

    phi_idx.append([C_prev, N_i, CA_i, C_i])
    psi_idx.append([N_i, CA_i, C_i, N_next])

    return np.array(phi_idx, dtype=int), np.array(psi_idx, dtype=int)

def phi_psi_calc(*args, out):
    if len(args) == 1:
        traj_file = args[0]
        m = md.load(traj_file)
    elif len(args) == 2:
        traj_file = args[0]
        top_file = args[1]
        m = md.load(traj_file, top=top_file)

    top = m.topology
    backbone_list = find_backbone_atoms(top)
    phi_idx, psi_idx = make_phi_psi_indices(backbone_list)

    phi_vals = np.degrees(md.compute_dihedrals(m, phi_idx))  
    psi_vals = np.degrees(md.compute_dihedrals(m, psi_idx))

    n_frames, n_residues_minus1 = phi_vals.shape
    with open(out, 'w') as f:
        f.write('phi,psi\n')
        for frame in range(n_frames):
            line = []
            for res in range(n_residues_minus1):
                line.append(f"{phi_vals[frame, res]}")
                line.append(f"{psi_vals[frame, res]}")
            f.write(",".join(line) + "\n")


def rama_plot(uaa_name, out):

    df = pd.read_csv(out)

    f = plt.figure(figsize=(4,4))
    ax1 = sns.kdeplot(x='phi', y='psi',data=df)
    ax1.set_title(f'Ramachandran plot of {uaa_name}')
    ax1.set_xlabel('Phi (degree)')
    ax1.set_ylabel('Psi (degree)')
    ax1.set_xlim(-180,180)
    ax1.set_ylim(-180,180)
    ax1.axhline(y=0, linewidth=1,color='gray',ls=':')
    ax1.axvline(x=0, linewidth=1,color='gray',ls=':')

    plt.xticks([-180,-135,-90,-45,0,45,90,135,180])
    plt.yticks([-180,-135,-90,-45,0,45,90,135,180])

    f.tight_layout()

    plt.savefig(out[:-4] + '.svg')
    plt.savefig(out[:-4] + '.png')