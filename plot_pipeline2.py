import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

def parse_unimog(filename):
    genomes = {}
    current_name = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_name = line[1:]
                genomes[current_name] = []
            elif current_name is not None:
                parts = [g for g in line.split() if g not in ['|', '']]
                genomes[current_name].extend(parts)
    return genomes

def plot_synteny(input_file, output_file):
    print(f"Читаю {input_file}")
    genomes = parse_unimog(input_file)
    names = list(genomes.keys())
    
    g1 = genomes[names[0]]
    g2 = genomes[names[1]]
    
    def parse_genes(gene_list):
        nums = []
        signs = []
        for g in gene_list:
            if g.startswith('-'):
                nums.append(int(g[1:]))
                signs.append(-1)
            else:
                nums.append(int(g))
                signs.append(1)
        return nums, signs
    
    nums1, signs1 = parse_genes(g1)
    nums2, signs2 = parse_genes(g2)
    
    pos1 = {g: i+1 for i, g in enumerate(nums1)}
    pos2 = {g: i+1 for i, g in enumerate(nums2)}
    
    common = set(nums1) & set(nums2)
    
    x = []
    y = []
    colors = []
    
    for gene in sorted(common):
        x.append(pos1[gene])
        y.append(pos2[gene])
        idx1 = nums1.index(gene)
        idx2 = nums2.index(gene)
        if signs1[idx1] == signs2[idx2]:
            colors.append('blue')
        else:
            colors.append('red')
        
    x = np.array(x)
    y = np.array(y)
    
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    
    ax.scatter(x, y, alpha=0.8, s=20, c=colors, edgecolors='white', linewidth=0.5)
    
    max_pos = max(len(nums1), len(nums2))
    ax.plot([0, max_pos], [0, max_pos], 'k--', alpha=0.4, linewidth=1.5)
    
    ax.set_title(
        f'Синтения X-хромосомы: {names[0]} vs {names[1]}\n'
        f'{len(common)} строгих 1:1 ортологов (DCJ = 31)',
        fontsize=12, fontweight='bold', pad=10
    )
    
    ax.set_xlabel(f'Позиция гена на X-хромосоме {names[0]}', fontsize=10)
    ax.set_ylabel(f'Позиция гена на X-хромосоме {names[1]}', fontsize=10)
    
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Синтенные блоки'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Инверсии'),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='upper left', framealpha=0.9)
    
    ax.grid(True, alpha=0.3, linestyle=':')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"График сохранен: {output_file}")
    plt.close()

plot_synteny('unimog_input_correct.txt', 'synteny_final.png')
