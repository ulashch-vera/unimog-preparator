import matplotlib.pyplot as plt
import numpy as np

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
                parts = [g.replace('|', '') for g in line.split() if g not in ['|', '']]
                genomes[current_name].extend(parts)
    return genomes

def plot_synteny_clean(input_file, output_file, reverse_genome2=False):
    print(f"Обрабатываю {input_file}")
    genomes = parse_unimog(input_file)
    names = list(genomes.keys())
    
    g1 = genomes[names[0]]
    g2 = genomes[names[1]]
    
    if reverse_genome2:
        g2 = g2[::-1]
        print("Порядок второго генома развернут")

    #Создаем числовые перестановки
    all_genes = list(set(g1 + g2))
    gene_to_num = {gene: i+1 for i, gene in enumerate(all_genes)}
    
    perm1 = [gene_to_num[g] for g in g1 if g in gene_to_num]
    perm2 = [gene_to_num[g] for g in g2 if g in gene_to_num]

    n_genes = len(perm1)
    print(f"  Найдено {n_genes} общих генов.")

    x = np.array(range(1, n_genes + 1))
    pos_in_g2 = {g: i+1 for i, g in enumerate(perm2)}
    y = np.array([pos_in_g2.get(g, 0) for g in perm1])
    
    mask = y > 0
    x, y = x[mask], y[mask]

    slope, intercept = np.polyfit(x, y, 1)
    trend_line = slope * x + intercept

    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    
    ax.scatter(x, y, alpha=0.6, s=15, color='blue', label='Orthologous genes')
    ax.plot(x, trend_line, color='red', linestyle='--', linewidth=1.5, label='Linear trend')

    ax.set_title(
        f'Синтения X-хромосомы: Homo sapiens vs Mustela\n'
        f'{n_genes} строгих 1:1 ортологов',
        fontsize=14,
        fontweight='bold',
        pad=15
    )
    
    ax.set_xlabel('Позиция гена на X-хромосоме человека', fontsize=12)
    ax.set_ylabel('Позиция гомологичного гена на X-хромосоме Mustela', fontsize=12)

    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle=':')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Сохранено как: {output_file}\n")
    plt.close()

plot_synteny_clean('unimog_input.txt', 'synteny_pipeline2_correct.png', reverse_genome2=False)

plot_synteny_clean('unimog_input.txt', 'synteny_pipeline2_wrong.png', reverse_genome2=True)

