#Читаем пары генов
with open("human_X_unimog.txt", "r") as f:
    human_ordered_long = [line.strip() for line in f if line.strip()]
with open("mustela_X_unimog.txt", "r") as f:
    mustela_ordered_long = [line.strip() for line in f if line.strip()]

#Создаём словарь: длинное имя гена -> целое число
id_to_int = {}
for i, (h_long, m_long) in enumerate(zip(human_ordered_long, mustela_ordered_long), start=1):
    id_to_int[h_long] = i
    id_to_int[m_long] = i

print(f"Создано сопоставление для {len(id_to_int)} уникальных генов.")

#Функция для извлечения порядка генов из FASTA файла
def get_chromosomal_order(fasta_file, id_map):
    order = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                long_id = line.split()[0].lstrip(">")
                if long_id in id_map:
                    order.append(id_map[long_id])
    return order

#Получаем порядок для каждого вида
human_order = get_chromosomal_order("homo_sapiens.faa", id_to_int)
mustela_order = get_chromosomal_order("mustela_nivalis.faa", id_to_int)

print(f"Генов в порядке хромосомы человека: {len(human_order)}")
print(f"Генов в порядке хромосомы ласки: {len(mustela_order)}")

# Проверяем, инвертирован ли порядок ласки
# Если первый ген ласки соответствует последнему гену человека — разворачиваем
if mustela_order[0] > human_order[0] and mustela_order[-1] < human_order[-1]:
    print("Обнаружена инверсия порядка")
    mustela_order = mustela_order[::-1]

#Записываем в формат UniMoG
with open("unimog_input.txt", "w") as f:
    f.write(">Human\n")
    f.write(" ".join(map(str, human_order)) + " |\n")
    f.write(">Mustela\n")
    f.write(" ".join(map(str, mustela_order)) + " |\n")

