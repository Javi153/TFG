import ast

best_score = 0
best_table = []
worst_score = 100
worst_table = []
score_mean = 0
tables = []
for i in range(43):
    archivo = "resultados_it_%d.txt" % i
    with open(archivo, "r") as file:
        aux = 0
        score_aux = 0
        while True:
            _ = file.readline()
            line1 = file.readline()[:-1]
            line2 = file.readline()[:-1]
            if not line1:
                break
            line1 = ast.literal_eval(line1)
            if line1 not in tables:
                aux += 1
                tables.append(line1)
                line2 = float(line2)
                score_aux += line2
                if line2 > 27:
                    print(line1, line2)
                if line2 > best_score:
                    best_score = line2
                    best_table = line1
                if line2 < worst_score:
                    worst_score = line2
                    worst_table = line1
    score_mean += score_aux / aux
score_mean /= 43
print("La mejor puntuación es", best_score, best_table)
print("La peor puntuación es", worst_score, worst_table)
print("La media de valores es", score_mean)
print(len(tables))
