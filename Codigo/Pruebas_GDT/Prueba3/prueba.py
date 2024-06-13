import ast

best_score = 0
best_table = []
for i in range(20):
    archivo = "resultados_it_%d.txt" % i
    with open(archivo, "r") as file:
        while True:
            _ = file.readline()
            line1 = file.readline()[:-1]
            line2 = file.readline()[:-1]
            if not line1:
                break
            line1 = ast.literal_eval(line1)
            line2 = float(line2)
            if line2 > best_score:
                best_score = line2
                best_table = line1
print(best_score, best_table)