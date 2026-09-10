import pandas as pd
import numpy as np
import os

N=20

COL_NAMES = ['ALA', 'VAL', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'THR', 'GLN', 'CYS', 'ASN', 'SER', 'TYR']

def array_to_latex(arr: np.array[int], score: float):
    df = pd.DataFrame(arr, COL_NAMES, COL_NAMES)
    latex_code = df.to_latex(float_format="%.0f", caption=f"Similarity score {score}\\%", column_format="|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|")
    latex_code = latex_code.replace("\\\\\n", "\\\\\n\\hline\n")
    latex_code = latex_code.replace("\\hline\n\\hline", "\\hline")
    latex_code = latex_code.replace("\\toprule\n", "\\hline\n")
    latex_code = latex_code.replace("\\bottomrule\n", "")
    latex_code = latex_code.replace("\\midrule\n", "")
    return latex_code

def matrix_to_np(matrix: list[int])->np.array[int]:
    arr = np.zeros((N, N))
    for i in range(N):
        arr[i,i:] = matrix[i]
        arr[i:,i] = matrix[i]
    return arr

def parse_matrix_over_score(file: str, score:float=28.6) -> np.array:
    f = open(file, "r")
    lines = f.read().splitlines()
    f.close()
    scores = [eval(x) for x in lines[2::3]]
    matrixes = [eval(lines[3*i+1]) for i in range(len(lines)//3) if scores[i] >= score]
    scores = [x for x in scores if x > score]
    matrixes = [array_to_latex(matrix_to_np(matrixes[i]), scores[i]) for i in range(len(matrixes))]
    return matrixes, scores

def main(out_file:str="out.txt"):
    files = os.listdir(".")
    files = [x for x in files if ".txt" in x and x != "final_result.txt"]
    f = open(out_file, "w")
    for file in files:
        matrixes, scores = parse_matrix_over_score(file)
        for i in range(len(matrixes)):
            f.write(matrixes[i])
            f.write("\n")
            f.write(str(scores[i]))
            f.write("\n")
    f.close()

if __name__ == "__main__":
    os.chdir("./Pruebas_GDT/Prueba6/results")
    main()