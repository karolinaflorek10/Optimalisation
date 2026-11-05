import csv
import math
import matplotlib.pyplot as plt

x1_hj, x2_hj = [], []
x1_r, x2_r = [], []

with open('arkusz_wykres.csv', mode='r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter=';')
    next(reader)  # pominięcie nagłówka
    for row in reader:
        if len(row) >= 3 and row[1] != "" and row[2] != "":
            x1_hj.append(float(row[1].replace(',', '.')))
            x2_hj.append(float(row[2].replace(',', '.')))
        if len(row) >= 5 and row[3] != "" and row[4] != "":
            x1_r.append(float(row[3].replace(',', '.')))
            x2_r.append(float(row[4].replace(',', '.')))

# Siatka do funkcji celu
N = 350
x1_vals = [-1.0 + 2.0 * i / (N - 1) for i in range(N)]
x2_vals = [-1.0 + 2.0 * i / (N - 1) for i in range(N)]

Z = []
for x2 in x2_vals:
    row_z = []
    for x1 in x1_vals:
        val = x1**2 + x2**2 - math.cos(2.5 * math.pi * x1) - math.cos(2.5 * math.pi * x2) + 2.0
        row_z.append(val)
    Z.append(row_z)

plt.figure(figsize=(8.5, 7.5))
cs = plt.contour(x1_vals, x2_vals, Z, levels=25, cmap='viridis')
plt.clabel(cs, inline=1, fontsize=8)

# Ścieżka Hooke'a-Jeevesa
if x1_hj:
    plt.plot(x1_hj, x2_hj, 'r.-', linewidth=2.0, markersize=7, label="Ścieżka Hooke'a-Jeevesa")
    plt.plot(x1_hj[0], x2_hj[0], 'go', markersize=9, label='Start $x^{(0)}$')
    plt.plot(x1_hj[-1], x2_hj[-1], 'b*', markersize=13, label='Optimum $x^*_{HJ}$')

# Opcjonalnie: ścieżka Rosenbrocka (przerywana pomarańczowa)
if x1_r:
    plt.plot(x1_r, x2_r, color='#e36c09', linestyle='--', marker='o', markersize=4, alpha=0.7, label='Ścieżka Rosenbrocka')
    plt.plot(x1_r[-1], x2_r[-1], color='#e36c09', marker='X', markersize=10, label='Optimum $x^*_{Rosen}$')

plt.title("Poziomice funkcji testowej z trajektorią optymalizacji", fontsize=12)
plt.xlabel('$x_1$', fontsize=11)
plt.ylabel('$x_2$', fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.savefig('wykres_poziomic_hj.png', dpi=300)
print("Wygenerowano zaktualizowany wykres zgodny z arkuszem: wykres_poziomic_hj.png")
