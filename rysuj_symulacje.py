import csv
import matplotlib.pyplot as plt

t_vals, alpha_vals, omega_vals = [], [], []

with open('symulacja.csv', mode='r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter=';')
    next(reader) # pominięcie nagłówka
    for row in reader:
        if len(row) >= 3:
            t_vals.append(float(row[0].replace(',', '.')))
            alpha_vals.append(float(row[1].replace(',', '.')))
            omega_vals.append(float(row[2].replace(',', '.')))

# Wykres położenia
plt.figure(figsize=(8, 4.5))
plt.plot(t_vals, alpha_vals, 'b-', linewidth=2, label=r'$\alpha(t)$ - położenie kątowe')
plt.axhline(3.14, color='r', linestyle='--', label=r'$\alpha_{ref} = 3.14$ rad')
plt.title('Przebieg położenia kątowego ramienia robota w czasie')
plt.xlabel('Czas t [s]')
plt.ylabel('Kąt obrotu [rad]')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig('wykres_polozenie.png', dpi=300)

# Wykres prędkości
plt.figure(figsize=(8, 4.5))
plt.plot(t_vals, omega_vals, 'g-', linewidth=2, label=r'$\omega(t)$ - prędkość kątowa')
plt.axhline(0.0, color='r', linestyle='--', label=r'$\omega_{ref} = 0.0$ rad/s')
plt.title('Przebieg prędkości kątowej ramienia robota w czasie')
plt.xlabel('Czas t [s]')
plt.ylabel('Prędkość kątowa [rad/s]')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig('wykres_predkosc.png', dpi=300)

print("Wygenerowano wykresy: wykres_polozenie.png oraz wykres_predkosc.png")
