import numpy as np
import matplotlib.pyplot as plt

x_data = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1])
y_data = np.array([-2.3026, -0.69315, -0.10536, 0.26236, 0.53063, 0.74194])


c1 = np.polyfit(x_data, y_data, 1)[::-1]   # [a0, a1]
c2 = np.polyfit(x_data, y_data, 2)[::-1]   # [a0, a1, a2]
c3 = np.polyfit(x_data, y_data, 3)[::-1]   # [a0, a1, a2, a3]

def poly(c, x):
    """c[0] + c[1]*x + c[2]*x^2 + ... — порядок как в C++."""
    y = np.zeros_like(x, dtype=float)
    xp = np.ones_like(x, dtype=float)
    for k in c:
        y += k * xp
        xp *= x
    return y

def fmt(c):
    terms = []
    for i, a in enumerate(c):
        sign = '+' if a >= 0 else '-'
        s = f"{abs(a):.4f}"
        if i > 0:
            s += "x" + (f"^{i}" if i > 1 else "")
        terms.append((sign, s))
    out = ("-" if terms[0][0] == '-' else "") + terms[0][1]
    for sign, s in terms[1:]:
        out += f" {sign} {s}"
    return out

def phi(c, x, y):
    return np.sum((poly(c, x) - y) ** 2)   # Φ из (3.15), как в методичке

x_plot = np.linspace(0.05, 2.15, 400)

plt.figure(figsize=(12, 6))
plt.plot(x_plot, poly(c1, x_plot), 'r-',
         label=f'P1: {fmt(c1)}  (e={phi(c1, x_data, y_data):.4f})')
plt.plot(x_plot, poly(c2, x_plot), 'g-',
         label=f'P2: {fmt(c2)}  (e={phi(c2, x_data, y_data):.4f})')
plt.plot(x_plot, poly(c3, x_plot), 'b-',
         label=f'P3: {fmt(c3)}  (e={phi(c3, x_data, y_data):.4f})')
plt.plot(x_data, y_data, 'ko', label='Исходные точки', zorder=5)

plt.grid(True, alpha=0.3)
plt.legend(fontsize=9)
plt.title('Аппроксимация полиномами 1-й, 2-й и 3-й степени (МНК)')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()

best = int(np.argmin([phi(c1, x_data, y_data),
                      phi(c2, x_data, y_data),
                      phi(c3, x_data, y_data)])) + 1
print(f"Лучшее приближение: F{best}(x)")