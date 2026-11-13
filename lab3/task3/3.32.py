import numpy as np
import matplotlib.pyplot as plt

x_data = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1])
y_data = np.array([-2.3026, -0.69315, -0.10536, 0.26236, 0.53063, 0.74194])

# np.polyfit возвращает [a_n, ..., a_1, a_0] — старшие первыми.
# [::-1] -> [a0, a1, ..., an], как в C++.
c1 = np.polyfit(x_data, y_data, 1)[::-1]
c2 = np.polyfit(x_data, y_data, 2)[::-1]
c3 = np.polyfit(x_data, y_data, 3)[::-1]
c5 = np.polyfit(x_data, y_data, 5)[::-1]

def poly(c, x):
    """c[0] + c[1]*x + ... — порядок как в C++ (свободный член первым)."""
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
    return np.sum((poly(c, x) - y) ** 2)   # Φ из (3.15)

x_plot = np.linspace(0.05, 2.15, 400)

plt.figure(figsize=(13, 7))
plt.plot(x_plot, poly(c1, x_plot), 'r-',  lw=1.5,
         label=f'F1: {fmt(c1)}  (Φ={phi(c1, x_data, y_data):.4f})')
plt.plot(x_plot, poly(c2, x_plot), 'g-',  lw=1.5,
         label=f'F2: {fmt(c2)}  (Φ={phi(c2, x_data, y_data):.4f})')
plt.plot(x_plot, poly(c3, x_plot), 'b-',  lw=1.5,
         label=f'F3: {fmt(c3)}  (Φ={phi(c3, x_data, y_data):.4f})')
plt.plot(x_plot, poly(c5, x_plot), 'm--', lw=1.5,
         label=f'F5: {fmt(c5)}  (Φ={phi(c5, x_data, y_data):.2e})')
plt.plot(x_data, y_data, 'ko', ms=7, zorder=5, label='Исходные точки')

plt.grid(True, alpha=0.3)
plt.legend(fontsize=8, loc='best')
plt.title('Аппроксимация полиномами 1-й, 2-й, 3-й и 5-й степени (МНК)')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()