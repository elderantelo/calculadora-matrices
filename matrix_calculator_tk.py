
import tkinter as tk
from tkinter import ttk, messagebox
import random

# -------------------- Matrix Utilities (pure Python) --------------------
def create_matrix(r, c, fill=0.0):
    return [[float(fill) for _ in range(c)] for __ in range(r)]

def clone_matrix(A):
    return [row[:] for row in A]

def resize_matrix(A, r, c):
    R = create_matrix(r, c, 0.0)
    for i in range(min(r, len(A))):
        for j in range(min(c, len(A[0]))):
            R[i][j] = A[i][j]
    return R

def add(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("Las dimensiones deben coincidir para A + B / A - B.")
    r, c = len(A), len(A[0])
    R = create_matrix(r, c)
    for i in range(r):
        for j in range(c):
            R[i][j] = A[i][j] + B[i][j]
    return R

def sub(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("Las dimensiones deben coincidir para A - B / B - A.")
    r, c = len(A), len(A[0])
    R = create_matrix(r, c)
    for i in range(r):
        for j in range(c):
            R[i][j] = A[i][j] - B[i][j]
    return R

def mul(A, B):
    if len(A[0]) != len(B):
        raise ValueError("Columnas de A deben igualar filas de B para A × B.")
    r, c, n = len(A), len(B[0]), len(A[0])
    R = create_matrix(r, c)
    for i in range(r):
        for j in range(c):
            s = 0.0
            for k in range(n):
                s += A[i][k] * B[k][j]
            R[i][j] = s
    return R

def transpose(A):
    r, c = len(A), len(A[0])
    R = create_matrix(c, r)
    for i in range(r):
        for j in range(c):
            R[j][i] = A[i][j]
    return R

def is_square(A):
    return len(A) == len(A[0])

def trace(A):
    if not is_square(A):
        raise ValueError("La traza requiere una matriz cuadrada.")
    return sum(A[i][i] for i in range(len(A)))

def determinant(Ain):
    n = len(Ain)
    if not is_square(Ain):
        raise ValueError("El determinante requiere una matriz cuadrada.")
    A = clone_matrix(Ain)
    det = 1.0
    for i in range(n):
        # Partial pivot
        max_row = i
        for r in range(i+1, n):
            if abs(A[r][i]) > abs(A[max_row][i]):
                max_row = r
        if abs(A[max_row][i]) < 1e-12:
            return 0.0
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            det *= -1.0
        det *= A[i][i]
        pivot = A[i][i]
        for r in range(i+1, n):
            factor = A[r][i] / pivot
            for c in range(i, n):
                A[r][c] -= factor * A[i][c]
    return det

def inverse(Ain):
    n = len(Ain)
    if not is_square(Ain):
        raise ValueError("La inversa requiere una matriz cuadrada.")
    A = clone_matrix(Ain)
    I = create_matrix(n, n, 0.0)
    for i in range(n):
        I[i][i] = 1.0
    for i in range(n):
        # Pivot
        max_row = i
        for r in range(i+1, n):
            if abs(A[r][i]) > abs(A[max_row][i]):
                max_row = r
        if abs(A[max_row][i]) < 1e-12:
            raise ValueError("La matriz es singular, no tiene inversa.")
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            I[i], I[max_row] = I[max_row], I[i]
        # Normalize row i
        pivot = A[i][i]
        for c in range(n):
            A[i][c] /= pivot
            I[i][c] /= pivot
        # Eliminate other rows
        for r in range(n):
            if r != i:
                factor = A[r][i]
                for c in range(n):
                    A[r][c] -= factor * A[i][c]
                    I[r][c] -= factor * I[i][c]
    return I

def rank(Ain):
    A = clone_matrix(Ain)
    rows, cols = len(A), len(A[0])
    row = 0
    rnk = 0
    for col in range(cols):
        if row >= rows:
            break
        pivot_row = row
        for i in range(row+1, rows):
            if abs(A[i][col]) > abs(A[pivot_row][col]):
                pivot_row = i
        if abs(A[pivot_row][col]) < 1e-12:
            continue
        # Swap
        A[row], A[pivot_row] = A[pivot_row], A[row]
        pivot = A[row][col]
        for j in range(col, cols):
            A[row][j] /= pivot
        for i in range(rows):
            if i != row:
                f = A[i][col]
                for j in range(col, cols):
                    A[i][j] -= f * A[row][j]
        row += 1
        rnk += 1
    return rnk

# -------------------- UI Widgets --------------------
class MatrixGrid(ttk.Frame):
    def __init__(self, master, title="Matriz", rows=3, cols=3, **kwargs):
        super().__init__(master, **kwargs)
        self.title = title
        self.rows = rows
        self.cols = cols

        self.header = ttk.Frame(self)
        self.header.pack(fill="x", pady=(0,4))

        self.title_label = ttk.Label(self.header, text=f"{self.title} ({self.rows}×{self.cols})", font=("Segoe UI", 10, "bold"))
        self.title_label.pack(side="left")

        self.btn_random = ttk.Button(self.header, text="Aleatorio", command=self.fill_random)
        self.btn_random.pack(side="right", padx=4)
        self.btn_clear = ttk.Button(self.header, text="Limpiar", command=self.clear)
        self.btn_clear.pack(side="right")

        sizebar = ttk.Frame(self)
        sizebar.pack(fill="x", pady=(0,6))
        ttk.Label(sizebar, text="Filas:").pack(side="left")
        self.spin_r = ttk.Spinbox(sizebar, from_=1, to=10, width=3, command=self.update_size)
        self.spin_r.set(self.rows)
        self.spin_r.pack(side="left", padx=(2,8))
        ttk.Label(sizebar, text="Columnas:").pack(side="left")
        self.spin_c = ttk.Spinbox(sizebar, from_=1, to=10, width=3, command=self.update_size)
        self.spin_c.set(self.cols)
        self.spin_c.pack(side="left", padx=(2,8))

        # Scrollable grid
        self.canvas = tk.Canvas(self, height=220)
        self.grid_frame = ttk.Frame(self.canvas)
        self.vscroll = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vscroll.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.vscroll.pack(side="right", fill="y")
        self.canvas.create_window((0,0), window=self.grid_frame, anchor="nw")
        self.grid_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

        self.entries = []
        self.build_entries()

    def update_size(self):
        try:
            r = int(self.spin_r.get())
            c = int(self.spin_c.get())
            r = max(1, min(10, r))
            c = max(1, min(10, c))
        except ValueError:
            return
        self.rows, self.cols = r, c
        self.title_label.config(text=f"{self.title} ({self.rows}×{self.cols})")
        self.build_entries()

    def build_entries(self):
        for w in self.grid_frame.winfo_children():
            w.destroy()
        self.entries = []
        for i in range(self.rows):
            row_entries = []
            for j in range(self.cols):
                e = ttk.Entry(self.grid_frame, width=10, justify="right")
                e.grid(row=i, column=j, padx=3, pady=3, sticky="nsew")
                e.insert(0, "0")
                row_entries.append(e)
            self.entries.append(row_entries)

    def get_matrix(self):
        A = create_matrix(self.rows, self.cols, 0.0)
        for i in range(self.rows):
            for j in range(self.cols):
                txt = self.entries[i][j].get().strip().replace(",", ".")
                try:
                    A[i][j] = float(txt) if txt != "" else 0.0
                except ValueError:
                    raise ValueError(f"Valor inválido en {self.title} [{i+1},{j+1}]: '{txt}'")
        return A

    def set_matrix(self, M):
        r, c = len(M), len(M[0])
        self.spin_r.set(r)
        self.spin_c.set(c)
        self.update_size()
        for i in range(r):
            for j in range(c):
                self.entries[i][j].delete(0, tk.END)
                self.entries[i][j].insert(0, f"{M[i][j]:.6g}")

    def clear(self):
        for i in range(self.rows):
            for j in range(self.cols):
                self.entries[i][j].delete(0, tk.END)
                self.entries[i][j].insert(0, "0")

    def fill_random(self):
        for i in range(self.rows):
            for j in range(self.cols):
                val = random.uniform(-5, 5)
                self.entries[i][j].delete(0, tk.END)
                self.entries[i][j].insert(0, f"{val:.2f}")

class ResultView(ttk.Frame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.scalar_var = tk.StringVar(value="")
        self.scalar_label = ttk.Label(self, textvariable=self.scalar_var, font=("Segoe UI", 11, "bold"))
        self.scalar_label.pack(anchor="w", pady=(0,4))

        self.canvas = tk.Canvas(self, height=240)
        self.grid_frame = ttk.Frame(self.canvas)
        self.vscroll = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vscroll.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.vscroll.pack(side="right", fill="y")
        self.canvas.create_window((0,0), window=self.grid_frame, anchor="nw")
        self.grid_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

    def show_scalar(self, value):
        self.scalar_var.set(f"Resultado (escalar): {value:.6g}")
        for w in self.grid_frame.winfo_children():
            w.destroy()

    def show_matrix(self, M):
        self.scalar_var.set("Resultado (matriz)")
        for w in self.grid_frame.winfo_children():
            w.destroy()
        r, c = len(M), len(M[0])
        for i in range(r):
            for j in range(c):
                lab = ttk.Label(self.grid_frame, text=f"{M[i][j]:.6g}", width=12, anchor="e", relief="solid")
                lab.grid(row=i, column=j, padx=2, pady=2, sticky="nsew")

# -------------------- Main App --------------------
class MatrixCalculatorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Calculadora de Matrices - ALGEBRA LINEAL UPDS")
        self.geometry("1100x720")
        self.minsize(900, 600)

        style = ttk.Style(self)
        try:
            self.call("source", "azure.tcl")  # If available
            style.theme_use("azure")
        except Exception:
            style.theme_use("clam")

        # Header
        header = ttk.Frame(self, padding=8)
        header.pack(fill="x")
        ttk.Label(header, text="Calculadora de Matrices - Algebra Lineal  Docente: Ing. LIGIA BASPINEIRO ZABALA  UPDS", font=("Segoe UI", 16, "bold")).pack(anchor="w")
        ttk.Label(header, text="Operaciones: suma, resta, multiplicación, transpuesta, traza, determinante, inversa y rango.").pack(anchor="w")

        # Grids
        grids = ttk.Frame(self, padding=8)
        grids.pack(fill="both", expand=True)
        self.gridA = MatrixGrid(grids, title="Matriz A", rows=3, cols=3)
        self.gridB = MatrixGrid(grids, title="Matriz B", rows=3, cols=3)
        self.gridA.pack(side="left", fill="both", expand=True, padx=(0,8))
        self.gridB.pack(side="left", fill="both", expand=True, padx=(8,0))

        # Sidebar operations
        ops = ttk.Frame(self, padding=8)
        ops.pack(fill="x")
        btns = [
            ("A + B", lambda: self._binary_op(add)),
            ("A − B", lambda: self._binary_op(sub)),
            ("B − A", lambda: self._binary_reverse_op(sub)),
            ("A × B", lambda: self._binary_op(mul)),
            ("Transponer A", lambda: self._unary_op_A(transpose)),
            ("Transponer B", lambda: self._unary_op_B(transpose)),
            ("Traza(A)", lambda: self._scalar_A(trace)),
            ("Traza(B)", lambda: self._scalar_B(trace)),
            ("Determinante(A)", lambda: self._scalar_A(determinant)),
            ("Determinante(B)", lambda: self._scalar_B(determinant)),
            ("Inversa(A)", lambda: self._unary_op_A(inverse)),
            ("Inversa(B)", lambda: self._unary_op_B(inverse)),
            ("Rango(A)", lambda: self._scalar_A(rank)),
            ("Rango(B)", lambda: self._scalar_B(rank)),
        ]

        rowframe = ttk.Frame(ops)
        rowframe.pack(fill="x")
        for i, (label, cmd) in enumerate(btns):
            b = ttk.Button(rowframe, text=label, command=cmd)
            b.grid(row=i//4, column=i%4, padx=6, pady=6, sticky="ew")
        for i in range(4):
            rowframe.grid_columnconfigure(i, weight=1)

        # Result View
        self.result = ResultView(self, padding=8)
        self.result.pack(fill="both", expand=True)

        # Status bar
        self.status = tk.StringVar(value="Listo.")
        statusbar = ttk.Label(self, textvariable=self.status, anchor="w", padding=4, relief="groove")
        statusbar.pack(fill="x")

    # Helper methods
    def _get_A(self):
        return self.gridA.get_matrix()
    def _get_B(self):
        return self.gridB.get_matrix()

    def _show_result(self, out):
        if isinstance(out, (int, float)):
            self.result.show_scalar(float(out))
        else:
            self.result.show_matrix(out)

    def _binary_op(self, fn):
        try:
            A = self._get_A()
            B = self._get_B()
            out = fn(A, B)
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

    def _binary_reverse_op(self, fn):
        try:
            A = self._get_A()
            B = self._get_B()
            out = fn(B, A)  # B - A
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

    def _unary_op_A(self, fn):
        try:
            A = self._get_A()
            out = fn(A)
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

    def _unary_op_B(self, fn):
        try:
            B = self._get_B()
            out = fn(B)
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

    def _scalar_A(self, fn):
        try:
            A = self._get_A()
            out = fn(A)
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

    def _scalar_B(self, fn):
        try:
            B = self._get_B()
            out = fn(B)
            self._show_result(out)
            self.status.set("Operación completada.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status.set("Error: " + str(e))

 # Status bar

if __name__ == "__main__":
    app = MatrixCalculatorApp()
    app.mainloop()

