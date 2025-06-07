"""
Software de Projeto de Filtros FIR
Desenvolvido seguindo a metodologia do Problema 03
"""

import tkinter as tk
import matplotlib

from FilterDesignApp import FilterDesignApp
matplotlib.use("TkAgg")
    

def main():
    """Função principal"""
    root = tk.Tk()
    
    app = FilterDesignApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()