
import tkinter as tk

class Filter:
    
    filter_parameters = None
    def __init__(self):
        self.window_parameters = {
                'Bartlett': {
                    'largura_transicao_normalizada': 2.3,
                    'ondulacao_banda_passante_db': 0.185,
                    'lobulo_principal_lateral_db': 25,
                    'atenuacao_banda_rejeicao_db': 25,
                    'expressao': "w(n) = 2 - 2n/M"
                },
                'Hanning': {
                    'largura_transicao_normalizada': 3.1,
                    'ondulacao_banda_passante_db': 0.0546,
                    'lobulo_principal_lateral_db': 31,
                    'atenuacao_banda_rejeicao_db': 44,
                    'expressao': "w(n) = 0.5 + 0.5*cos(2πn/N)"
                },
                'Hamming': {
                    'largura_transicao_normalizada': 3.3,
                    'ondulacao_banda_passante_db': 0.0194,
                    'lobulo_principal_lateral_db': 41,
                    'atenuacao_banda_rejeicao_db': 53,
                    'expressao': "w(n) = 0.54 + 0.46*cos(2πn/N)"
                },
                'Blackman': {
                    'largura_transicao_normalizada': 5.5,
                    'ondulacao_banda_passante_db': 0.0017,
                    'lobulo_principal_lateral_db': 57,
                    'atenuacao_banda_rejeicao_db': 75,
                    'expressao': "w(n) = 0.42 + 0.5*cos(2πn/N) + 0.08*cos(4πn/(N-1))"
                },
                'Kaiser_beta_4.54': {
                    'largura_transicao_normalizada': 2.93,
                    'ondulacao_banda_passante_db': 0.0274,
                    'lobulo_principal_lateral_db': None,
                    'atenuacao_banda_rejeicao_db': 50,
                    'beta': 4.54,
                    'expressao': "I₀(β√(1-(2n/(N-1))²))/I₀(β)"
                }
            }
        
        self.filters = ["Passa-Baixa", "Passa-Alta", "Passa-Banda", "Rejeita-Banda"]
        
        # Variáveis de controle
        self.filter_type_var = tk.StringVar(value="Passa-Baixa")
        self.fs_var = tk.StringVar(value="8000")  # Frequência de amostragem
        self.fpassband1_var = tk.StringVar(value="1500")  # Frequência da borda da banda passante
        self.fpassband2_var = tk.StringVar(value="2500")  # Frequência da borda da banda passante
        self.transition_width_var = tk.StringVar(value="400")  # Largura de transição
        self.stopband_atten_var = tk.StringVar(value="50")  # Atenuação na banda de rejeição

        # Resultados
        self.filter_coeffs = None
        self.calculated_order = None
        self.window_function = None
        self.ideal_response = None
        
        # Variável para janela selecionada
        self.selected_window_var = tk.StringVar(value="Blackman")

    @staticmethod
    def get_instance():
        if Filter.filter_parameters is None:
            Filter.filter_parameters = Filter()
        return Filter.filter_parameters