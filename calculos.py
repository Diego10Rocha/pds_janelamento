import numpy as np
from tkinter import  messagebox
from scipy import signal
from parametros import Filter



def design_filter():
    
    filter_parameters = Filter.get_instance()
    """Projeta o filtro com as especificações fornecidas"""
    try:
        # Obter especificações
        fs = float(filter_parameters.fs_var.get())
        fp = float(filter_parameters.fp_var.get())
        transition_width = float(filter_parameters.transition_width_var.get())
        stopband_atten = float(filter_parameters.stopband_atten_var.get())
        filter_type = filter_parameters.filter_type_var.get()
        window_name = filter_parameters.selected_window_var.get()
        
        # Validações
        if fp >= fs/2:
            raise ValueError("Frequência da banda passante deve ser menor que Fs/2")
        
        # Calcular frequência de corte baseada no tipo de filtro
        if filter_type == "Passa-Baixa":
            fs_freq = fp + transition_width  # Frequência de stopband
            fc1 = (fp + fs_freq) / 2  # Frequência de corte centrada na transição
        elif filter_type == "Passa-Alta":  # Passa-Alta
            fs_freq = fp - transition_width  # Frequência de stopband
            fc1 = (fs_freq + fp) / 2  # Frequência de corte centrada na transição
            
            if fs_freq <= 0:
                raise ValueError("Para passa-alta: fp - largura_transição deve ser > 0")
        elif filter_type == "Passa-Banda":
            f1 = fp / 2
            f2 = fp + transition_width / 2
            if f1 <= 0 or f2 >= fs / 2:
                raise ValueError("Frequências do passa-banda devem estar no intervalo válido.")
            fc1 = f1 / (fs / 2)
            fc2 = f2 / (fs / 2)

        elif filter_type == "Rejeita-Banda":
            f1 = fp / 2
            f2 = fp + transition_width / 2
            if f1 <= 0 or f2 >= fs / 2:
                raise ValueError("Frequências do rejeita-banda devem estar no intervalo válido.")
            fc1 = f1 / (fs / 2)
            fc2 = f2 / (fs / 2)
        else:
            raise ValueError("Tipo de filtro desconhecido.")
        
        # Calcular ordem do filtro
        delta_f_norm = transition_width / fs
        window_params = filter_parameters.window_parameters[window_name]
        
        # Usar fator da janela para calcular ordem
        if 'Kaiser' in window_name:
            # Para Kaiser, usar fórmula específica
            A = stopband_atten
            beta = window_params['beta']
            M = (A - 8) / (2.285 * 2 * np.pi * delta_f_norm)
            order = int(np.ceil(M))
        else:
            # Para outras janelas, usar fator de largura de transição
            factor = window_params['largura_transicao_normalizada']
            N_calc = factor / delta_f_norm
            order = int(np.ceil(N_calc))
        
        # Garantir ordem ímpar para simetria
        if order % 2 == 0:
            order += 1
        
        # Limitar ordem
        order = max(11, min(501, order))
        
        # Projetar filtro ideal
        nyquist = fs / 2
        #fc_norm = fc1 / nyquist
        
        if filter_type == "Passa-Baixa":
            h_ideal = ideal_lowpass(order, fc1)
        elif filter_type == "Passa-Alta":  # Passa-Alta
            h_ideal = ideal_highpass(order, fc1)
        elif filter_type == "Passa-Banda":
            h_ideal = ideal_bandpass(order, fc1, fc2)
        elif filter_type == "Rejeita-Banda":
            h_ideal = ideal_bandstop(order, fc1, fc2)
        else:
            raise ValueError("Tipo de filtro desconhecido.")
        
        # Criar janela
        if 'Kaiser' in window_name:
            beta = window_params['beta']
            window = signal.get_window(('kaiser', beta), order)
        elif window_name == 'Retangular':
            window = signal.get_window('boxcar', order)
        elif window_name == 'Bartlett':
            window = signal.get_window('bartlett', order)
        elif window_name == 'Hanning':
            window = signal.get_window('hann', order)
        elif window_name == 'Hamming':
            window = signal.get_window('hamming', order)
        elif window_name == 'Blackman':
            window = signal.get_window('blackman', order)
        
        # Aplicar janelamento
        h_windowed = h_ideal * window
        
        # Armazenar resultados
        filter_parameters.filter_coeffs = h_windowed
        filter_parameters.calculated_order = order
        filter_parameters.window_function = window
        filter_parameters.ideal_response = h_ideal
        
        # Atualizar visualizações
        #self.update_all_plots(h_windowed, window, h_ideal, fs, fc)
        
        # Mostrar cálculos detalhados
        #filter_parameters.show_calculations(fs, fp, transition_width, stopband_atten, 
        #                        filter_type, window_name, order, fc, 
        #                        h_windowed, delta_f_norm)
        
        #messagebox.showinfo("Sucesso", "Filtro projetado com sucesso!")
        return (h_windowed, window, h_ideal, fs, fc1, fc2)
        
    except ValueError as e:
        messagebox.showerror("Erro de Entrada", str(e))
    except Exception as e:
        messagebox.showerror("Erro", f"Erro no projeto: {e}")

def ideal_lowpass(N, fc_norm):
    """Calcula filtro passa-baixa ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)
    
    for i in range(N):
        n_centered = n[i] - M/2
        
        if abs(n_centered) < 1e-10:
            h[i] = 2 * fc_norm
        else:
            omega_c = fc_norm * np.pi
            h[i] = 2 * fc_norm * np.sin(omega_c * n_centered) / (omega_c * n_centered)
    
    return h

def ideal_highpass(N, fc_norm):
    """Calcula filtro passa-alta ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)
    
    for i in range(N):
        n_centered = n[i] - M/2
        
        if abs(n_centered) < 1e-10:
            h[i] = 1.0 - 2 * fc_norm
        else:
            # Impulso delta menos passa-baixa
            delta_impulse = np.sin(np.pi * n_centered) / (np.pi * n_centered)
            omega_c = fc_norm * np.pi
            lowpass_term = 2 * fc_norm * np.sin(omega_c * n_centered) / (omega_c * n_centered)
            h[i] = delta_impulse - lowpass_term
    
    return h

def ideal_bandpass(N, fc1, fc2):
    """Calcula filtro passa-banda ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)

    for i in range(N):
        n_centered = n[i] - M / 2
        if abs(n_centered) < 1e-10:
            h[i] = 2 * (fc2 - fc1)
        else:
            h[i] = (np.sin(2 * np.pi * fc2 * n_centered) - np.sin(2 * np.pi * fc1 * n_centered)) / (np.pi * n_centered)

    return h

def ideal_bandstop(N, fc1, fc2):
    """Calcula filtro rejeita-banda ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)

    for i in range(N):
        n_centered = n[i] - M / 2
        if abs(n_centered) < 1e-10:
            h[i] = 1 - 2 * (fc2 - fc1)
        else:
            pass_all = np.sin(np.pi * n_centered) / (np.pi * n_centered)
            stop_band = (np.sin(2 * np.pi * fc2 * n_centered) - np.sin(2 * np.pi * fc1 * n_centered)) / (np.pi * n_centered)
            h[i] = pass_all - stop_band

    return h
