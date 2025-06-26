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
        fpassband1 = float(filter_parameters.fpassband1_var.get())
        fpassband2 = float(filter_parameters.fpassband2_var.get())
        transition_width = float(filter_parameters.transition_width_var.get())
        stopband_atten = float(filter_parameters.stopband_atten_var.get())
        filter_type = filter_parameters.filter_type_var.get()
        window_name = filter_parameters.selected_window_var.get()

        fc1, fc2 = 100, 1000  # Inicializar frequências de corte

        # Calcular frequência de corte baseada no tipo de filtro
        if filter_type == "Passa-Baixa":
            fstopband = fpassband1 + transition_width  # Frequência de stopband
            fc1 = (fpassband1 + fstopband) / 2  # Frequência de corte centrada na transição
        elif filter_type == "Passa-Alta":  # Passa-Alta
            fstopband = fpassband1 - transition_width  # Frequência de stopband
            fc1 = (fpassband1 + fstopband) / 2  # Frequência de corte
            if fstopband <= 0:
                raise ValueError("Para passa-alta: fp - largura_transição deve ser > 0")
        elif filter_type == "Passa-Banda":
            if fpassband1 >= fpassband2:
                raise ValueError("Frequência inferior deve ser menor que a superior")
            if fpassband2 >= fs/2:
                raise ValueError("Frequência superior deve ser menor que Fs/2")
            if fpassband1 <= 0:
                raise ValueError("Frequência inferior deve ser positiva")
            
            fstopband = fpassband1 - transition_width  # Frequência de stopband inferior
            fstopband2 = fpassband2 + transition_width # Frequência de stopband superior
            fc1 = (fpassband1 + fstopband) / 2  # Frequência de corte inferior
            fc2 = (fpassband2 + fstopband2) / 2  # Frequência de corte superior

            # Verificar se as frequências de corte estão dentro do intervalo válido
            if fpassband1 <= 0:
                raise ValueError("Para passa-banda: fp1 - largura_transição deve ser > 0")
            if fpassband2 >= fs/2:
                raise ValueError("Para passa-banda: fp2 + largura_transição deve ser < Fs/2")

        elif filter_type == "Rejeita-Banda":
            if fpassband1 >= fpassband2:
                raise ValueError("Frequência inferior deve ser menor que a superior")
            if fpassband2 >= fs/2:
                raise ValueError("Frequência superior deve ser menor que Fs/2")
            if fpassband1 <= 0:
                raise ValueError("Frequência inferior deve ser positiva")
                
            fstopband = fpassband1 + transition_width
            fstopband2 = fpassband2 - transition_width
            if fstopband <= 0:
                raise ValueError("Para rejeita-banda: fp1 - largura_transição deve ser > 0")
            if fpassband2 >= fs/2:
                raise ValueError("Para rejeita-banda: fp2 + largura_transição deve ser < Fs/2")
                
            fc1 = (fstopband + fpassband1) / 2
            fc2 = (fstopband2 + fpassband2) / 2
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
        
        # Projetar filtro ideal normalizado em termos de Nyquist para garantir que ele foi aplicado 
        # corretamente e que não haverá aliasing ou distorção de fase
        nyquist = fs / 2
        fc1 = fc1 / nyquist
        fc2 = fc2 / nyquist if filter_type in ["Passa-Banda", "Rejeita-Banda"] else None
        
        # Calcular resposta ideal do filtro
        
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
            window = signal.get_window(('kaiser', beta), order, fftbins=False)
        elif window_name == 'Bartlett':
            window = signal.get_window('bartlett', order, fftbins=False)
        elif window_name == 'Hanning':
            window = signal.get_window('hann', order, fftbins=False)
        elif window_name == 'Hamming':
            window = signal.get_window('hamming', order, fftbins=False)
        elif window_name == 'Blackman':
            window = signal.get_window('blackman', order, fftbins=False)
        
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

#O filtro é centrado em M/2 para simetria tornando o filtro causal
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

#O filtro é centrado em M/2 para simetria tornando o filtro causal
def ideal_highpass(N, fc_norm):
    """Calcula filtro passa-alta ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)
    
    for i in range(N):
        n_centered = n[i] - M/2
        
        if abs(n_centered) < 1e-10:
            h[i] = 1.0 - fc_norm
        else:
            # Impulso delta menos passa-baixa
            delta_impulse = np.sin(np.pi * n_centered) / (np.pi * n_centered)
            sinc_arg = fc_norm * n_centered
            lowpass = fc_norm * np.sin(np.pi * sinc_arg) / (np.pi * sinc_arg)
            h[i] = delta_impulse - lowpass
    
    return h

#O filtro é centrado em M/2 para simetria tornando o filtro causal
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
            omega_c1 = fc1 * np.pi
            omega_c2 = fc2 * np.pi
            h[i] = (2 * fc2 * np.sin(omega_c2 * n_centered) / (omega_c2 * n_centered) - 
                    2 * fc1 * np.sin(omega_c1 * n_centered) / (omega_c1 * n_centered))

    return h

#O filtro é centrado em M/2 para simetria tornando o filtro causal
def ideal_bandstop(N, fc1, fc2):
    """Calcula filtro rejeita-banda ideal"""
    M = N - 1
    n = np.arange(N)
    h = np.zeros(N)

    for i in range(N):
            n_centered = n[i] - M / 2
            
            if abs(n_centered) < 1e-10:
                h[i] = 1 - (fc2 - fc1)
            else:
                
                pass_all = np.sin(np.pi * n_centered) / (np.pi * n_centered)
                stop_band = (np.sin(np.pi * fc2 * n_centered) - 
                           np.sin(np.pi * fc1 * n_centered)) / (np.pi * n_centered)
                h[i] = pass_all - stop_band

    return h
