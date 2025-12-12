"""
PROJETO 2 - TRANSCAL 2: ANÁLISE DE EVAPORADORES DE FILME DESCENDENTE
Versão final consolidada - UFRJ
Autor: Gabriel José Costa Teles
Data: Novembro/2025

DESCRIÇÃO GERAL:
Este script simula o comportamento térmico e hidrodinâmico de filmes descendentes.
Baseia-se nas equações e resultados numéricos apresentados por:
    Fang, J., Li, K., & Diao, M. (2019). "Establishment of the falling film evaporation 
    model and correlation of the overall heat transfer coefficient."

ITENS GERADOS:
 - Item 4A: Esquema visual das camadas limite (Térmica vs Hidrodinâmica) com perfis de T e u.
 - Item 4B: Gráficos de espessura do filme variando Re (Rastejante, Laminar, Turbulento).
 - Item 4C: Mapa de contorno mostrando Pr em função de Re e Nu.
 - Item 4D: Relação entre resistência térmica (espessura) e Nusselt.

ARQUIVOS DE SAÍDA:
    Item4A_with_Tc.png, Item4B_3subplots.png, Item4C_Pr_Re_vs_Nu.png, Item4D_Nu_vs_delta.png
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker

# SciPy griddata (usado em Item 4C para interpolação da superfície de contorno)
from scipy.interpolate import griddata

# ------------------------
# CONSTANTES GLOBAIS
# ------------------------
# Parâmetros geométricos e propriedades baseados na Tabela 2 e texto de Fang et al. (2019)
L_PLATE = 2.0           # Comprimento do tubo/placa [m]
D_0 = 0.081             # Diâmetro externo [m]
K_WALL = 16.0           # Condutividade térmica da parede [W/m·K]

# Constante da correlação de Nusselt (Coeficiente 'a' na Eq. 3.8/3.9 do artigo)
CONST_NU = 5.32e4       
Z_AXIS = np.linspace(0, L_PLATE, 200)  # Eixo vertical discretizado para os esquemas visuais

# ------------------------
# TEMPERATURAS (°C) conforme Fang et al. (usadas em 4A)
# ------------------------
# Lado da Evaporação (Interno) - Temperaturas mais altas
Tw_evap = 150.0   # Parede quente [°C]
Ts_evap = 100.0   # Interface líquido-vapor (Saturação) [°C]

# Lado da Condensação (Externo) - Temperaturas mais baixas
Tw_cond = 40.0    # Parede fria [°C]
Ts_cond = 100.0   # Interface líquido-vapor (Saturação) [°C]

# Limites globais para a colorbar (garante que a escala de cor cubra todo o range térmico)
T_vmin = min(Tw_cond, Ts_evap, Ts_cond)
T_vmax = max(Tw_evap, Ts_evap, Ts_cond)

# ------------------------
# FUNÇÕES FÍSICAS
# ------------------------

def calculate_nusselt_real(Re, Pr):
    """
    Calcula o Número de Nusselt para condensação externa.
    Baseado na Eq. 3.10 de Fang et al. (2019):
    Nu = 5.32e4 * Re^(-0.1418) * Pr^(-3.1975) * (D_0 / K_WALL)
    
    Args:
        Re: Número de Reynolds do filme
        Pr: Número de Prandtl do condensado
    Returns:
        Nu: Número de Nusselt
    """
    geo_factor = D_0 / K_WALL
    return CONST_NU * (Re**-0.1418) * (Pr**-3.1975) * geo_factor


def calculate_delta_ratio_laminar(Re):
    """
    Relação adimensional da espessura da camada limite para regime Laminar.
    Teoria clássica de Nusselt: delta/x ~ Re^-0.5
    """
    return 5.0 / np.sqrt(Re)


def calculate_delta_ratio_turbulent(Re):
    """
    Relação adimensional da espessura da camada limite para regime Turbulento.
    Correlação empírica típica: delta/x ~ Re^-0.2
    """
    return 0.37 / (Re**0.2)


def calculate_film_thickness_estimate(Re):
    """
    Estimativa da espessura física do filme em metros.
    Usa a relação de Nusselt (delta ~ Re^(1/3)) escalada para visualização.
    """
    return (Re)**(1/3) * 0.001

# ------------------------
# PLOT ITEM 4A (ESQUEMÁTICO)
# ------------------------

def plot_item_4a_physics(ax, Pr, title):
    """
    Gera o esquema lateral comparando Evaporação (Esq) e Condensação (Dir).
    
    FÍSICA IMPLEMENTADA:
    1. Evaporação: Perfil de temperatura Cúbico (simula efeitos convectivos/ebulição).
    2. Condensação: Perfil de temperatura Linear (hipótese de condução pura, seção 2.1 do artigo).
    3. Velocidade: Perfil parabólico (laminar) com vetores (quiver).
    4. Mapa de Cores: 'plasma' (alto contraste).
    """
    # -- Geometria visual qualitativa --
    # Evaporação: Filme afina (perda de massa)
    delta_in  = 1.0 * ((1.05 - Z_AXIS / L_PLATE)**0.5)      
    # [cite_start]Condensação: Filme engrossa (ganho de massa por condensação) [cite: 131]
    delta_out = 1.0 * ((Z_AXIS / L_PLATE)**0.25 + 0.02)     
    
    # Relação entre espessura térmica (delta_t) e hidrodinâmica (delta)
    # Escala com Pr^(-1/3)
    ratio_Pr = Pr**(-1/3)
    delta_t_in  = delta_in * ratio_Pr
    delta_t_out = delta_out * ratio_Pr

    # Configuração dos eixos (remove bordas para focar no esquema)
    max_w = max(delta_out.max(), delta_t_out.max()) * 1.6
    ax.set_ylim(L_PLATE, 0)
    ax.set_xlim(-max_w, max_w)
    ax.axis('off')
    ax.set_title(f"{title}\n(Pr = {Pr})", fontsize=11)

    # Desenha a parede central do tubo
    ax.axvline(0, color='black', lw=3, zorder=10)

    # ==========================================
    # LADO ESQUERDO: Evaporação (Parede Quente)
    # ==========================================
    x_vals_in = np.linspace(-max_w, 0, 250)
    XX_in, ZZ_in = np.meshgrid(x_vals_in, Z_AXIS)
    dist_wall_in = np.abs(XX_in)

    # Interpolação das espessuras na malha
    d_h_local_in = np.interp(ZZ_in, Z_AXIS, delta_in)
    d_t_local_in = np.interp(ZZ_in, Z_AXIS, delta_t_in)
    y_norm_in = np.clip(dist_wall_in / (d_t_local_in + 1e-12), 0, 1)

    # Perfil térmico CÚBICO (Evaporação)
    # Representa o gradiente mais acentuado na parede devido à convecção
    T_field_in = Ts_evap + (Tw_evap - Ts_evap) * (1 - y_norm_in)**3
    
    # Mascara para plotar apenas dentro do filme
    T_plot_in = np.where(dist_wall_in <= d_h_local_in, T_field_in, np.nan)

    # Plot Temperatura (PLASMA)
    ax.pcolormesh(XX_in, ZZ_in, T_plot_in, cmap='plasma', shading='auto',
                  vmin=T_vmin, vmax=T_vmax, alpha=0.9, zorder=1)

    # Linhas de contorno (Preto=Hidro, Azul=Térmica)
    ax.plot(-delta_in, Z_AXIS, 'k-', lw=2, zorder=8)       
    ax.plot(-delta_t_in, Z_AXIS, 'b--', lw=1.5, zorder=8)  

    # --- PERFIL DE VELOCIDADE (ESQUERDA) ---
    z_prof = L_PLATE * 0.3
    w_loc = np.interp(z_prof, Z_AXIS, delta_in)
    y_vect = np.linspace(0, w_loc, 24)
    y_n = y_vect / w_loc
    
    # [cite_start]Perfil parabólico laminar: u ~ 2y - y^2 [cite: 120]
    u_mag = (2*y_n - y_n**2) * 0.25

    # Desenha a linha do perfil e os vetores (setas)
    ax.plot(-y_vect, z_prof + u_mag, 'k-', lw=1.6, zorder=12)
    ax.quiver(-y_vect[::3], np.full_like(y_vect[::3], z_prof),
              np.zeros_like(y_vect[::3]), u_mag[::3],
              angles='xy', scale_units='xy', scale=1,
              width=0.003, headwidth=4, headlength=5,
              color='k', zorder=12)

    # ==========================================
    # LADO DIREITO: Condensação (Parede Fria)
    # ==========================================
    x_vals_out = np.linspace(0, max_w, 250)
    XX_out, ZZ_out = np.meshgrid(x_vals_out, Z_AXIS)
    dist_wall_out = XX_out

    d_h_local_out = np.interp(ZZ_out, Z_AXIS, delta_out)
    d_t_local_out = np.interp(ZZ_out, Z_AXIS, delta_t_out)
    y_norm_out = np.clip(dist_wall_out / (d_t_local_out + 1e-12), 0, 1)

    # Perfil térmico LINEAR (Condensação)
    # [cite_start]Hipótese 6 do artigo: condução pura através do filme [cite: 81, 143]
    T_field_out = Tw_cond + (Ts_cond - Tw_cond) * y_norm_out
    T_plot_out = np.where(dist_wall_out <= d_h_local_out, T_field_out, np.nan)

    # Plot Temperatura (PLASMA)
    ax.pcolormesh(XX_out, ZZ_out, T_plot_out, cmap='plasma', shading='auto',
                  vmin=T_vmin, vmax=T_vmax, alpha=0.9, zorder=1)

    # Linhas de contorno
    ax.plot(delta_out, Z_AXIS, 'k-', lw=2, zorder=8)
    ax.plot(delta_t_out, Z_AXIS, 'b--', lw=1.5, zorder=8)

    # --- PERFIL DE VELOCIDADE (DIREITA) ---
    z_prof_r = L_PLATE * 0.7
    w_loc_r = np.interp(z_prof_r, Z_AXIS, delta_out)
    y_vect_r = np.linspace(0, w_loc_r, 24)
    y_nr = y_vect_r / w_loc_r
    
    # Perfil parabólico laminar
    u_mag_r = (2*y_nr - y_nr**2) * 0.25

    ax.plot(y_vect_r, z_prof_r + u_mag_r, 'k-', lw=1.6, zorder=12)
    ax.quiver(y_vect_r[::3], np.full_like(y_vect_r[::3], z_prof_r),
              np.zeros_like(y_vect_r[::3]), u_mag_r[::3],
              angles='xy', scale_units='xy', scale=1,
              width=0.003, headwidth=4, headlength=5,
              color='k', zorder=12)

    # ==========================================
    # ROTULAGEM
    # ==========================================
    # Adiciona legenda indicando a relação entre delta_t e delta
    if Pr == 1.0:
        rel, col = r"$\delta_t = \delta$", 'black'
    elif Pr < 1.0:
        rel, col = r"$\delta_t > \delta$", 'blue'
    else:
        rel, col = r"$\delta_t < \delta$", 'red'

    ax.text(0, L_PLATE * 0.92, rel, color=col, fontsize=11, ha='center', fontweight='bold', 
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

    ax.text(-max_w*0.5, L_PLATE*1.02, "Evaporação", ha='center', fontsize=9, fontweight='bold')
    ax.text(max_w*0.5, L_PLATE*1.02, "Condensação", ha='center', fontsize=9, fontweight='bold')


# ------------------------
# ITEM 4B – 3 subplots (Influência de Re)
# ------------------------

def plot_item_4b_three_subplots():
    """
    Gera 3 gráficos comparando a evolução da camada limite para diferentes Reynolds.
    - Re << 1 (Rastejante)
    - Re = 1 (Laminar Padrão)
    - Re >> 1 (Transição/Turbulento)
    
    NOTA: Trava o eixo Y dos dois primeiros casos na mesma escala para evidenciar 
    que o aumento de Re diminui a espessura adimensional do filme.
    """
    Re_cases = [0.1, 1.0, 3e5] 
    titles = [r"$Re \ll 1$ (Rastejante)", r"$Re = 1$ (Laminar)", r"$Re \gg 1$ (Transição/Turb)"]
    Pr_ref = 0.7 
    
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"Item 4B: Influência de Re na Camada Limite (Pr={Pr_ref})", fontsize=14)
    
    x_dim = np.linspace(0.001, 1.0, 500)
    Re_crit_local = 2e5 # Critério visual para transição

    summary = []
    
    # Pré-cálculo para pegar o limite máximo do primeiro gráfico (Re=0.1)
    # Isso serve de referência para travar a escala dos plots laminares
    Re_lowest = Re_cases[0]
    delta_max_ref = 5.0 * 1.0 / np.sqrt(Re_lowest * 1.0) 
    y_lim_laminar = delta_max_ref * 1.1 # +10% de margem

    for ax, ReL, title in zip(axs, Re_cases, titles):
        Re_x = ReL * x_dim
        
        # Correlações teóricas para espessura
        delta_lam_dim = 5.0 * x_dim / np.sqrt(Re_x)
        delta_tur_dim = 0.37 * x_dim / (Re_x**0.2)
        
        # Lógica de Transição (Visual)
        if ReL < 100:
            delta_final = delta_lam_dim
        else:
            # Cria um "degrau" visual onde o escoamento transiciona
            mask_lam = Re_x < Re_crit_local
            delta_final = np.where(mask_lam, delta_lam_dim, delta_tur_dim)

        # Espessura térmica aproximada
        delta_t_final = delta_final * (Pr_ref**(-1/3))
        
        # Plotagem (com preenchimento 'fill_between' para realçar o filme)
        ax.fill_between(x_dim, 0, delta_final, color='skyblue', alpha=0.4, label=r'Filme $\delta$')
        ax.plot(x_dim, delta_final, 'b-', lw=2, label=r'$\delta$ (Hidro)')
        ax.plot(x_dim, delta_t_final, 'r--', lw=1.5, label=r'$\delta_t$ (Térmica)')
        
        # --- Lógica de Escala Fixa ---
        if ReL <= 1.0:
            # Trava a escala no máximo do caso Re=0.1
            ax.set_ylim(0, y_lim_laminar)
            ax.text(0.05, y_lim_laminar*0.9, "Escala Travada", fontsize=8, color='gray')
        else:
            # Para Re turbulento (filme muito fino), usa zoom automático fixo
            ax.set_ylim(0, 0.05)
            # Linha vertical indicando onde ocorre a transição
            if np.any(Re_x > Re_crit_local):
                x_trans = x_dim[np.searchsorted(Re_x, Re_crit_local)]
                ax.axvline(x_trans, color='k', linestyle=':', alpha=0.5)

        ax.set_xlabel('x/L')
        ax.set_ylabel(r'$\delta/L$')
        ax.set_title(f"{title}\n$Re_L={ReL:g}$")
        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)

        ratio_avg = np.mean(delta_final / delta_t_final)
        summary.append((ReL, ratio_avg))

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("Item4B_3subplots.png", dpi=300)


# ------------------------
# ITEM 4C – Contorno de Pr em Re × Nu
# ------------------------

def plot_item_4c_contour_pr():
    """
    Gera um gráfico de contorno (mapa de cores) para o Número de Prandtl.
    
    Procedimento:
      1) Cria uma malha de cálculo (Re, Pr) e calcula Nu usando Eq. 3.10.
      2) Inverte o mapeamento interpolando Pr em uma nova malha (Re, Nu).
      3) Plota o contorno de Pr sobre os eixos Re e Nu.
    """

    # 1. Malha paramétrica original (RE, PR) -> NU
    re_vals = np.linspace(100, 2000, 300)
    pr_vals = np.linspace(0.6, 12.0, 200)
    REm, PRm = np.meshgrid(re_vals, pr_vals)
    NU_m = calculate_nusselt_real(REm, PRm)

    # 2. Malha de destino para o gráfico (Re no eixo X, Nu no eixo Y Logarítmico)
    nu_vals = np.logspace(np.log10(max(NU_m.min(), 1e-6)), np.log10(NU_m.max()), 300)
    RE2, NU2 = np.meshgrid(re_vals, nu_vals)

    # Prepara pontos para interpolação: (Re, Nu) -> Pr
    pts = np.column_stack((REm.flatten(), NU_m.flatten()))
    vals = PRm.flatten()

    # Interpolação Linear (Griddata)
    PR_interp = griddata(pts, vals, (RE2, NU2), method='linear')
    PR_interp = np.ma.masked_invalid(PR_interp)

    # 3. Plotagem
    fig, ax = plt.subplots(figsize=(9, 6))
    levels = np.linspace(pr_vals.min(), pr_vals.max(), 20)
    
    # Contorno Preenchido (contourf)
    cf = ax.contourf(RE2, NU2, PR_interp, levels=levels, cmap='viridis', extend='both')
    
    # Linhas de Contorno (contour) - Pedido na especificação
    cl = ax.contour(RE2, NU2, PR_interp, levels=levels, colors='k', linewidths=0.6)
    ax.clabel(cl, fmt='%.2f', fontsize=8)

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label('Prandtl (Pr)', rotation=270, labelpad=15)

    ax.set_xlabel('Reynolds (Re)')
    ax.set_ylabel('Nusselt (Nu) [Escala Log]')
    ax.set_yscale('log')
    ax.set_title('Item 4C – Contorno de Pr em Re × Nu')
    ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.savefig("Item4C_Pr_Re_vs_Nu.png", dpi=300)


# ------------------------
# ITEM 4D – Nu vs Espessura
# ------------------------

def plot_item_4d():
    """
    Gráfico relacionando o Número de Nusselt com a espessura estimada do filme.
    Mostra que filmes mais espessos (maior resistência térmica) levam a menores coeficientes de troca (Nu).
    """
    re_vals = np.linspace(100, 1500, 200)
    # Calcula espessura estimada para eixo X (delta ~ Re^(1/3))
    delta_sim = calculate_film_thickness_estimate(re_vals)

    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plota curvas para diferentes Pr (Ar, Água, Óleo)
    for pr_val, col in zip([1.0, 5.0, 10.0], ['purple', 'green', 'blue']):
        nu_val = calculate_nusselt_real(re_vals, pr_val)
        ax.plot(delta_sim, nu_val, color=col, lw=2.2, label=f'Pr={pr_val}')

    ax.set_yscale('log')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Espessura estimada do filme [m]')
    ax.set_ylabel('Nusselt (Nu)')
    ax.legend(title='Prandtl')
    ax.set_title('Item 4D – Nu vs Espessura do Filme')
    plt.tight_layout()
    plt.savefig("Item4D_Nu_vs_delta.png", dpi=300)


# ------------------------
# MAIN – Execução e Geração das Figuras
# ------------------------

if __name__ == "__main__":
    print(">>> Iniciando geração de figuras...")

    # ---- ITEM 4A: Esquema das Camadas Limite ----
    # Gera 3 subplots para diferentes Pr (0.7, 1.0, 100)
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 6))
    fig1.suptitle('Item 4A: Camadas Limite (Temperatura em °C - Mapa Plasma)', fontsize=16)

    PR_list = [0.7, 1.0, 100.0]   
    titles = ["Gases / Ar (Pr=0.7)", "Fluido Ideal / Água (Pr=1.0)", "Óleos Pesados (Pr=100)"]

    for i, (ax, Pr, title) in enumerate(zip(axes1, PR_list, titles)):
        plot_item_4a_physics(ax, Pr, title)

        # Colorbar individual para cada subplot
        norm = plt.Normalize(vmin=T_vmin, vmax=T_vmax)
        sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
        sm.set_array([])  
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Temperatura (°C)", fontsize=8)
        cbar.ax.tick_params(labelsize=8)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("Item4A_with_Tc.png", dpi=300)
    print(" [OK] Item 4A salvo (arquivo: Item4A_with_Tc.png).")

    # ---- ITEM 4B: Comparação de Regimes (Re) ----
    plot_item_4b_three_subplots()
    print(" [OK] Item 4B salvo (arquivo: Item4B_3subplots.png).")

    # ---- ITEM 4C: Mapa de Contorno ----
    plot_item_4c_contour_pr()
    print(" [OK] Item 4C salvo (arquivo: Item4C_Pr_Re_vs_Nu.png).")

    # ---- ITEM 4D: Nu vs Delta ----
    plot_item_4d()
    print(" [OK] Item 4D salvo (arquivo: Item4D_Nu_vs_delta.png).")

    print("\n>>> Todas as figuras foram geradas com sucesso.")
