# Projeto 2: Análise de Evaporadores de Filme Descendente - Transferência de Calor (UFRJ)

**Autor:** Gabriel José Costa Teles  
**Data:** Novembro/2025

## Descrição
Este projeto realiza a modelagem computacional dos perfis térmicos e hidrodinâmicos em um evaporador de filme descendente. O código foi desenvolvido em Python e baseia-se nas correlações teóricas e experimentais apresentadas no artigo de *Fang et al. (2019)*.

## Estrutura do Código
O script `simulacao_evaporador.py` está dividido nas seguintes etapas:

1.  **Definição de Parâmetros:** Estabelecimento das constantes geométricas ($L=2.0m$, $D=0.081m$) e propriedades térmicas baseadas no artigo.
2.  **Modelagem Física:** Implementação das funções para cálculo do Número de Nusselt (Eq. 3.10 de Fang et al.) e estimativas de espessura de filme laminar e turbulento.
3.  **Geração Gráfica (Item 4):**
    * **Item 4A:** Visualização dos perfis de temperatura usando mapas de cor (`plasma`) para diferenciar a parede quente (evaporação) da fria (condensação).
    * **Item 4B:** Análise de sensibilidade do Número de Reynolds, demonstrando a transição de regimes e afinamento do filme.
    * **Item 4C:** Mapeamento do Número de Prandtl através de curvas de contorno.
    * **Item 4D:** Validação da relação inversa entre espessura do filme e eficiência térmica ($Nu$).

## Como Executar
O código requer as bibliotecas `numpy`, `matplotlib` e `scipy`.
