# SPH

## Abstract

Este repositório apresenta a implementação do método Smoothed Particle Hydrodynamics (SPH) no contexto de simulação numérica de fluidos. O projeto tem como objetivo explorar a modelagem computacional de meios contínuos por meio de uma abordagem Lagrangiana baseada em partículas.

A implementação busca reproduzir o comportamento físico de fluidos através da discretização do domínio em partículas, utilizando funções de suavização (kernels) para aproximar grandezas físicas e resolver equações fundamentais da dinâmica dos fluidos.

---

## 1. Introdução

A simulação de fluidos é um problema central na física computacional e na engenharia. Métodos tradicionais baseados em malhas apresentam limitações em cenários com grandes deformações.

O método Smoothed Particle Hydrodynamics (SPH) surge como uma alternativa mesh-free, na qual o fluido é representado por partículas que carregam propriedades físicas e evoluem ao longo do tempo.

Este projeto tem como objetivo implementar e analisar esse método, explorando seus aspectos computacionais e numéricos.

---

## 2. Fundamentação Teórica

O método SPH baseia-se na aproximação de funções contínuas por meio de interpolação com kernels:

- Cada partícula contribui para o valor de uma grandeza em sua vizinhança
- Interações locais determinam propriedades globais
- O comportamento do fluido emerge das interações entre partículas

As equações fundamentais envolvidas incluem:

- Conservação de massa
- Conservação de momento (equação de Navier-Stokes)
- Equações de estado

A discretização é realizada por meio de somatórios sobre partículas vizinhas.

---

## 3. Estrutura do Repositório
'''sph/
│── src/ # Código-fonte principal
│── simulations/ # Configurações ou cenários
│── README.md


A organização segue uma separação entre lógica do simulador e definição de cenários experimentais.

---

## 4. Arquitetura do Sistema

### 4.1 Sistema de Partículas

Responsável por representar o fluido como um conjunto discreto de partículas.

Funcionalidades:
- Armazenamento de posição, velocidade e massa
- Atualização temporal das partículas
- Estruturas de dados para acesso eficiente

---

### 4.2 Cálculo de Vizinhança

Módulo responsável por identificar partículas próximas.

Importância:
- Redução da complexidade computacional
- Base para cálculo de densidade e forças

---

### 4.3 Kernels de Suavização

Implementa funções de interpolação utilizadas no SPH.

Função:
- Aproximar funções contínuas
- Determinar influência entre partículas

---

### 4.4 Cálculo de Densidade e Pressão

Responsável por calcular propriedades físicas locais.

Inclui:
- Soma ponderada das contribuições vizinhas
- Equação de estado para pressão

---

### 4.5 Integração Temporal

Atualiza o sistema ao longo do tempo.

Métodos possíveis:
- Euler explícito
- Runge Kutta
- Leapfrog

---

### 4.6 Cálculo de Forças

Responsável pela dinâmica do sistema.

Inclui:
- Forças de pressão
- Viscosidade
- Forças externas (ex: gravidade)

---

## 5. Metodologia

O desenvolvimento seguiu as etapas:

1. Modelagem matemática do problema
2. Discretização do fluido em partículas
3. Implementação incremental dos componentes do SPH
4. Validação qualitativa dos resultados
5. Ajuste de parâmetros (densidade, viscosidade, smoothing length)

---

## 6. Resultados e Discussão

A implementação permite observar:

- Comportamento emergente do fluido
- Sensibilidade a parâmetros físicos
- Impacto do número de partículas na estabilidade e precisão

Observa-se que:

- A escolha do kernel influencia diretamente a estabilidade
- O custo computacional cresce com o número de partículas
- Estruturas eficientes de vizinhança são essenciais

---

## 7. Limitações

- Alto custo computacional para grandes sistemas
- Sensibilidade numérica a parâmetros
- Necessidade de otimizações (ex: paralelismo)

---

## 8. Conclusão

O projeto demonstra a viabilidade do método SPH como ferramenta para simulação de fluidos, destacando sua flexibilidade e capacidade de lidar com domínios complexos.

Além disso, evidencia os desafios computacionais associados, especialmente em termos de desempenho e estabilidade numérica.

---


