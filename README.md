# Beamforming adaptativo em arranjo circular uniforme

Simulações MATLAB de beamforming adaptativo para sinais QAM em um arranjo circular uniforme (UCA) de dipolos com acoplamento mútuo. Os scripts comparam Wiener, LMS, RLS, CMA, CMA-DD concorrente e CMA-SDD.

## Cenários e modelo

- **Estático:** sinal desejado (DS) e interferentes (IS) com posições angulares fixas.
- **Dinâmico:** DS com movimento circular de velocidade angular constante em torno do arranjo; interferentes fixos.

Os dipolos são paralelos ao eixo z, com centros no plano xy. O modelo de recepção usa a aproximação de campo distante e banda estreita. As impedâncias próprias e mútuas são calculadas pelas expressões analíticas das notas técnicas de Fernando de Castro, e as cargas dos receptores entram na resposta do arranjo.

Os sinais QAM são gerados com potência média unitária. A entrada de cada antena resulta da soma das contribuições do DS e dos IS, seguida da adição de ruído gaussiano complexo independente por antena. As potências dos interferentes são ajustadas para obter a SINR de entrada especificada, dividindo igualmente a potência interferente recebida entre as fontes.

A resposta usada para gerar as amostras é normalizada e não inclui explicitamente a perda com a distância TX–RX. O parâmetro `r` dos scripts é usado no cálculo do campo para os diagramas. Os diagramas dos filtros mostram o campo normalizado em dB, no corte horizontal (θ = 90°); o plot de projeto do arranjo apresenta ganho em dBi.

## Estrutura

```text
Codes/
  UCA_Static_Beamforming.m    Cenário estático
  UCA_Dynamic_Beamforming.m   Cenário dinâmico
  Auxiliar/                  Modelagem, filtros, plots e animações
Figuras/                     Figuras de resultados
```

## Requisitos e execução

O ambiente utilizado é o **MATLAB R2025b**, com Communications Toolbox (`qammod`, `qamdemod`) e Symbolic Math Toolbox (`sinint`, `cosint`, `eulergamma`). A função `physconst` também deve estar disponível na instalação. O cálculo de impedâncias implementado não utiliza Antenna Toolbox.

Com a raiz do repositório como **Current Folder**, execute o cenário escolhido:

```matlab
% Cenário estático
run(fullfile('Codes', 'UCA_Static_Beamforming.m'))

% Ou cenário dinâmico
run(fullfile('Codes', 'UCA_Dynamic_Beamforming.m'))
```

Cada script limpa o workspace e fecha as figuras anteriores. Os resultados permanecem no workspace após a execução. Algumas chamadas de plots podem estar comentadas; descomente as que desejar visualizar.

Os principais parâmetros ficam no início e nas seções de configuração dos scripts:

| Parâmetros | Significado |
| --- | --- |
| `f`, `K`, `L`, `b`, `R`, `ZL` | Frequência, número de dipolos, comprimento, raio do fio, raio da UCA e cargas |
| `Dstream_size`, `Dmod_order` | Quantidade de símbolos e ordem QAM do DS |
| `phiDS` ou `phiDS_0` | Azimute fixo ou inicial do DS, em radianos |
| `IS_config` | Tabela de interferentes: azimute em graus e ordem QAM; use `[]` para nenhum IS |
| `SNR_RX`, `SINR_target_dB` | SNR e SINR nominais de entrada, em dB |
| `SR`, `ang_vel` | Taxa de símbolos [símbolos/s] e velocidade angular [rad/s], no cenário dinâmico |
| `stepLMS`, `lambdaRLS`, `deltaRLS` | Parâmetros dos filtros supervisionados |
| `stepCMA` e passos do Concurrent | Parâmetros dos algoritmos cegos |
| `stepCMASDDCMA`, `stepCMASDDSDD`, `rhoCMASDD` | Passos CMA/SDD e variância dos clusters por eixo real |
| `patternFloor_dB`, `mseWindow` | Limite inferior dos diagramas e janela de média do MSE |

Com interferentes presentes, escolha uma SINR alvo menor que a SNR para que a potência interferente calculada seja positiva. Sem interferentes, a SINR coincide com a SNR. A configuração atual não fixa a semente aleatória; execuções diferentes geram realizações diferentes dos símbolos e do ruído.

## Resultados e critérios de avaliação

As figuras permitem acompanhar MSE em dB, constelação, coeficientes e diagrama normalizado, além de comparar LMS × RLS e CMA × Concurrent. Nos dois cenários QAM, a SER de cada método conta todos os símbolos de dados do frame, identificados pela mesma `validationMask`. O MSE é exibido durante todo o registro. Os gráficos de constelação que usam `postStart` continuam mostrando o trecho final.

O cenário estático reserva os primeiros 20% do frame para treinamento. No dinâmico, a mesma fração é reservada no início de cada um dos dez subframes. Wiener estima os pesos em lote usando os símbolos de treinamento e aplica uma solução fixa ao frame inteiro. LMS e RLS atualizam os coeficientes durante o treinamento e os mantêm fixos durante os símbolos de dados. A matriz `P` do RLS também é mantida fixa nesses intervalos.

CMA, Concurrent e CMA-SDD são adaptados sem a sequência desejada. Para avaliar suas saídas, a correção de ambiguidade é feita **offline**, testando as rotações de 0°, 90°, 180° e 270°. A rotação é escolhida pela menor SER nos `ambiguitySymbols` símbolos iniciais (até 1000 símbolos do primeiro intervalo de treinamento) e aplicada ao registro inteiro.

O CMA-SDD usa quatro símbolos por região e atualiza os dois termos a partir da mesma saída. `rhoCMASDD` está na escala da constelação unitária: o valor inicial é `0.4/(2*(Dmod_order-1)/3)`, equivalente a `0.04` para 16-QAM e aproximadamente `0.009524` para 64-QAM. Cada cenário gera uma janela própria de resultados do método.

## Animações do cenário dinâmico

Após executar o script dinâmico, mantendo seus resultados no workspace, rode pela Command Window:

```matlab
Animate_Wiener_Results();
Animate_LMS_Results();
Animate_RLS_Results();
Animate_CMA_Results();
Animate_Concurrent_Results();
```

Cada função abre uma janela com MSE, SER por janela, constelação recente e diagrama com a posição atual do DS. A reprodução usa os resultados existentes. Os controles permitem pausar, reiniciar, avançar no tempo, alterar a velocidade e salvar um GIF.

```matlab
Animate_LMS_Results('Duration', 30, 'WindowSize', 5000);

Animate_Concurrent_Results( ...
    'GIFFile', fullfile(pwd, 'Figuras', 'Concurrent.gif'), ...
    'Duration', 10, 'FPS', 15, 'AutoPlay', false);
```

| Opção | Efeito | Padrão |
| --- | --- | --- |
| `Duration` | Duração da reprodução/exportação, em segundos | 20 |
| `FPS` | Quadros por segundo, entre 1 e 60 | 20 |
| `WindowSize` | Número de símbolos da janela da SER e da constelação | 2000 |
| `AutoPlay` | Iniciar automaticamente a reprodução | `true` |
| `GIFFile` | Caminho para exportar o GIF ao abrir a janela | Vazio |

A SER animada usa escala logarítmica; janelas sem erros aparecem em uma faixa identificada como `0 errors`. O GIF contém os quatro painéis, sem os controles, e também pode ser exportado pelo botão **Save GIF**.

## Referências

- Müller e de Castro, *Adaptive Blind Antenna Array Beam Shaping* (2025): [artigo no SBrT](https://biblioteca.sbrt.org.br/articlefile/4837.pdf).
- Fernando de Castro, [notas técnicas sobre arranjos de antenas e beamforming](https://www.fccdecastro.com.br/download.html).
