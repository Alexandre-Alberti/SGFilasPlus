# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 12:01:26 2024

@author: alexa
"""

import numpy as np
from numpy import random as rd
import streamlit as stlt

stlt.title('SGFilas')
stlt.text('')
stlt.subheader('Análise de sistemas de fila com esquema de priorização')
stlt.text('')
stlt.text('Insira as informações solicitadas nos campos abaixo ')

#n_at = stlt.number_input("Número de servidores:")
#n_at = int(n_at)

n_grupos = stlt.number_input("Número de grupos para segmentação:")
n_grupos = int(n_grupos)

taxas = np.zeros(n_grupos) #taxa de chegadas de cada tipo de demanda

NS = np.zeros(n_grupos) #nível de serviço desejado para cada grupo, em % 

tol = np.zeros(n_grupos) #limite de tempo tolerado

T_zero = np.zeros(n_grupos) #tempo de espera para deslocamento da demanda para o grupo Zero (frente da fila)
#T_zero é a solução que eu busco

tipos = []

parametros = np.zeros((n_grupos,6))

opcoes_dist = ["Exponencial","Normal","Triangular"]

for i in range(0,n_grupos):
    
    stlt.write(f'<font color="blue" size="5"><b>Grupo {i+1}</b></font>', unsafe_allow_html=True)
    
    taxas[i] = stlt.number_input(f'Taxa de chegada de demandas (grupo {i+1}):')
    
    #NS[i] = stlt.number_input(f'Nível de serviço desejado (grupo {i+1}):')
    
    tol[i] = stlt.number_input(f'Tempo máximo de espera tolerado (grupo {i+1}):')
    
    NS = stlt.number_input(f'Nível de serviço desejado (grupo {i+1}):')
    NS[i] = NS/100    
    
    tipo_dist = stlt.selectbox(f'Distribuição de probabilidade para o tempo de atendimento (grupo {i+1}):', opcoes_dist)
    
    if tipo_dist == 'Exponencial':
        tipos.append("Exp")
        parametros[i,0] = stlt.number_input(f'Tempo médio de atendimento (grupo {i+1}):')
    if tipo_dist == 'Normal':
        tipos.append("Norm")
        parametros[i,1] = stlt.number_input(f'Média do tempo de atendimento (grupo {i+1}):')
        parametros[i,2] = stlt.number_input(f'Desvio-padrão (grupo {i+1}):')
    if tipo_dist == 'Triangular':
        tipos.append("Tri")
        parametros[i,3] = stlt.number_input(f'Tempo mínimo (grupo {i+1}):')
        parametros[i,4] = stlt.number_input(f'Moda (grupo {i+1}):')
        parametros[i,5] = stlt.number_input(f'Tempo máximo (grupo {i+1}):')


#função de apoio
def soma_parcial(a,n):
    soma = 0
    for i in range(0,n):
        soma = soma + a[i]        
    return soma

def sorteio_classe(a):
    p = rd.rand()
    n = len(a)
    if p < a[0]:
        classe = 1
    else:
        for i in range(1,n+1):
            if (p < soma_parcial(a,i)) and (p > soma_parcial(a,i-1)):
                classe = i
    return classe

def transpor(B):
    linhas = len(B)
    colunas = len(B[0])
    B_transposto = np.zeros((colunas,linhas))
    for i in range(0,linhas):
        for j in range(0,colunas):
            B_transposto[j][i] = B[i][j]
    return B_transposto

def ordenar_linha(B,n):
    B = transpor(B)
    B = sorted(B, key=lambda x: x[n])
    B = transpor(B)
    return B



#SIMULAÇÃO
def fila(n_atendentes, T_zero, horizonte):
    #transformando as taxas em percentual
    taxa_global = sum(taxas)
    taxa_global_inverso = 1/taxa_global
    taxas_percentual = np.zeros(n_grupos)    
    
    for i in range(0,n_grupos):
        taxas_percentual[i] = taxas[i]/taxa_global
    
    #definindo o periodo de tempo para simulação
    menor_taxa = min(taxas)
    
    #
    #
    #
    #janela_obs = 10000*(1/menor_taxa)
    janela_obs = horizonte*(1/menor_taxa) 
    
    #simular as chegadas de demandas ao longo do tempo
    #vou utilizar um laço tipo while, até alcançar a extensão da janela de observação
    parada = 0 #critério de parada
    tempo = 0 #inicio da contagem
    vetor_chegadas = ()
    
    while parada == 0:
        X = rd.exponential(taxa_global_inverso,1)
        tempo = tempo + X
        vetor_chegadas = np.append(vetor_chegadas,tempo)
        if tempo >= janela_obs:
            parada = 1
     
    num_chegadas = len(vetor_chegadas)

    #vou criar uma matriz, que concentre as informações de que preciso para uma primeira análise
    matriz_dados = np.zeros((8,num_chegadas))

    # primeiro linha, ordem de chegada
    # segunda linha, momento da chegada
    # terceira linha, grupo de origem da demanda
    # quarta linha por ora fica vazia, onde serão registrados o tempo de saída
    # quinta linha, 1 para demanda não atendida, 0 para demanda atendida
    # sexta linha, tempo de chegada + t_x
    # sétima linha, 1 para demanda atendida dentro do prazo tolerado, 0 caso contrario
    # oitava linha, tempo de atendimento
    
    for i in range(0,num_chegadas):
        
        matriz_dados[0][i] = i+1
        
        tempo_chegada = vetor_chegadas[i]
        matriz_dados[1][i] = tempo_chegada
        
        grupo_origem = sorteio_classe(taxas_percentual)
        matriz_dados[2][i] = grupo_origem
        
        #matriz_dados[3][i] - campos vazios por ora
        matriz_dados[4][i] = 1 # 1 denota demanda não atendida, para demandas atendidas se utilizará zero
        
        tempo_muda = T_zero[grupo_origem-1]
        matriz_dados[5][i] = tempo_chegada + tempo_muda
        
        #matriz_dados[6][i] - campos vazios por ora
        
        #gravando os tempos de atendimento
        if tipos[grupo_origem-1] == 'Exp':
            tempo_serv = rd.exponential(parametros[grupo_origem-1][0])
        if tipos[grupo_origem-1] == 'Norm':
            tempo_serv = rd.normal(parametros[grupo_origem-1][1],parametros[grupo_origem-1][2])
        if tipos[grupo_origem-1] == 'Tri':
            tempo_serv = rd.triagular(parametros[grupo_origem-1][3],parametros[grupo_origem-1][4],parametros[grupo_origem-1][5])
        matriz_dados[7][i] = tempo_serv
            
    #simulando atendimentos
    
    #para fins de controle
    
    grupos_ordem_chegada = [[] for _ in range(n_grupos)]
    grupos_tempo_chegada = [[] for _ in range(n_grupos)]
    
    for i in range(0, num_chegadas):
        nn = matriz_dados[2][i] #identificação do grupo
        nn = int(nn)
        grupos_ordem_chegada[nn-1] = np.append(grupos_ordem_chegada[nn-1],matriz_dados[0][i])
        grupos_tempo_chegada[nn-1] = np.append(grupos_tempo_chegada[nn-1],matriz_dados[1][i])
    #print('grupos_ordem',grupos_ordem_chegada)
    #print('grupos_tempo',grupos_tempo_chegada)
    
    for i in range(0, n_grupos):
        if len(grupos_ordem_chegada[i]) == 0:
            grupos_ordem_chegada[i] = [100*janela_obs]
            grupos_tempo_chegada[i] = [100*janela_obs]
    
    #Preciso ordenar os grupos de acordo com a prioridade, isso vai ser definido pelos valores de T (do menor para o maior)
    ordem_relevancia = np.zeros((3,n_grupos))
    for i in range(0,n_grupos):
        ordem_relevancia[1][i] = i+1  #numero de identificação do grupo
        ordem_relevancia[2][i] = T_zero[i]
    
    ordem_relevancia = ordenar_linha(ordem_relevancia, 2)
    for i in range(0,n_grupos):
        ordem_relevancia[0][i] = i
    
    #print('ORDEM RELEVANCIA',ordem_relevancia)
    #simulando os atendimentos
    saidas = [0]*n_atendentes
    criterio_parar = 0
    parada_1 = 0
    n_contagem = 0
    
    while parada_1 == 0:
        n_contagem = n_contagem+1
        liberacao = min(saidas)
        atendente = saidas.index(liberacao)
        #print('liberacao',liberacao,atendente)
        
        #identifico os primeiros de cada grupo
        primeiros = ()
        chegadas_em_fila = ()
        for i in range(0,n_grupos):
            chegada_grupo_i = grupos_tempo_chegada[i][0]
            primeiros = np.append(primeiros,chegada_grupo_i)
            if chegada_grupo_i <= liberacao: #pega fila, já vou marcar a hora que a demanda seria deslocada para o grupo ZERO
                chegadas_em_fila = np.append(chegadas_em_fila, chegada_grupo_i + T_zero[i])
                #chegadas_em_fila = np.append(chegadas_em_fila, chegada_grupo_i + T_zero[i])
            else:
                chegadas_em_fila = np.append(chegadas_em_fila, 100*janela_obs)
            
        primeiros = primeiros.tolist()
        chegadas_em_fila = chegadas_em_fila.tolist()
        t_minimo = min(primeiros)
        tz_minimo = min(chegadas_em_fila)
        #print('tz_min',tz_minimo)
        #print('primeiros_1',primeiros)
        #print('chegadas_em_fila',chegadas_em_fila)
        
        if t_minimo > liberacao: #o cara é atendido na hora que chega
            
            #print('primeiros', primeiros)
            grupo_atendido = primeiros.index(t_minimo) + 1
            demanda_atendida = grupos_ordem_chegada[grupo_atendido-1][0]
            demanda_atendida = int(demanda_atendida)
            #print('grupo',grupo_atendido,'demanda',demanda_atendida)
            tempo_na_saida = matriz_dados[1][demanda_atendida-1] + matriz_dados[7][demanda_atendida-1]
            saidas[atendente] = tempo_na_saida
            matriz_dados[3][demanda_atendida-1] = tempo_na_saida
            grupos_ordem_chegada[grupo_atendido-1] = np.delete(grupos_ordem_chegada[grupo_atendido-1],0)
            grupos_tempo_chegada[grupo_atendido-1] = np.delete(grupos_tempo_chegada[grupo_atendido-1],0)
            if len(grupos_ordem_chegada[grupo_atendido-1]) == 0:
                #grupos_ordem_chegada[grupo_atendido-1] = [100*janela_obs]
                grupos_tempo_chegada[grupo_atendido-1] = [100*janela_obs]
                criterio_parar = criterio_parar+1
            
            #print(n_contagem,'CASO 1', demanda_atendida)
            
        if tz_minimo < liberacao: #o cara com tz_minimo é atendido, pois chegou primeiro no grupo z
            
            grupo_atendido = chegadas_em_fila.index(tz_minimo) + 1
            #print('cef', chegadas_em_fila)
            #print('grupo no caso 2', grupo_atendido)
            demanda_atendida = grupos_ordem_chegada[grupo_atendido-1][0]
            demanda_atendida = int(demanda_atendida)
            tempo_na_saida = liberacao + matriz_dados[7][demanda_atendida-1]
            saidas[atendente] = tempo_na_saida
            matriz_dados[3][demanda_atendida-1] = tempo_na_saida
            grupos_ordem_chegada[grupo_atendido-1] = np.delete(grupos_ordem_chegada[grupo_atendido-1],0)
            grupos_tempo_chegada[grupo_atendido-1] = np.delete(grupos_tempo_chegada[grupo_atendido-1],0)
            if len(grupos_ordem_chegada[grupo_atendido-1]) == 0:
                #grupos_ordem_chegada[grupo_atendido-1] = [100*janela_obs]
                grupos_tempo_chegada[grupo_atendido-1] = [100*janela_obs]
                criterio_parar = criterio_parar+1
            
            #print(n_contagem,'CASO 2', demanda_atendida)
            
        if (t_minimo <= liberacao) and (tz_minimo >=liberacao): #vai ser atendido aquele que está no grupo prioritário
            #buscar o grupo prioritário com demanda presente
            parada_2 = 0
            jj = -1 #indice para busca
            while parada_2 == 0:
                jj = jj+1
                grupo_consultado = ordem_relevancia[1][jj]
                grupo_consultado = int(grupo_consultado)
                tempo_t = grupos_tempo_chegada[grupo_consultado-1][0]
                #print('tempo_t',tempo_t)
                if tempo_t <= liberacao: # esse cara é atendido
                    demanda_atendida = grupos_ordem_chegada[grupo_consultado-1][0]
                    #print('demanda atendida', demanda_atendida)
                    demanda_atendida = int(demanda_atendida)
                    tempo_na_saida = liberacao + matriz_dados[7][demanda_atendida-1]        
                    saidas[atendente] = tempo_na_saida
                    matriz_dados[3][demanda_atendida-1] = tempo_na_saida
                    #print('qual grupo', grupo_consultado)
                    #print('grupos antes', grupos_ordem_chegada)
                    grupos_ordem_chegada[grupo_consultado-1] = np.delete(grupos_ordem_chegada[grupo_consultado-1],0)
                    grupos_tempo_chegada[grupo_consultado-1] = np.delete(grupos_tempo_chegada[grupo_consultado-1],0)
                    #print('grupos depois', grupos_ordem_chegada,'len',len(grupos_ordem_chegada[grupo_consultado-1]))
                    if len(grupos_ordem_chegada[grupo_consultado-1]) == 0:
                        #print('VI O IF')
                        #grupos_ordem_chegada[grupo_consultado-1] = [100*janela_obs]
                        grupos_tempo_chegada[grupo_consultado-1] = [100*janela_obs]
                        #print(grupos_tempo_chegada)
                        criterio_parar = criterio_parar+1    
                    parada_2 = 1
                
            #print('grupo AGORA', grupos_tempo_chegada)
            #print(n_contagem,'CASO 3', demanda_atendida)
        
        #print('chegadas_em_fila FINAL', chegadas_em_fila)
        
        if n_contagem == num_chegadas:
            parada_1 = 1
    
    n_1 = np.zeros(n_grupos) #registrar numero de demandas de cada grupo
    n_2 = np.zeros(n_grupos) #registrar numero de demandas de cada grupo atendido dentro do limite de tempo
    n_3 = np.zeros(n_grupos) #vai somar os tempos de espera + atendimento
    nivel_servico = np.zeros(n_grupos)
    tempo_medio = np.zeros(n_grupos)
    
    for i in range(0,num_chegadas):
        grupo_obs = matriz_dados[2][i]
        grupo_obs = int(grupo_obs)
        n_1[grupo_obs-1] = n_1[grupo_obs-1]+1
        n_3[grupo_obs-1] = n_3[grupo_obs-1] + (matriz_dados[3][i] - matriz_dados[1][i])
        if (matriz_dados[3][i] - matriz_dados[1][i]) <= tol[grupo_obs-1]: 
            n_2[grupo_obs-1] = n_2[grupo_obs-1]+1
    
    for i in range(0,n_grupos):
        nivel_servico[i] = n_2[i]/n_1[i]
        tempo_medio[i] = n_3[i]/n_1[i]
        
    #terminei a simulação
    return(nivel_servico, tempo_medio)

#Definir número mínimo de atendentes
n_at_min = 0
medias = np.zeros(n_grupos)
for i in range(0,n_grupos):
    if tipos[i] == 'Exp':
        media = parametros[i][0]
    if tipos[i] == 'Norm':
        media = parametros[i][1]
    if tipos[i] == 'Tri':
        media = (parametros[i][3]+parametros[i][4]+parametros[i][5])/3
    n_at_min = n_at_min + (taxas[i]*media)
    medias[i] = media

n_at_min = int(n_at_min) + 1

#Otimização por GA com número fixo de atendentes
pop = 10
pres = 2
n_ger = 2
taxa_mut = 0.3

def fitness(n_att, T_zero_teste, horizonte):
    fila_sim = fila(n_att, T_zero_teste, horizonte)[0]
    tempos_medios = fila_sim[1]
    NSr = fila_sim[0]
    NSd = NS
    valor = 0
    adequacao = 0
    #print(NSr, NSd)
    for i in range(0,n_grupos):
        if NSr[i] >= NSd[i]:
            valor = valor + (NSr[i] - NSd[i])
            adequacao = adequacao + 1
        else:
            valor = valor + 1000*(NSd[i] - NSr[i])
        
    if adequacao == n_grupos:
        valida_sol = 1
    else:
        valida_sol = 0
    
    return NSr, valor, valida_sol, tempos_medios
  
def tempos_transicao(n_att):
    T_zero_testes = np.zeros((pop,2*n_grupos+2))
    
    #primeira geracao
    for i in range(0,pop):
        
        T_zero_obs = np.zeros(n_grupos)
        for j in range(0,n_grupos):
            T_zero_sorteio = rd.uniform(0,tol[j]-medias[j])
            T_zero_obs[j] = T_zero_sorteio
            T_zero_testes[i][j] = T_zero_sorteio
        
        if n_att == n_at_min:
            fitness_ = fitness(n_att,T_zero_obs,50)
        else:
            fitness_ = fitness(n_att,T_zero_obs,10)
            
        for j in range(n_grupos, 2*n_grupos):
            T_zero_testes[i][j] = fitness_[0][j-n_grupos] #dados dos níveis de serviço
        
        T_zero_testes[i][2*n_grupos] = fitness_[1] #valor de análise
        
        T_zero_testes[i][2*n_grupos + 1] = fitness_[2] #validador
        
    T_zero_testes = T_zero_testes[T_zero_testes[:,2*n_grupos].argsort()]
    #print('geracao 1')
    print(T_zero_testes)
    
    #novas gerações
    for i in range(0,n_ger):   
        for j in range(pres,pop):
            #gerando novo individuo
            #PAI 1 
            p1 = rd.randint(0,pop-1)
            #PAI 2 
            p2 = rd.randint(0,pop-1)
            #compondo o filho
            n_p1 = rd.randint(0,n_grupos-1) #numero máximo de genes do pai 1
            np1 = 0
            np2 = 0
            npx = 0
            T_zero_obs = np.zeros(n_grupos)
            while npx < n_grupos:
                
                moeda = rd.randint(0,1)
                mutacao = rd.uniform(0,1)
                if ((moeda == 0) and (np1 < n_p1)) or (np2 == n_grupos -1): #pai elite
                    T_zero_filho = T_zero_testes[p1][npx]
                    np1 = np1 + 1
                else:
                    T_zero_filho = T_zero_testes[p2][npx]
                    np2 = np2 + 1
                if mutacao <= taxa_mut:
                    T_zero_filho = rd.uniform(0,tol[npx]-medias[npx])
                T_zero_obs[npx] = T_zero_filho
                T_zero_testes[j][npx] = T_zero_filho
                npx = npx+1
            
            if n_att == n_at_min:
                fitness_ = fitness(n_att,T_zero_obs,50)
            else:
                fitness_ = fitness(n_att,T_zero_obs,10)
                
            for x in range(n_grupos, 2*n_grupos):
                T_zero_testes[j][x] = fitness_[0][x-n_grupos] #dados dos níveis de serviço
                
            T_zero_testes[j][2*n_grupos] = fitness_[1] #valor de análise
        
            T_zero_testes[j][2*n_grupos + 1] = fitness_[2] #validador

        T_zero_testes = T_zero_testes[T_zero_testes[:,2*n_grupos].argsort()]
        #print('geracao', i+2)
        print(T_zero_testes)
        
        print('solucao atual')
        print(n_att, T_zero_testes[0])
        
    return T_zero_testes[0]

#testando o número de atendentes

def recomenda():
    parada_3 = 0
    n_atendentes = n_at_min-1

    while parada_3 == 0:
        n_atendentes = n_atendentes+1
        print('numero de atendentes:', n_atendentes)
    
        resultado_n = tempos_transicao(n_atendentes)
    
        tempos_T = np.zeros(n_grupos)
        for i in range(0,n_grupos):
            tempos_T[i] = resultado_n[i]
    
        fila_sim = fitness(n_atendentes, tempos_T, 10)
        niveis_servico = fila_sim[0]
        tempos_medios = fila_sim[3]
    
        print('proposta para tempos de transição', tempos_T)
        print('níveis de serviço estimados', niveis_servico)
        print('tempos médios', tempos_medios)
        if fila_sim[2] == 0:
            print('solução não adequada, avaliar a possibilidade de incluir mais atendentes')
        else:
            print('solução adequada')
            parada_3 = 1

    return (tempos_T, niveis_servico, tempos_medios)

# Botão para executar o teste
if stlt.button("Obter recomendação"):
    with stlt.spinner("Aguarde, o cálculo pode levar bastante tempo"):
        solucao = recomenda()
        stlt.write(f"Solução: {solucao}")
        #TESTE = fila(n_at, n_grupos, taxas, tipos, parametros, T_zero, tol)
    
    # Exibir os resultados
    #stlt.subheader("Resultados da Simulação")
    #for i in range(0,n_grupos):
        #stlt.write(f"Nível de serviço alcançado para o Grupo {i + 1}: {TESTE[0][i]*100}%")
        #stlt.write(f"Tempo médio de espera mais atendimento para o Grupo {i + 1}: {TESTE[1][i]}") 
    
