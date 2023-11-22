import pandas as pd
import gurobipy as gp
from gurobipy import  GRB
import numpy as np
from src.utils import operacoes_dfs_dicts, gera_df_solution


def core(
        file,
        up_consec = True, # Carregamentos nas UPs ocorre na forma de uma série consecutiva em dias 
        fz_consec = True, # Carregamentos nas fazendas ocorre na forma de uma série consecutiva em dias
        transp_fazenda_unico = True, # Cada transportador só pode estar em uma fazenda a cada momento
        up_transp_unico = True, # Cada UP é carregada por somente uma transportadora
        gruas_max = True, # Limita a quantidade de frentes de acordo com o número de gruas máximo
        rsp = True, # Limita o intervalo de RSP médio possível diário
        demanda_fabrica = True, # Volume total diário limitado pela fábrica
        volume_up_min = True, # Volume total carregado de cada UP maior igual a uma porcentagem do volume da UP (caso a UP seja carregada)
        prct_min_transp = True, # Cada frente de trabalho deve ter uma percentual mínimo de caminhões definido pela transportadora    
        frota = True # Soma de caminhões limitadas pela frota_min e frota_max em cada transportadora, em cada dia
    ):

    # Carregando inputs
    df_horizonte = pd.read_excel(file, sheet_name= 'HORIZONTE')
    df_frota = pd.read_excel(file, sheet_name= 'FROTA')
    df_grua = pd.read_excel(file, sheet_name= 'GRUA')
    df_bd_up = pd.read_excel(file, sheet_name= 'BD_UP')
    df_fabrica = pd.read_excel(file, sheet_name= 'FABRICA')
    df_rota = pd.read_excel(file, sheet_name= 'ROTA')

    # Parâmetros do modelo
    big_M = 100000000
    eps = 0.01
    prct_min_vol_up = min(.95, df_fabrica['DEMANDA_MIN'].sum()/df_bd_up['VOLUME'].sum()) # Porcentagem mínima do volume de uma UP que deve ser carregado, caso essa UP seja carregada
    
    qtd_fazendas = df_bd_up['FAZENDA'].nunique()
    qtd_ups = len(df_bd_up['UP'])
    qtd_transportadoras = df_frota['TRANSPORTADOR'].nunique()
    qtd_dias = len(df_horizonte)
    
    fazendas = range(qtd_fazendas)
    ups = range(qtd_ups)
    transportadoras = range(qtd_transportadoras)
    dias = range(qtd_dias)
    
    min_caminhoes = df_frota['FROTA_MIN'].min()*df_grua['PORCENTAGEM_VEICULOS_MIN'].min()
    max_caminhoes = df_frota['FROTA_MAX'].max()
    df_rota['qtd_ciclos'] = df_rota['TEMPO_CICLO'].apply(lambda x: np.floor(x)) # Simplificação da variável TEMPO_CICLO

    # Operações nos dataframes para criar variáveis que terão utilidade na construção das variáveis e restrições
    df_bd_up, df_rota, dict_tr_up, dict_uptr_crgcic,\
        dict_fz_ups, arr_min_frota, arr_max_frota  = operacoes_dfs_dicts(df_bd_up, df_grua, df_rota, 
                                                                         df_frota, ups, transportadoras, 
                                                                         qtd_dias)

    # Criação do modelo
    m = gp.Model('model')

    # Criação das variáveis dos modelos

    # Variável de decisão principal que, para cada transportadora, define quantos caminhões terão em cada UP em cada dia
    caminhoes = m.addVars( 
    qtd_ups,
    qtd_dias,
    qtd_transportadoras,
    vtype=GRB.SEMIINT,
    name="CAMINHOES",
    lb = arr_min_frota,
    ub = arr_max_frota) 

    # x: Variáveis de controle que dizem se haverá caminhoes na UP(fazenda) ou não
    # z: Variáveis de controle que indicam se haverá caminhoes na UP(fazenda) no dia atual mas que não havia no dia anterior
    # Ex.:
    # x = [[0, 0, 0, 1, 1, 1, 0, 0], 
    #      [0, 1, 1, 1, 1, 1, 0, 0]]
    # z = [[0, 0, 0, 1, 0, 0, 0, 0],
    #      [0, 1, 0, 0, 0, 0, 0, 0]]
    x = m.addVars(
        qtd_ups,
        qtd_dias,
        qtd_transportadoras,
        vtype=GRB.BINARY,
        name="CONTROLE_CAMINHOES_UP"
    ) 
    
    z = m.addVars(
        qtd_ups,
        qtd_dias,
        qtd_transportadoras,
        vtype=GRB.BINARY,
        name="CONTROLE_CAMINHOES_UP_D-1"
    ) 
    
    x_fz = m.addVars(
        qtd_fazendas,
        qtd_dias,
        qtd_transportadoras,
        vtype=GRB.BINARY,
        name="CONTROLE_CAMINHOES_FZ"
    ) 
    
    z_fz = m.addVars(
        qtd_fazendas,
        qtd_dias,
        qtd_transportadoras,
        vtype=GRB.BINARY,
        name="CONTROLE_CAMINHOES_FZ_D-1"
    ) 
    
    # Variavel auxiliar para executar uma cláusula if-else na restrição de que as fazendas tem de ter série consecutiva de carregamentos
    b = m.addVars(
        qtd_fazendas,
        qtd_dias,
        qtd_transportadoras,
        vtype=GRB.BINARY,
        name="CONTROLE_CAMINHOES_FZ_AUX"
    ) 
    
    # Variáveis auxiliar da restrição de RSP, utilizada para executar uma operação de divisão
    vol_tot_inv = m.addVars(
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="1/VOLUME_TOTAL_DIARIOS"
    )

    # Variáveis auxiliares para executar o cálculo da função objetivo
    bd_max = m.addVars(
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="BD_MAX"
    )

    bd_min = m.addVars(
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="BD_MIN"
    )

    bd_max_aux = m.addVars(
        qtd_ups,
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="BD_DIA_CALC_MAX"
    )

    bd_min_aux = m.addVars(
        qtd_ups,
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="BD_DIA_CALC_MIN"
    )

    var_bd_aux = m.addVars(
        qtd_ups,
        qtd_dias,
        vtype=GRB.BINARY,
        name="VAR_BD_AUX"
    )

    var_bd = m.addVars(
        qtd_dias,
        vtype=GRB.CONTINUOUS,
        name="VAR_BD"
    )

    # Adiciona as restrições ao modelo.     
    params = (
        m, caminhoes, x, z, x_fz, z_fz, b, vol_tot_inv, ups, dias, transportadoras, fazendas,
        df_bd_up, df_grua, df_fabrica, df_frota, dict_fz_ups, dict_tr_up, dict_uptr_crgcic,
        max_caminhoes, min_caminhoes, prct_min_vol_up, eps, big_M, up_consec, fz_consec, 
        transp_fazenda_unico, up_transp_unico,  gruas_max, rsp, demanda_fabrica, 
        volume_up_min, prct_min_transp, frota
    )
    m = adiciona_restricoes(*params)

    # Adiciona a função objetivo ao modelo
    params = (
        m, bd_max, bd_min, bd_max_aux, bd_min_aux, var_bd, var_bd_aux, x,
        df_bd_up, transportadoras, ups, dias, big_M, eps
    )
    m = set_f_obj(*params)

    # Resolvendo o modelo
    m.update()
    m.printStats()
    #m.Params.NonConvex = 2
    m.optimize()
    status = m.Status
    if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
        print('O modelo não pode ser resolvido')
    if status != GRB.INF_OR_UNBD and status != GRB.INFEASIBLE:
        print('A otimização terminou com status %d' % status)
    if status == GRB.OPTIMAL:
        print('\nO valor ótimo da função objetivo é %g' % m.ObjVal)
        df_solution = gera_df_solution(caminhoes, df_bd_up, df_grua, df_horizonte, df_rota)
        file = 'solution_output.xlsx'
        df_solution.to_excel(file)
        print(f'\nArquivo {file} gerado com sucesso')
    

def adiciona_restricoes(
        m, caminhoes, x, z, x_fz, z_fz, b, vol_tot_inv, ups, dias, transportadoras, fazendas,
        df_bd_up, df_grua, df_fabrica, df_frota, dict_fz_ups, dict_tr_up, dict_uptr_crgcic,
        max_caminhoes, min_caminhoes, prct_min_vol_up, eps, big_M, up_consec, fz_consec,
        transp_fazenda_unico, up_transp_unico, gruas_max, rsp, demanda_fabrica,
        volume_up_min,prct_min_transp, frota
    ):
    
    # A variável de caminhoes só será diferente de zero quando a variável de controle for igual a 1
    # e deve ser nula quando a variável de controle for igual a 0
    m.addConstrs(caminhoes[up, d, tr] <= x[up, d, tr] * max_caminhoes for up in ups for d in dias for tr in transportadoras)
    m.addConstrs(caminhoes[up, d, tr] >= x[up, d, tr] * min_caminhoes for up in ups for d in dias for tr in transportadoras)

    # Para cada fazenda, para cada dia, para cada transportadora, se alguma UP está sendo carregada (x[up,d,tr]==1) então z[fz,d,tr] == 1
    # Modelando x > 0 --> b = 1 else b = 1 
    #           b = 1 --> z = 1 else z = 0    
    m.addConstrs(gp.quicksum(x[up,d,tr] for up in dict_fz_ups[fz]) >= 0 + eps- big_M* (1 - b[fz,d,tr]) for d in dias for fz in fazendas for tr in transportadoras)
    m.addConstrs(gp.quicksum(x[up,d,tr] for up in dict_fz_ups[fz]) <= 0 + big_M* b[fz,d,tr] for d in dias for fz in fazendas for tr in transportadoras)
    m.addConstrs((b[fz,d,tr] == 1) >> (x_fz[fz,d,tr] == 1) for d in dias for fz in fazendas for tr in transportadoras)
    m.addConstrs((b[fz,d,tr] == 0) >> (x_fz[fz,d,tr] == 0) for d in dias for fz in fazendas for tr in transportadoras)

    # Transportadoras especificas podem ser fazer o serviço em cada fazenda, de acordo com o df_rota
    print(dict_tr_up)
    for up in ups:
        for tr in set(transportadoras).difference(dict_tr_up[up]):
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] for d in dias) == 0)
    
    # Cada UP tem de ter somente uma série consecutiva de carregamentos
    # z[up,d, tr] == 1 se há caminhoes da transportadora tr na up no dia d, mas não no dia d-1
    if up_consec:
        m.addConstrs((x[up, d2, tr] - x[up, d1, tr] <= z[up, d2, tr]) for up in ups for d1,d2 in zip(dias, dias[1:]) for tr in transportadoras)
        m.addConstrs(x[up, 0, tr] == z[up, 0, tr] for up in ups for tr in transportadoras)
        m.addConstrs(z[up, d, tr] <= x[up, d, tr] for d in dias for up in ups for tr in transportadoras)
        m.addConstrs(gp.quicksum(z[up,d,tr] for d in dias) <= 1 for up in ups for tr in transportadoras)
    
    # Cada fazenda tem de ter somente uma série consecutiva de carregamentos
    # z[fz,d,tr] == 1 se há caminhoes da transportadora tr na fazenda no dia d, mas não no dia d-1
    if fz_consec:
        m.addConstrs((x_fz[fz, d2, tr] - x_fz[fz, d1, tr] <= z_fz[fz, d2, tr]) for fz in fazendas for d1,d2 in zip(dias, dias[1:]) for tr in transportadoras)
        m.addConstrs(x_fz[fz, 0, tr] == z[fz, 0, tr] for fz in fazendas for tr in transportadoras)
        m.addConstrs(gp.quicksum(z_fz[fz,d,tr] for d in dias) <= 1 for fz in fazendas for tr in transportadoras)
        m.addConstrs(z_fz[fz, d, tr] <= x_fz[fz, d, tr] for d in dias for fz in fazendas for tr in transportadoras)

    # Em cada dia, se caminhoes[UP, dia] != 0, caminhoes[UP, dia] >= %min[transportadora] * (sum_{UP} caminhoes[UP, dia] tq UP IN transportadora)
    if prct_min_transp:
        for d in dias:
            for tr in transportadoras:
                m.addConstrs(caminhoes[up, d, tr] >= x[up, d, tr] * df_grua.loc[tr,'PORCENTAGEM_VEICULOS_MIN']*gp.quicksum(caminhoes[up2, d, tr] 
                                                                                                                        for up2 in ups) for up in ups)

    # Cada transportador só pode estar em uma fazenda a cada momento
    if transp_fazenda_unico:
        m.addConstrs(gp.quicksum(x_fz[fz,d,tr] for fz in fazendas) <= 1 for tr in transportadoras for d in dias)

    # Cada up é carregada por somente uma transportadora
    if up_transp_unico:
        m.addConstrs(gp.quicksum(x[up,d,tr] for tr in transportadoras ) <= 1 for up in ups for d in dias)
    
    # Em cada transportadora, em cada dia, sum_{up} x[up,dia] < QTD_gruas[transportadora]
    if gruas_max:
        m.addConstrs(gp.quicksum(x[up,d,tr] for up in ups) <= df_grua.loc[tr,'QTD_GRUAS'] for d in dias for tr in transportadoras)

    # Em cada dia, sum_{up, tr} caminhoes[up, dia, tr]*capacidade[up,transportadora]*qtd_ciclos[up,transportadora] >= demanda_fabrica_min[dia]
    if demanda_fabrica:
        for d in dias:
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] 
                                    for up in ups for tr in transportadoras) >= df_fabrica.loc[d,'DEMANDA_MIN'])
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] 
                                    for up in ups for tr in transportadoras) <= df_fabrica.loc[d,'DEMANDA_MAX'])
            
    # Em cada UP, sum_{dias} caminhoes[up,dias,transportadoras]*qtd_ciclos[up,transporadora]*capacidade[up,transportadora] <= VOLUME_TOTAL_UP    
    for up in ups:
        m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] 
                                for d in dias for tr in transportadoras) <= df_bd_up.loc[up,'VOLUME'])
    
    # Em cada UP que tenha carregamento, o volume total carregado tem de ser pelo menos prct_min_vol_up*VOLUME_UP
    if volume_up_min:
        for up in ups:    
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] 
                                    for d in dias for tr in transportadoras) >= gp.quicksum(x[up,d,tr] for d in dias for tr in transportadoras) * prct_min_vol_up*df_bd_up.loc[up,'VOLUME'])

    # Em cada dia, a média dos RSPs de cada UP, ponderada pelos volumes transportados, deve estar dentro dos limites dados pela Fábrica.
    # volume[up,d] = sum_{tr} caminhoes[up,d,tr]*capacidade[up,tr]*ciclos[up,tr]
    # sum_{up} volume[up,d]*RSP[up]/volume_total >= limite_min_rsp[d]
    # sum_{up} volume[up,d]*RSP[up]/volume_total <= limite_max_rsp[d]
    if rsp:
        for d in dias:
            volume_tot_dia = gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] for tr in transportadoras for up in ups)
            m.addConstr(vol_tot_inv[d] * volume_tot_dia == 1)
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] * df_bd_up.loc[up, 'RSP']
                                    for tr in transportadoras for up in ups) * vol_tot_inv[d] >= df_fabrica.loc[d,'RSP_MIN'])
            m.addConstr(gp.quicksum(caminhoes[up, d, tr] * dict_uptr_crgcic[f'{up}{tr}'][1] * dict_uptr_crgcic[f'{up}{tr}'][0] * df_bd_up.loc[up, 'RSP']
                                    for tr in transportadoras for up in ups) * vol_tot_inv[d] <= df_fabrica.loc[d,'RSP_MAX'])

    # Em cada dia, em cada transportadora, a quantidade de caminhões alocados nas UPs deve estar limitado aos valores frota_min e frota_max
    if frota:
        m.addConstrs(gp.quicksum(caminhoes[up,d,tr] for up in ups) >= df_frota.loc[tr,'FROTA_MIN'] for tr in transportadoras for d in dias)
        m.addConstrs(gp.quicksum(caminhoes[up,d,tr] for up in ups) <= df_frota.loc[tr,'FROTA_MAX'] for tr in transportadoras for d in dias)
    
    return m


def set_f_obj(
        m, bd_max, bd_min, bd_max_aux, bd_min_aux, var_bd, var_bd_aux, x,
        df_bd_up, transportadoras, ups, dias, big_M, eps):
    # Restrições para calcular o máximo e o mínimo do DB diário. 
    # em cada dia, em cada UP, bd_max_aux recebe o BD caso tenha havido algum carregamento naquela UP naquele dia e recebe 0 caso não tenha havido carregamento
    # Posteriormente é calculado o máximo em cada dia e armazenado em db_max
    # em cada dia, em cada UP, bd_min_aux recebe o BD caso tenha havido algum carregamento naquela UP naquele dia e recebe big_M caso não tenha havido carregamento
    # Posteriormente é calculado o mínimo em cada dia e armazenado em db_min
    for d in dias:
        for up in ups:
            m.addConstr(bd_max_aux[up,d] == gp.quicksum(x[up,d,tr] for tr in transportadoras)*df_bd_up.loc[up, 'DB'])
            m.addConstr(gp.quicksum(x[up,d,tr] for tr in transportadoras)*df_bd_up.loc[up, 'DB'] >= 0 + eps - big_M* (1 - var_bd_aux[up,d]))
            m.addConstr(gp.quicksum(x[up,d,tr] for tr in transportadoras)*df_bd_up.loc[up, 'DB'] <= 0 + big_M* var_bd_aux[up,d])
            m.addConstr((var_bd_aux[up,d] == 1) >> (bd_min_aux[up,d] == gp.quicksum(x[up,d,tr] for tr in transportadoras)*df_bd_up.loc[up, 'DB']))
            m.addConstr((var_bd_aux[up,d] == 0) >> (bd_min_aux[up,d] == big_M))
        m.addConstr(bd_max[d] == gp.max_([bd_max_aux[up,d] for up in ups]))
        m.addConstr(bd_min[d] == gp.min_([bd_min_aux[up,d] for up in ups]))
    m.addConstrs((var_bd[d] == bd_max[d] - bd_min[d]) for d in dias)
    
    # Criação da função objetivo
    f_obj = gp.quicksum(var_bd[d]/len(dias) for d in dias)
    m.setObjective(f_obj, GRB.MINIMIZE)
    return m