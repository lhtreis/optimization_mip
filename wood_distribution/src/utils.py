import pandas as pd
import numpy as np

def operacoes_dfs_dicts(df_bd_up, df_grua, df_rota, df_frota, ups, transportadoras, qtd_dias):
    # Cria dicionários que associa as fazendas aos respectivos indices
    ar = df_bd_up['FAZENDA'].unique()
    dict_fazendas = dict(zip(range(len(ar)), ar))
    dict_fazendas_inv = {v: k for k, v in dict_fazendas.items()}

    # Cria dicionários que associa as transportadoras aos respectivos indices
    ar = df_grua['TRANSPORTADOR'].unique()
    dict_transportadoras = dict(zip(range(len(ar)), ar))
    dict_transportadoras_inv = {v: k for k, v in dict_transportadoras.items()}

    # Adiciona os informações dos indices no df_rota
    df_rota['Fazenda_idx'] = df_rota['Fazenda'].map(dict_fazendas_inv)
    df_rota['TRANSPORTADOR_idx'] = df_rota['TRANSPORTADOR'].map(dict_transportadoras_inv)
    #dict_fz_tr = df_rota[['Fazenda_idx','TRANSPORTADOR_idx']].drop_duplicates().groupby('Fazenda_idx')['TRANSPORTADOR_idx'].unique().to_dict()

    # Cria dicionário que associa os idx das fazendas aos idx das UPs
    df_bd_up['up_idx'] = df_bd_up.index
    df_bd_up['Fazendas_idx'] = df_bd_up['FAZENDA'].map(dict_fazendas_inv)
    dict_fz_ups = df_bd_up[['up_idx', 'Fazendas_idx']].groupby('Fazendas_idx')['up_idx'].unique().to_dict()
    
    # Adiciona os informações dos indices no db_bd_up
    df_fazendas_idx = pd.DataFrame(list(dict_fazendas_inv.items()), columns=['FAZENDA', 'idx'])
    df_bd_up = df_bd_up.merge(df_fazendas_idx, how = 'left', on='FAZENDA')

    # Cria dicionário que associa os idx das transportadoras aos idx das UPs
    df_rota = df_rota.merge(df_bd_up[['UP', 'up_idx']], how = 'left', left_on='ORIGEM', right_on='UP')
    dict_tr_up = df_rota[['up_idx', 'TRANSPORTADOR_idx']].groupby('up_idx')['TRANSPORTADOR_idx'].unique().to_dict()

    # Cria dicionário que associa o par UP-TR aos respectivos tempos de ciclo e caixa_carga
    dict_uptr_crgcic = {}
    for i, row in df_rota[['up_idx', 'TRANSPORTADOR_idx', 'qtd_ciclos', 'CAIXA_CARGA']].iterrows():
        up = int(row['up_idx'])
        tr = int(row['TRANSPORTADOR_idx'])
        dict_uptr_crgcic[f'{up}{tr}'] = [row['qtd_ciclos'], row['CAIXA_CARGA']]
    for up in ups:
        for tr in transportadoras:
            try:
                dict_uptr_crgcic[f'{up}{tr}']
            except KeyError:
                dict_uptr_crgcic[f'{up}{tr}'] = [0,0]

    # Cria arrays com as quantidades minimas e máximas de caminhoes que 
    arr_min_frota = np.zeros((len(ups),qtd_dias,len(transportadoras)))
    arr_max_frota = np.zeros((len(ups),qtd_dias,len(transportadoras)))
    for tr in transportadoras:
        arr_min_frota[:,:,tr] = np.floor(df_frota.loc[tr,'FROTA_MIN'] * df_grua.loc[tr,'PORCENTAGEM_VEICULOS_MIN'])
        arr_max_frota[:,:,tr] = df_frota.loc[tr,'FROTA_MAX']

    return df_bd_up, df_rota, dict_tr_up, dict_uptr_crgcic, \
            dict_fz_ups, arr_min_frota, arr_max_frota


def gera_df_solution(caminhoes, df_bd_up, df_grua, df_horizonte, df_rota):
    up_sol = []
    tr_sol = []
    dia_sol = []
    cm_sol = []
    features_output =  ['UP', 'FAZENDA', 'TRANSPORTADOR', 'DIA', 'DB', 'RSP', 'QTE_VEICULOS', 'VOLUME'] 

    # Cria o df_solution a partir dos valores na variavel 'caminhoes'
    solution = { k : v.X for k,v in caminhoes.items() if v.X > 0}
    for chave, valor in solution.items():
        up, dia, transportador = chave
        up_sol.append(up)
        dia_sol.append(dia)
        tr_sol.append(transportador)
        cm_sol.append(valor)
    df_solution = pd.DataFrame({
        'UP_IDX': up_sol,
        'DIA_IDX': dia_sol,
        'TRANSPORTADOR_IDX': tr_sol,
        'QTE_VEICULOS': cm_sol
        })
    
    # Merge com os df's de input
    df_solution = df_solution.merge(df_bd_up[['UP', 'FAZENDA', 'DB', 'RSP' ,'up_idx']], how = 'left', left_on= 'UP_IDX', right_on= 'up_idx')
    df_solution = df_solution.merge(df_grua[['TRANSPORTADOR']], how='left', left_on='TRANSPORTADOR_IDX', right_index=True)
    df_solution = df_solution.merge(df_horizonte[['DIA']], how='left', left_on='DIA_IDX', right_index=True)
    df_solution = df_solution.merge(df_rota[['ORIGEM', 'TRANSPORTADOR','qtd_ciclos', 'CAIXA_CARGA' ]], how = 'left', left_on= ['UP','TRANSPORTADOR'], right_on=['ORIGEM', 'TRANSPORTADOR'])
    df_solution['VOLUME'] = df_solution['QTE_VEICULOS'] * df_solution['qtd_ciclos'] * df_solution['CAIXA_CARGA']

    return df_solution[features_output].sort_values(by='DIA').reset_index(drop= True)
