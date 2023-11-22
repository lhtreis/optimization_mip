from src.core import core

def main():
    # Os parâmetros opcionais são variáveis booleanas associadas às restrições existentes no modelo.
    # Caso passe algum dos argumentos opcionais como 'False', a restrição em questão não será adicionada ao modelo
        core(
            filename,
            up_consec = False, # Carregamentos nas UPs ocorre na forma de uma série consecutiva em dias 
            fz_consec = False, # Carregamentos nas fazendas ocorre na forma de uma série consecutiva em dias 
            transp_fazenda_unico = False, # Cada transportador só pode estar em uma fazenda a cada momento
            up_transp_unico = False, # Cada UP é carregada por somente uma transportadora em cada dia
            #gruas_max = False, # Limita a quantidade de frentes de acordo com o número de gruas máximo
            rsp = False, # Limita o intervalo de RSP médio possível diário
            #demanda_fabrica = False, # Volume total diário limitado pela fábrica
            volume_up_min = False, # Volume total carregado de cada UP maior igual a uma porcentagem do volume da UP (caso a UP seja carregada)
            prct_min_transp = False, # Cadça frente de trabalho deve ter uma percentual mínimo de caminhões definido pela transportadora
            #frota = False # Soma de caminhões limitadas pela frota_min e frota_max em cada transportadora, em cada dia
        )

if __name__ == "__main__":
    filename = 'generic_input_case.xlsx' 
    main()
    print('\nExecução encerrada!\n')