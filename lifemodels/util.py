def prob_to_code(v):
    '''
    Convert P values to significance codes:
    less than 0.001: ***
    (0.01, 0.001]  : **
    (0.05, 0.01]   : *
    (0.1, 0.05]    : +

    prob_to_code(v)
    v: float
    '''
    if   0.1  >v>=0.05:
        return '+'
    elif 0.05 >v>=0.01:
        return '*'
    elif 0.01 >v>=0.001:
        return '**'
    elif 0.001>v>=0:
        return '***'
    else:
        return ''
