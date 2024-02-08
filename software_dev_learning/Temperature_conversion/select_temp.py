def getConvertTo():
    ''' 
    This function provides two options to choose from. 
    01: F for Fahrenheit temperature conversion 
    02: C for Celsius temperature conversion 
    '''
    which = input('Enter selection: ')
    while which != 'F' and which != 'C':
        which = input('Enter selection: ')
    return which 