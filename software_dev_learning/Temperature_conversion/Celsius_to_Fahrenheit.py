def displayCelsiusToFahren(start, end):
     ''' 
    This function calculates temperature conversion from Celsisus to Fahrenheit 
    '''
    print('\n Degrees',' Degrees')
    print('Celsius', 'Fahrenheit')
    
    for temp in range(start, end + 1):
        converted_temp = (9/5 * temp) + 32
        print(' ', format(temp, '4.1f'), ' ', format(converted_temp, '4.1f'))
