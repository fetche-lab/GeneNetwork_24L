def displayFahrenToCelsius(start, end):
    ''' 
    This function calculates temperature conversion from Fahrenheit to Celsisus 
    '''
    print('\n Degrees',' Degrees')
    print('Fahrenheit', 'Celsius')
    
    for temp in range(start, end + 1):
        converted_temp = (temp -32) * 5/9
        print(' ', format(temp, '4.1f'), ' ', format(converted_temp, '4.1f'))
        