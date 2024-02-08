# ---- main 

# Display program welcome 
displayWelcome()

# Get which conversion from user
which = getConvertTo()

# Get range of temperatures to convert 
temp_start = int(input('Enter starting temperature to convert: '))
temp_end = int(input('Enter ending temperature to convert: '))

# Display range of converted temperatures
if which == 'F':
    displayFahrenToCelsius(temp_start, temp_end)
else:
    displayCelsiusToFahren(temp_start, temp_end)
    