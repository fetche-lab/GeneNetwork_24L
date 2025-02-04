## Basics on the fundamental concepts on Object Oriented Programming with Python 
### The following provides a quick introduction to OOP principles, which are (absraction, encapsulation, inheritance, and polymorphism)
### For more information and learning, visit the following sites; 
`https://jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/`
`https://developer.mozilla.org/en-US/Learn/Python/Quickly_Learn_Object_Oriented_Programming`
`http://www.tutorialspoint.com/python/python_classes_objects.htm`

### A) Inheritance 
With inheritance, a class (child class) inherits attributes and methods of another class (parent class). 
```python
## Inheritance in OOP 
## GeneNetwork Webservice group 

#parent class 
class GeneNetwork_group: 
    def __init__(self, department, project): 
        self.__department = department 
        self.__project = project 
        
    def __str__(self):
        return f'{self.__project} project is under the {self.__department} department in GN' 

#child class 01 
class GN_backend(GeneNetwork_group):
    
    def backend_desc(self):
        return f'These engineers are among the web service dev ops team to maintain GN databases' 
         
#child class 02 
class data_engineers(GeneNetwork_group):
    def data_desc(self):
        return f'This part of the team, deals with processing and uploading datasets to GN databases' 
    
team1 = GN_backend('{Genetics, Genomics, and Informatics}', '{GN_database modelling and scaling}')
print(team1)
print(team1.backend_desc()) 
print()
team2 = data_engineers('{Genetics, Genomics and Informatics}', '{GN data wrangling}')
print(team2)
print(team2.data_desc())

```

### B) Encapsulation 
With encapsulation, the data in classes is protected or hidden from outside access. This involves keeping some attributes private, usually by applying the double underscore before them.
```python
## Encapsulation in OOP 
#parent class 
class GeneNetwork_group: 
    def __init__(self, department, project, project_ID): 
        self.department = department # public variable 
        self.project = project # public variable
        self.__project_ID = project_ID ## private variable 
        
    def __str__(self):
        return f'{self.project} project, with ID number {self.__project_ID}, is under the {self.department} department in GN' 

#child class 01 
class GN_backend(GeneNetwork_group):
    
    def backend_desc(self):
        return f'These engineers are among the web service dev ops team to maintain GN databases' 
         
#child class 02 
class data_engineers(GeneNetwork_group):
    def data_desc(self):
        return f'This part of the team, deals with processing and uploading datasets to GN databases' 
    
team1 = GN_backend('{Genetics, Genomics, and Informatics}', '{GN_database modelling and scaling}', 101)
print(team1)
print(team1.backend_desc()) 
print()
team2 = data_engineers('{Genetics, Genomics and Informatics}', '{GN data wrangling}', 102)
print(team2)
print(team2.data_desc())

## in encapsulation, one can access the protected variable directly using name mangling as follows; 
print(team1._GeneNetwork_group__project_ID)
print(team2._GeneNetwork_group__project_ID) 

```

### C) Polymorphism 
With polymorphism, different classes use the same method/function but with different functionalities, with respect to each type of class 
```python
## Polymorphism in OOP 
## same functions but having different meaning in each of the available classes 
class data_engineers: 
    def description(self):
        return f' This describes a team of engineers working with data wrangling for GN' 
    
class software_engineers: 
    def description(self):
        return f' This describes a team of engineers working on developing softwares for data wrangling in GN' 
    
data = data_engineers()
software = software_engineers()

for information in (data, software):
    print(information.description()) 

```

### D) Abstraction 
With abstraction, we tend to hide the internal details or implementations of a function, and only show its functionalities. 
Encapsulation and abstraction can sometimes be confused. In short, abstraction provides a blueprint of what is to be encapsuled. 
```python
## Abstraction in OOP 
## abstraction is all about hiding internal details or implementations of a function and only to show its fucntionalities 
## Any class with at least one abstraction function is an abstract class 
## The ABC (Abstract Base Class) from the abc module is used to create abstract classes 
from abc import ABC, abstractmethod 

#parent class 
class GN_data_mode:
    def __init__(self, category):
        self.category = category 
    
    def description(self):
        return f'There are different data formats to be uploaded in GN, and this is one of them' 
    
    @abstractmethod 
    def curation_threshold(self, x):
        pass 
    
#child class 01 
class mode1(GN_data_mode):
    def curation_threshold(self, x):
        return f'{x} is the curation threshold for {self.category} data format' 
    
#child class 02
class mode2(GN_data_mode):
    def curation_threshold(self, x):
        return f'{x} is the curation threshold for {self.category} data format' 
    
data1 = mode1('Expression dataset') 
print(data1.description())
print(data1.curation_threshold(0.5))
print()
data2 = mode2('Rqtl2 data bundle') 
print(data2.description())
print(data2.curation_threshold(2.5))   

```
