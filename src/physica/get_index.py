def getIndex(val: float, lst: list):
    """procedure to get index position of a list from a value

    it checks whether the value is between `i` and `i+1` element of 
    the list. If it is true, it will return the `i` value.
    
    Parameters
    ----------
    val: float
        the value to be searched
    lst: list
        the ist
    
    Return
    ------
    i:
        index of the lowest boundary

    Note
    ----
    Have not implemented for the case of value outside the list     
    """
    for i in range(len(lst)):
        if lst[i] <= val < lst[i+1]:
            return i
    return "The value is not in the range list"
