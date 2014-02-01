"""Performs statistics operations for ``mapmuts`` package.

Written by Jesse Bloom.


List of functions
----------------------
`Median` : returns the median of a list of numbers.

`Mean` : returns the mean of a list of numbers.


Documentation for individual functions
----------------------------------------
Documentation for individual functions is provided in their
definitions below.

"""


def Median(xlist):
    """Returns the median of a list of numbers.
   
    `xlist` : a list of numbers for which to calculate the median. Note
    that the entries in this list will be sorted upon completion of
    this function.

    Examples:

    >>> xlist = [1, 2, 3, 4.1]
    >>> print "%.3f" % Median(xlist)
    2.500
    """
    if not xlist:
        raise ValueError, "Empty list."
    xlist.sort()
    n = len(xlist)
    if n % 2 == 0: # even length list, average two middle entries
        med = xlist[n / 2] + xlist[n / 2 - 1]
        med = float(med) / 2.0
        return med
    else: # odd length list, get middle entry
        return xlist[n / 2]


def Mean(numlist):
    """Returns the mean of a list of numbers.
    
    `numlist` : a list of numbers.
    
    Example:
    
    >>> print "%.3f" % Mean([1.0, 1.5, 3])
    1.833
    """
    if not numlist:
        raise ValueError("Empty list.")
    mean = 0.0
    for x in numlist:
        mean += x
    return mean / float(len(numlist))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
