"""Shared utilities and helper functions."""

import collections


def partition(n, k):
    """Split an integer n into k groups as evenly as possible.
    
    Example
    -------
    >>> partition(13, 5)
    [3, 3, 3, 2, 2]
    """
    sizes = [n // k for _ in range(k)]
    for i in range(n % k):
        sizes[i] += 1
    return sizes


def distribute(items, max_size, even=False):
    """Split items into groups of at most max_size items.
    
    Inputs
    ------
    items: list
        List (or list-like object) of items to distribute.
    max_size: int
        Maximum number of elements in each group.
    even: bool, default=False
        If True, distribute the items evenly among the groups. If False,
        all groups except the final group will have `max_size` elements.
        
    Returns
    -------
    A list of lists with at most `max_size` elements. Only the final group
    will have fewer than `max_size` elements. If `items` is empty an empty
    list is returned.
    
    Example
    -------
    >>> distribute([1,2,3,4,5,6,7], 3)
    [[1, 2, 3], [4, 5, 6], [7]]
    
    >>> distribute([1,2,3,4,5,6,7], 3, even=True)
    [[1, 2, 3], [4, 5], [6, 7]]
    """
    if not items:
        return []

    n = len(items)
    k = n // max_size
    groups = []

    if even:
        if k * max_size < n:
            k += 1
        sizes = partition(n, k)
        idx = 0
        for size in sizes:
            groups.append(items[idx : idx + size])
            idx += size
    else:
        for i in range(k):
            groups.append(items[i * max_size : (i + 1) * max_size])
        if n % max_size > 0:
            groups.append(items[k * max_size :])

    return groups


def assert_list(items, none_as_empty=True):
    """Ensure items is a list; if not, make it one.
    
    If `none_as_empty` is True, calling with items=None returns [].
    
    >>> assert_list([1,2])
    [1,2]
    >>> assert_list('abc')
    ['abc']
    >>> assert_list(None)
    []
    >>> assert_list(None, none_as_empty=False)
    [None]
    """
    if not isinstance(items, list):
        if isinstance(items, collections.abc.KeysView):
            return list(items)
        elif items is None and none_as_empty:
            return []
        else:
            return [items]
    else:
        return items


def flatten(list_of_lists):
    """Flatten a list of lists into a single list.
    
    >>> flatten([[2, 3], [4], [5, 6, 7], []])
    [2, 3, 4, 5, 6, 7]
    """
    return [i for sub in list_of_lists for i in sub]


def not_none(l):
    """Filter a list by removing None elements."""
    return [x for x in l if x is not None]
