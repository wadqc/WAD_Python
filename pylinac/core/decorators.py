
# The following is adapted from: http://code.activestate.com/recipes/578809-decorator-to-check-method-param-types/
# Another type checking decorator: http://code.activestate.com/recipes/454322-type-checking-decorator/
from abc import ABCMeta
from functools import wraps
# from inspect import signature
import time

import numpy as np


def timethis(func):
    """Report execution time of function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, 'took {:3.2f}s'.format(end-start))
        return result
    return wrapper

def type_accept(*type_args, **type_kwargs):
    """Decorator to check function/method input types. Based on Python Cookbook 3rd ed. #9.7."""
    def decorate(func):

        # Map function argument names to supplied types
        sig = signature(func)
        bound_types = sig.bind_partial(*type_args, **type_kwargs).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound_values = sig.bind(*args, **kwargs)
            # Enforce type assertions across supplied arguments
            for name, value in bound_values.arguments.items():
                if name in bound_types:
                    if type(bound_types[name]) in (type, ABCMeta):  # Single-type comparisons
                        if not isinstance(value, bound_types[name]):
                            raise TypeError("Argument '{}' must be {}".format(name, bound_types[name]))
                    else:
                        if type(value) not in bound_types[name]:
                            if value not in bound_types[name]:
                                raise TypeError("Argument '{}' must be {}".format(name, bound_types[name]))
            return func(*args, **kwargs)
        return wrapper
    return decorate

def value_accept(*value_args, **value_kwargs):
    """Decorator to check function/method input types. Based on Python Cookbook 3rd ed. #9.7."""
    def decorate(func):
        sig = signature(func)
        # convert any dictionary value acceptances to tuples
        vkw = convert_dictvals2tuple(value_kwargs)
        # Map function argument names to supplied types
        bound_values = sig.bind_partial(*value_args, **vkw).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            passed_values = sig.bind(*args, **kwargs)
            # Enforce value assertions across supplied arguments
            for name, value in passed_values.arguments.items():
                if name in bound_values:
                    if type(value) in (float, int, np.float64):
                        # value must be within a number range
                        if not bound_values[name][0] <= value <= bound_values[name][1]:
                            raise ValueError("Argument '{}' needs to be between {:f} and {:f}".format(name,
                                                                                                    bound_values[name][0],
                                                                                                    bound_values[name][1]))
                    else:
                        # value is a str and must be one of the accepted str values
                        if value not in bound_values[name]:
                            raise ValueError("Argument '{}' must be one of {}".format(name, bound_values[name]))
            return func(*args, **kwargs)
        return wrapper
    return decorate

def convert_dictvals2tuple(args):
    """Convert from dictionary to tuple of dictionary values."""
    for arg in args:
        if type(args[arg]) == dict:
            args[arg] = tuple(args[arg].values())
    return args

def unwrap_func(wrapped_func, unwraps=0):
    """Return an unwrapped function from the wrapped function.

    Functions or methods can use one or more decorators. If the decorator needs access to the base function (e.g. its arguments), it
    won't get it if more than 1 decorator is used. This function allows one to get to the base function, or a function wrapped a certain
    number of times. The function must be wrapped with the stdlib functools.wraps for proper unwrapping.

    Discussion of this issue: http://bugs.python.org/issue17482

    :param unwraps: The number of times to unwrap a function, if wrapped with multiple decorators. If 0 or negative, will return the base
        function.
    :type unwraps: int
    """
    func = wrapped_func
    # if we want the base function
    if unwraps <= 0:
        at_base = False
        while not at_base:
            try:
                func = func.__wrapped__
            except AttributeError:
                at_base = True

    # else we want the function unwrapped n times
    else:
        for unwrap in range(unwraps):
            try:
                func = func.__wrapped__
            except AttributeError:
                break

    return func


from collections import namedtuple
from functools import update_wrapper
from threading import RLock

_CacheInfo = namedtuple("CacheInfo", ["hits", "misses", "maxsize", "currsize"])


class _HashedSeq(list):
    __slots__ = 'hashvalue'

    def __init__(self, tup, hash=hash):
        self[:] = tup
        self.hashvalue = hash(tup)

    def __hash__(self):
        return self.hashvalue


def _make_key(args, kwds, typed,
              kwd_mark=(object(),),
              fasttypes={int, str, frozenset, type(None)},
              sorted=sorted, tuple=tuple, type=type, len=len):
    'Make a cache key from optionally typed positional and keyword arguments'
    key = args
    if kwds:
        sorted_items = sorted(kwds.items())
        key += kwd_mark
        for item in sorted_items:
            key += item
    if typed:
        key += tuple(type(v) for v in args)
        if kwds:
            key += tuple(type(v) for k, v in sorted_items)
    elif len(key) == 1 and type(key[0]) in fasttypes:
        return key[0]
    return _HashedSeq(key)


def lru_cache(maxsize=100, typed=False):
    """Least-recently-used cache decorator.

    If *maxsize* is set to None, the LRU features are disabled and the cache
    can grow without bound.

    If *typed* is True, arguments of different types will be cached separately.
    For example, f(3.0) and f(3) will be treated as distinct calls with
    distinct results.

    Arguments to the cached function must be hashable.

    View the cache statistics named tuple (hits, misses, maxsize, currsize) with
    f.cache_info().  Clear the cache and statistics with f.cache_clear().
    Access the underlying function with f.__wrapped__.

    See:  http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

    """

    # Users should only access the lru_cache through its public API:
    #       cache_info, cache_clear, and f.__wrapped__
    # The internals of the lru_cache are encapsulated for thread safety and
    # to allow the implementation to change (including a possible C version).

    def decorating_function(user_function):

        cache = dict()
        stats = [0, 0]  # make statistics updateable non-locally
        HITS, MISSES = 0, 1  # names for the stats fields
        make_key = _make_key
        cache_get = cache.get  # bound method to lookup key or return None
        _len = len  # localize the global len() function
        lock = RLock()  # because linkedlist updates aren't threadsafe
        root = []  # root of the circular doubly linked list
        root[:] = [root, root, None, None]  # initialize by pointing to self
        nonlocal_root = [root]  # make updateable non-locally
        PREV, NEXT, KEY, RESULT = 0, 1, 2, 3  # names for the link fields

        if maxsize == 0:

            def wrapper(*args, **kwds):
                # no caching, just do a statistics update after a successful call
                result = user_function(*args, **kwds)
                stats[MISSES] += 1
                return result

        elif maxsize is None:

            def wrapper(*args, **kwds):
                # simple caching without ordering or size limit
                key = make_key(args, kwds, typed)
                result = cache_get(key, root)  # root used here as a unique not-found sentinel
                if result is not root:
                    stats[HITS] += 1
                    return result
                result = user_function(*args, **kwds)
                cache[key] = result
                stats[MISSES] += 1
                return result

        else:

            def wrapper(*args, **kwds):
                # size limited caching that tracks accesses by recency
                key = make_key(args, kwds, typed) if kwds or typed else args
                with lock:
                    link = cache_get(key)
                    if link is not None:
                        # record recent use of the key by moving it to the front of the list
                        root, = nonlocal_root
                        link_prev, link_next, key, result = link
                        link_prev[NEXT] = link_next
                        link_next[PREV] = link_prev
                        last = root[PREV]
                        last[NEXT] = root[PREV] = link
                        link[PREV] = last
                        link[NEXT] = root
                        stats[HITS] += 1
                        return result
                result = user_function(*args, **kwds)
                with lock:
                    root, = nonlocal_root
                    if key in cache:
                        # getting here means that this same key was added to the
                        # cache while the lock was released.  since the link
                        # update is already done, we need only return the
                        # computed result and update the count of misses.
                        pass
                    elif _len(cache) >= maxsize:
                        # use the old root to store the new key and result
                        oldroot = root
                        oldroot[KEY] = key
                        oldroot[RESULT] = result
                        # empty the oldest link and make it the new root
                        root = nonlocal_root[0] = oldroot[NEXT]
                        oldkey = root[KEY]
                        oldvalue = root[RESULT]
                        root[KEY] = root[RESULT] = None
                        # now update the cache dictionary for the new links
                        del cache[oldkey]
                        cache[key] = oldroot
                    else:
                        # put result in a new link at the front of the list
                        last = root[PREV]
                        link = [last, root, key, result]
                        last[NEXT] = root[PREV] = cache[key] = link
                    stats[MISSES] += 1
                return result

        def cache_info():
            """Report cache statistics"""
            with lock:
                return _CacheInfo(stats[HITS], stats[MISSES], maxsize, len(cache))

        def cache_clear():
            """Clear the cache and cache statistics"""
            with lock:
                cache.clear()
                root = nonlocal_root[0]
                root[:] = [root, root, None, None]
                stats[:] = [0, 0]

        wrapper.__wrapped__ = user_function
        wrapper.cache_info = cache_info
        wrapper.cache_clear = cache_clear
        return update_wrapper(wrapper, user_function)

    return decorating_function
