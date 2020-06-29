"""Create unique identifiers for PlatePlan objects."""

import datetime
import time
import hashlib


ID_SEP = "$"


def unique_id_from_time(length=None, salt=None):
    """Create a hex hash string from the current time.
    
    Return a string representing the SHA1 hex digest of the current time.
    Current time is precise to the microsecond, so rapid calls to this
    function are not necessarily unique. For truly unique IDs, use the
    `unique_id` function.
    
    length: int, default=None
        Number of characters of the ID to return. If None, the entire 40
        character ID is returned.
    salt: any object with a str method, default=None
        A salt to add before hashing.
    """
    time_str = datetime.datetime.now().isoformat(" ")
    if salt is not None:
        time_str += str(salt)
    id_ = hashlib.sha1(time_str.encode("utf-8")).hexdigest()
    if length is not None:
        id_ = id_[:length]
    return id_


def unique_id(prefix=""):
    """Return a unique DeepPhenotyping ID.
    
    IDs have the format 
    
        <prefix>$<hex>
        
    where <prefix> is an (optional) user-supplied prefix and <hex> is a
    unique 40 character hex string. The prefix is intended to be interpretable
    by humans while the hex string guarantees uniqueness.
    
    Example
    -------
    >>> unique_id("Plate 1")
    'Plate 1$df4efb10bcf0e53704b2d35254e7b9fa150b6c11'
    >>> unique_id()
    '$de71863bee378ecfb831a91c3884ca28b19327e3'
    """
    delta = 0.000001
    time.sleep(delta)
    id_ = unique_id_from_time()
    id2 = id_
    while id2 == id_:
        time.sleep(delta)
        id2 = unique_id_from_time()
    return prefix + ID_SEP + id_


def prefix(id_):
    """Return the (non-unique) prefix of an ID."""
    return id_.partition(ID_SEP)[0]


def trim_id(id_, length=4, dots=False):
    """Trim the hex string of an ID for display.
    
    If dots=True, adds an ellipsis at the end."""
    parts = id_.partition(ID_SEP)
    if dots and parts[2]:
        suffix = "..."
    else:
        suffix = ""
    return parts[0] + parts[1] + parts[2][:length] + suffix


if __name__ == "__main__":
    print(unique_id_from_time())
    print(unique_id_from_time())
    print(unique_id_from_time(length=8, salt="test"))

    print(unique_id())
    id_ = unique_id("Plate 1")
    print(id_)
    print(prefix(id_))
    print(trim_id(id_))

