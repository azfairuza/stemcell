"""contain: get_time"""

from datetime import datetime


def time_format(time: datetime):
    """function to get a formated time time"""
    if isinstance(time, datetime):
        year = time.strftime("%Y")
        month = time.strftime("%m")
        day = time.strftime("%d")
        hour = time.strftime("%H")
        minute = time.strftime("%M")
        return f"{year}-{month}-{day}-{hour}{minute}"
    return "debug"
