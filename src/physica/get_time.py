from datetime import datetime

def timeFormat(time: datetime):
    """function to get a formated time time"""
    if isinstance(time, datetime):
        year = time.strftime('%Y')
        month = time.strftime('%m')
        day = time.strftime('%d')
        hour = time.strftime('%H')
        minute = time.strftime('%M')
        return f'{year}-{month}-{day}-{hour}{minute}'
    else:
        return 'debug'