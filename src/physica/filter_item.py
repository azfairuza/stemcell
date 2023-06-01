"""contain: filter_item"""


def filter_item(input_lst: list):
    """procedure to remove empty item in a list and remove the '\n' character."""
    new_lst = []
    for item in input_lst:
        new_item = item.rstrip()
        if new_item:
            new_lst.append(new_item)
    return new_lst
