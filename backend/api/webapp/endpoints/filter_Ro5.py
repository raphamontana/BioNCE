# Filter for Rule of Five
def rule_of_five(values):
    """
    Filter for Rule of Five
    input = list with 8 descriptors from descriptors.py
    output = True if molecule passes all
    """
    # List, each value will correspond to a violation
    violations = [0, 0, 0, 0]
    if values[0] > 5:  # logP (octanol/water partition coef) > 5
        violations[0] = 1
    if values[1] > 500:  # molecular weight > 500
        violations[1] = 1
    if values[2] > 10:  # number HBond acceptors > 10
        violations[2] = 1
    if values[3] > 5:  # number HBond donors > 5
        violations[3] = 1
    if sum(violations):
        return False  # molecule does not pass Ro5 filter
    return True
