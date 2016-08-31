"""Error condition that will be reported to user with a message and cause
the execution to terminate with an error code.

Attributes:
    message -- explanation of the error
    value -- return value (default -1)
"""
class N88ReportedError (Exception):

    def __init__(self, message, value=-1):
        self.message = message
        self.value = value

    def __str__(self):
        return repr(self.message)
