from dataclasses import fields
from typing import Optional
from venv import logger

from tavi.tavi_model.FileSystem.tavi_class_factory import Scan


class Filter:
    """
    Arg:
        scans: the class that holds the loaded scan data, meta data, ubconf etc. Defined in TaviProject
        operations: a string of filter operations. keyword + operation + value.
                    example: "title+contains+temp", "sample_temp + < + 100". If the operations is: contains, notcontain,
                    then we look in scans.metadata. If the operations is: <, >, <=, >=, ==, !=, then we look in scans.data.
    """

    def __init__(self, scans: dict[str, Scan], operations: Optional[list[str]] = None, and_or: Optional[str] = None):
        self.scans = scans
        self.operations = operations
        self.and_or = and_or
        self.output = []

    def filter_data(self):
        if self.operations:
            for operation in self.operations:
                keyword, action, value = operation.split("+")
                match action:
                    case "contains":
                        self.output.append(self._contains(keyword, value))
                    case "notcontain":
                        self.output.append(self._notcontain(keyword, value))
                    case _:
                        logger.error("Filter operation not supported!")
            # print(self.output)
            # set union or intersection
            match self.and_or:
                case "or":
                    return sorted(set().union(*self.output))
                case "and":
                    return sorted(set.intersection(*map(set, self.output)))
                case _:
                    logger.error("Logic operation not accepted!")

    def _contains(self, keyword, value):
        tmp_output = set()
        for filename, scan in self.scans.items():
            for att in fields(scan.metadata):
                if keyword in att.name and value in getattr(scan.metadata, att.name):
                    tmp_output.add(filename)
        return tmp_output

    def _notcontain(self, keyword, value):
        tmp_output = set()
        for filename, scan in self.scans.items():
            for att in fields(scan.metadata):
                if keyword in att.name and value not in getattr(scan.metadata, att.name):
                    tmp_output.add(filename)
        return tmp_output

    def _is_not(self):
        pass

    def _is(self):
        pass

    def _less_than(self):
        pass

    def _less_than_equal_to(self):
        pass

    def _greater_than(self):
        pass

    def _greater_than_equal_to(self):
        pass

    def _equal(self):
        pass

    def _not_equal(self):
        pass
