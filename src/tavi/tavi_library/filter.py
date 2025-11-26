from enum import Enum
from typing import Optional
from venv import logger

from tavi.tavi_model.FileSystem.tavi_class_factory import Scan


class Operations(Enum):
    CONTAINS = "contains"
    NOTCONTAIN = "notcontain"
    IS = "is"
    ISNOT = "isnot"
    EQUAL = "=="
    NOTEQUAL = "!="
    LESS = "<"
    LESSEQUAL = "<="
    GREATER = ">"
    GREATEREUQAL = ">="


class Logic(Enum):
    AND = "and"
    OR = "OR"


class Category(Enum):
    DATA = "data"
    METADATA = "metadata"


class Filter:
    """
    Arg:
        scans: the class that holds the loaded scan data, meta data, ubconf etc. Defined in TaviProject
        conditions: a string of filter conditions. keyword + operation + value.
                    example: "title+contains+temp", "sample_temp + < + 100". If the conditions is: contains, notcontain,
                    then we look in scans.metadata. If the conditions is: <, >, <=, >=, ==, !=, then we look in scans.data.
    """

    def __init__(
        self,
        scan_list: dict[str, Scan],
        conditions: Optional[list[Operations]] = None,
        and_or: Optional[Logic] = None,
        tol: float = 0.01,  # this can be put into a TAVI config json file as filter equal tolerance
    ):
        self.scan_list = scan_list
        self.conditions = conditions
        self.and_or = and_or
        self.tol = tol
        self.output = []

    def filter_data(self):
        """
        Decide where to look for data (data, meta data, ubconf) and perform and/or operation.
        """
        if self.conditions:
            for condition in self.conditions:
                keyword, action, value = condition
                if action in [Operations.CONTAINS, Operations.NOTCONTAIN, Operations.IS, Operations.ISNOT]:
                    self.output.append(self._condition_factory(keyword, value, action, Category.METADATA))
                else:
                    self.output.append(self._condition_factory(keyword, value, action, Category.DATA))

            match self.and_or:
                case Logic.OR:
                    return sorted(set().union(*self.output))
                case Logic.AND:
                    return sorted(set.intersection(*map(set, self.output)))
                case _:
                    logger.error("Logic operation not accepted!")

    def _condition_factory(self, keyword, value, condition, category):
        """
        Abstract factory that returns the filtered results based on the keyword, value, condition and category.
        """
        tmp_output = set()
        if category == Category.METADATA:
            for filename, scan in self.scan_list.items():
                if hasattr(scan.metadata, keyword):
                    att = scan.metadata
                elif hasattr(scan.ubconf, keyword):
                    att = scan.ubconf
                else:
                    logger.log("No matching entry with", keyword)
                # search meta data and ubconf
                match condition:
                    case Operations.CONTAINS:
                        if value in getattr(att, keyword):
                            tmp_output.add(filename)
                    case Operations.NOTCONTAIN:
                        if value not in getattr(att, keyword):
                            tmp_output.add(filename)
                    case Operations.IS:
                        if value == getattr(att, keyword):
                            tmp_output.add(filename)
                    case Operations.ISNOT:
                        if value != getattr(att, keyword):
                            tmp_output.add(filename)
        elif category == Category.DATA:
            value = float(value)
            for filename, scan in self.scan_list.items():
                if not hasattr(scan.data, keyword):
                    logger.log("No matching entry with", keyword)
                else:
                    match condition:
                        case Operations.EQUAL:
                            if value + self.tol >= max(getattr(scan.data, keyword)) and value - self.tol <= min(
                                getattr(scan.data, keyword)
                            ):
                                tmp_output.add(filename)
                        case Operations.NOTEQUAL:
                            if value + self.tol < min(getattr(scan.data, keyword)) or value - self.tol > max(
                                getattr(scan.data, keyword)
                            ):
                                tmp_output.add(filename)
                        case Operations.GREATER:
                            if min(getattr(scan.data, keyword)) > value:
                                tmp_output.add(filename)
                        case Operations.GREATEREUQAL:
                            if min(getattr(scan.data, keyword)) >= value:
                                tmp_output.add(filename)
                        case Operations.LESS:
                            if max(getattr(scan.data, keyword)) < value:
                                tmp_output.add(filename)
                        case Operations.LESSEQUAL:
                            if min(getattr(scan.data, keyword)) <= value:
                                tmp_output.add(filename)
        return tmp_output
