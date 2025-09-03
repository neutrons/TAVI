from dataclasses import fields
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
        scans: dict[str, Scan],
        conditions: Optional[list[Operations]] = None,
        and_or: Optional[Logic] = None,
        tol: float = 0.01,
    ):
        self.scans = scans
        self.conditions = conditions
        self.and_or = and_or
        self.tol = tol
        self.output = []

    def filter_data(self):
        if self.conditions:
            for condition in self.conditions:
                keyword, action, value = condition
                match action:
                    case Operations.CONTAINS:
                        self.output.append(self._contains(keyword, value))
                    case Operations.NOTCONTAIN:
                        self.output.append(self._notcontain(keyword, value))
                    case Operations.IS:
                        self.output.append(self._is(keyword, value))
                    case Operations.ISNOT:
                        self.output.append(self._is_not(keyword, value))
                    case Operations.EQUAL:
                        self.output.append(self._equal(keyword, value, tol=self.tol))
                    case Operations.NOTEQUAL:
                        self.output.append(self._not_equal(keyword, value, tol=self.tol))
                    case Operations.LESS:
                        self.output.append(self._less_than(keyword, value))
                    case Operations.LESSEQUAL:
                        self.output.append(self._less_than_equal_to(keyword, value))
                    case Operations.GREATER:
                        self.output.append(self._greater_than(keyword, value))
                    case Operations.GREATEREUQAL:
                        self.output.append(self._greater_than_equal_to(keyword, value))
                    case _:
                        logger.error("Filter operation not supported!")

            match self.and_or:
                case Logic.OR:
                    return sorted(set().union(*self.output))
                case Logic.AND:
                    return sorted(set.intersection(*map(set, self.output)))
                case _:
                    logger.error("Logic operation not accepted!")

    def _contains(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.CONTAINS, category="metadata")

    def _notcontain(self, keyword, value):
        return self.condition_factory(
            keyword=keyword, value=value, condition=Operations.NOTCONTAIN, category="metadata"
        )

    def _is(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.IS, category="metadata")

    def _is_not(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.ISNOT, category="metadata")

    def _equal(self, keyword, value, tol):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.EQUAL, category="data")

    def _not_equal(self, keyword, value, tol):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.NOTEQUAL, category="data")

    def _less_than(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.LESS, category="data")

    def _less_than_equal_to(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.LESSEQUAL, category="data")

    def _greater_than(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.GREATER, category="data")

    def _greater_than_equal_to(self, keyword, value):
        return self.condition_factory(keyword=keyword, value=value, condition=Operations.GREATEREUQAL, category="data")

    def condition_factory(self, keyword, value, condition, category):
        tmp_output = set()
        if category == "metadata":
            for filename, scan in self.scans.items():
                for att in fields(scan.metadata):
                    if keyword == att.name:
                        match condition:
                            case Operations.CONTAINS:
                                if value in getattr(scan.metadata, att.name):
                                    tmp_output.add(filename)
                            case Operations.NOTCONTAIN:
                                if value not in getattr(scan.metadata, att.name):
                                    tmp_output.add(filename)
                            case Operations.IS:
                                if value == getattr(scan.metadata, att.name):
                                    tmp_output.add(filename)
                            case Operations.ISNOT:
                                if value != getattr(scan.metadata, att.name):
                                    tmp_output.add(filename)
        elif category == "data":
            for filename, scan in self.scans.items():
                for att in fields(scan.data):
                    if keyword == att.name:
                        match condition:
                            case "==":
                                if (
                                    abs(value - max(getattr(scan.data, att.name))) <= self.tol
                                    and abs(value - min(getattr(scan.data, att.name))) <= self.tol
                                ):
                                    tmp_output.add(filename)
                            case "!=":
                                if (
                                    abs(value - max(getattr(scan.data, att.name))) >= self.tol
                                    and abs(value - min(getattr(scan.data, att.name))) >= self.tol
                                ):
                                    tmp_output.add(filename)
                            case ">":
                                if min(getattr(scan.data, att.name)) > value:
                                    tmp_output.add(filename)
                            case ">=":
                                if min(getattr(scan.data, att.name)) >= value:
                                    tmp_output.add(filename)
                            case "<":
                                if max(getattr(scan.data, att.name)) < value:
                                    tmp_output.add(filename)
                            case "<=":
                                if min(getattr(scan.data, att.name)) <= value:
                                    tmp_output.add(filename)

        return tmp_output
