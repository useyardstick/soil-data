from enum import Enum
from typing import TypeVar

# TODO: replace this with typing.Self after upgrading to python 3.11
_DepthType = TypeVar("_DepthType", bound="DepthEnum")


class DepthEnum(Enum):
    """
    Provides helpful methods for enumerating depth ranges.

    Values should be (start_depth, end_depth) 2-tuples.
    """

    @classmethod
    def select_between(
        cls: type[_DepthType], start_depth: int, end_depth: int
    ) -> list[_DepthType]:
        start_depths = {depth.start_depth for depth in cls}
        if start_depth not in start_depths:
            raise Exception(f"start_depth {start_depth} must be one of {start_depths}")

        end_depths = {depth.end_depth for depth in cls}
        if end_depth not in end_depths:
            raise Exception(f"end_depth {end_depth} must be one of {end_depths}")

        return cls.select_including(start_depth, end_depth)

    @classmethod
    def select_including(
        cls: type[_DepthType], start_depth: int, end_depth: int
    ) -> list[_DepthType]:
        max_depth = max(depth.end_depth for depth in cls)
        if start_depth < 0 or end_depth > max_depth:
            raise Exception(f"Maximum depth range: 0 - {max_depth}")

        if end_depth <= start_depth:
            raise Exception(
                f"end_depth {end_depth} must be greater than start_depth {start_depth}"
            )

        selected_depths = (
            depth
            for depth in cls
            if start_depth < depth.end_depth and end_depth > depth.start_depth
        )
        return sorted(selected_depths, key=lambda depth: depth.start_depth)

    @property
    def start_depth(self) -> int:
        start_depth, _ = self.value
        return start_depth

    @property
    def end_depth(self) -> int:
        _, end_depth = self.value
        return end_depth

    @property
    def thickness(self) -> int:
        start_depth, end_depth = self.value
        return end_depth - start_depth
