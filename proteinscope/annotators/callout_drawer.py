"""Callout registry — manages numbered callout assignment and tracking."""

from __future__ import annotations

from core.models import NodeAnnotation


class CalloutRegistry:
    """Assigns sequential callout numbers to NodeAnnotation objects
    and maintains a mapping from number → annotation for footnote building."""

    def __init__(self) -> None:
        self._counter = 1
        self._map: dict[int, NodeAnnotation] = {}

    def register(self, annotation: NodeAnnotation) -> int:
        """Assign the next callout number to this annotation.

        Returns the assigned callout number.
        """
        n = self._counter
        annotation.callout_number = n
        self._map[n] = annotation
        self._counter += 1
        return n

    def all_annotations(self) -> list[NodeAnnotation]:
        """Return all registered annotations in callout-number order."""
        return [self._map[i] for i in sorted(self._map)]

    def __len__(self) -> int:
        return len(self._map)
