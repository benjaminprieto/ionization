"""
base.py - Base interface for protonation engines
==================================================
Abstract class that defines the contract every engine must follow.

To add a new engine:
    1. Create a new file in engines/ (e.g. qupkake_engine.py)
    2. Create a class that inherits from BaseEngine
    3. Implement all abstract methods
    4. Register it in ionizer.py ENGINE_REGISTRY

Project: ionprofile
"""

from abc import ABC, abstractmethod
from typing import Optional


class BaseEngine(ABC):
    """
    Abstract base class for protonation engines.

    Every engine must answer two questions:
        1. What is the formal charge of this molecule at pH X?
        2. What does this molecule look like (SMILES) at pH X?
    """

    @abstractmethod
    def calculate_charge_at_ph(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[int]:
        """Calculate formal charge at a specific pH."""
        pass

    @abstractmethod
    def get_protonated_smiles(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[str]:
        """Get the protonated SMILES at a specific pH."""
        pass

    @abstractmethod
    def is_available(self) -> bool:
        """Check if this engine's dependencies are installed."""
        pass

    @abstractmethod
    def name(self) -> str:
        """Return engine display name (e.g. 'Dimorphite-DL')."""
        pass
