"""
ionizer.py - Engine dispatcher
================================
Selects the configured protonation engine and returns it ready to use.

The engines live in profiling/engines/. This file only knows about
the registry — it doesn't contain any protonation logic itself.

Usage:
    engine = get_engine("dimorphite")
    charge = engine.calculate_charge_at_ph("CCO", ph=7.4)

Adding a new engine:
    1. Create the engine class in engines/ (inherit from BaseEngine)
    2. Add one line to ENGINE_REGISTRY below

Project: ionprofile
"""

import logging
from typing import Dict

from ionprofile.profiling.engines.base import BaseEngine
from ionprofile.profiling.engines.dimorphite_engine import DimorphiteEngine
from ionprofile.profiling.engines.openbabel_engine import OpenBabelEngine
from ionprofile.profiling.engines.qupkake_engine import QupKakeEngine
from ionprofile.profiling.rdkit_utils import is_rdkit_available

logger = logging.getLogger(__name__)


# =========================================================================
# ENGINE REGISTRY
# =========================================================================
# To add a new engine, import it above and add one line here.

ENGINE_REGISTRY: Dict[str, type] = {
    "dimorphite": DimorphiteEngine,   # fast empirical, ms/mol
    "openbabel": OpenBabelEngine,     # fast empirical, ms/mol
    "qupkake": QupKakeEngine,         # ML+QM pKa, RMSE 0.5-0.8, min/mol
}


def get_engine(engine_name: str = "dimorphite") -> BaseEngine:
    """
    Get the configured protonation engine.

    Args:
        engine_name: Key from ENGINE_REGISTRY (e.g. "dimorphite").

    Returns:
        Initialized engine instance, ready to use.

    Raises:
        ValueError: If engine_name is not registered.
    """
    if engine_name not in ENGINE_REGISTRY:
        available = ", ".join(sorted(ENGINE_REGISTRY.keys()))
        raise ValueError(
            f"Unknown engine: '{engine_name}'. "
            f"Available: {available}"
        )

    engine = ENGINE_REGISTRY[engine_name]()

    if not engine.is_available():
        logger.warning(
            f"Engine '{engine.name()}' dependencies not installed. "
            f"Charge calculations will return None."
        )

    logger.info(f"Protonation engine: {engine.name()} "
                f"(available: {engine.is_available()})")

    return engine


def check_dependencies(engine_name: str = "dimorphite") -> dict:
    """
    Check all dependencies for a given engine.

    Returns:
        Dict with availability status.
    """
    if engine_name not in ENGINE_REGISTRY:
        return {"engine": engine_name, "available": False,
                "error": "Unknown engine"}

    instance = ENGINE_REGISTRY[engine_name]()
    return {
        "engine": engine_name,
        "engine_name": instance.name(),
        "rdkit": is_rdkit_available(),
        "engine_available": instance.is_available(),
        "fully_operational": instance.is_available(),
    }


def list_engines() -> dict:
    """
    List all registered engines and their status.

    Returns:
        Dict of engine_name -> {name, available}
    """
    result = {}
    for key, engine_class in ENGINE_REGISTRY.items():
        instance = engine_class()
        result[key] = {
            "name": instance.name(),
            "available": instance.is_available(),
        }
    return result