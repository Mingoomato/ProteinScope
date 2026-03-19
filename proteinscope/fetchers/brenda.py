"""BRENDA enzyme database fetcher for temperature and pH optima.

Uses BRENDA SOAP API with credentials stored in BRENDA_USERNAME / BRENDA_PASSWORD.
Falls back to UniProt BIOPHYSICOCHEMICAL PROPERTIES comment when BRENDA is unavailable
or when the protein is not an enzyme.

Register at: https://www.brenda-enzymes.org/register.php
"""

from __future__ import annotations

import hashlib
import os
import re
from typing import Optional

import httpx

BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_server.php"


def _brenda_creds() -> tuple[Optional[str], Optional[str]]:
    return os.getenv("BRENDA_USERNAME"), os.getenv("BRENDA_PASSWORD")


def _sha256(password: str) -> str:
    return hashlib.sha256(password.encode()).hexdigest()


async def fetch_brenda_temperature(ec_number: str, organism: str) -> Optional[float]:
    """Query BRENDA SOAP for temperature optimum. Returns degrees Celsius or None."""
    username, password = _brenda_creds()
    if not username or not password:
        return None
    try:
        # Build a minimal SOAP envelope
        password_hash = _sha256(password)
        soap_body = f"""<?xml version="1.0" encoding="UTF-8"?>
<soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
                  xmlns:bren="https://www.brenda-enzymes.org">
  <soapenv:Body>
    <bren:getTemperatureOptimum>
      <bren:parameters>{username},{password_hash},{ec_number},*,*,{organism},*</bren:parameters>
    </bren:getTemperatureOptimum>
  </soapenv:Body>
</soapenv:Envelope>"""
        headers = {
            "Content-Type": "text/xml; charset=utf-8",
            "SOAPAction": "",
        }
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.post(BRENDA_WSDL, content=soap_body.encode(), headers=headers)
            r.raise_for_status()
            # Parse temperature value from XML response
            match = re.search(r"<return>(.*?)</return>", r.text, re.DOTALL)
            if match:
                content = match.group(1)
                # BRENDA returns entries like "temperatureOptimum*45#..."
                temp_match = re.search(r"temperatureOptimum\*(\d+(?:\.\d+)?)", content)
                if temp_match:
                    return float(temp_match.group(1))
    except Exception:
        pass
    return None


async def fetch_brenda_ph(ec_number: str, organism: str) -> Optional[float]:
    """Query BRENDA SOAP for pH optimum. Returns pH value or None."""
    username, password = _brenda_creds()
    if not username or not password:
        return None
    try:
        password_hash = _sha256(password)
        soap_body = f"""<?xml version="1.0" encoding="UTF-8"?>
<soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
                  xmlns:bren="https://www.brenda-enzymes.org">
  <soapenv:Body>
    <bren:getPhOptimum>
      <bren:parameters>{username},{password_hash},{ec_number},*,*,{organism},*</bren:parameters>
    </bren:getPhOptimum>
  </soapenv:Body>
</soapenv:Envelope>"""
        headers = {
            "Content-Type": "text/xml; charset=utf-8",
            "SOAPAction": "",
        }
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.post(BRENDA_WSDL, content=soap_body.encode(), headers=headers)
            r.raise_for_status()
            match = re.search(r"<return>(.*?)</return>", r.text, re.DOTALL)
            if match:
                ph_match = re.search(r"phOptimum\*(\d+(?:\.\d+)?)", match.group(1))
                if ph_match:
                    return float(ph_match.group(1))
    except Exception:
        pass
    return None


def extract_temperature_from_uniprot(entry: dict) -> Optional[float]:
    """Fallback: parse temperature from UniProt BIOPHYSICOCHEMICAL PROPERTIES comment."""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "BIOPHYSICOCHEMICAL PROPERTIES":
            temp_deps = comment.get("temperatureDependence", {})
            texts = temp_deps.get("texts", [])
            for t in texts:
                val = t.get("value", "")
                match = re.search(r"optimum.*?(\d+(?:\.\d+)?)\s*[°℃C]", val, re.IGNORECASE)
                if match:
                    return float(match.group(1))
    return None


def extract_ph_from_uniprot(entry: dict) -> Optional[float]:
    """Fallback: parse pH from UniProt BIOPHYSICOCHEMICAL PROPERTIES comment."""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "BIOPHYSICOCHEMICAL PROPERTIES":
            ph_dep = comment.get("phDependence", {})
            texts = ph_dep.get("texts", [])
            for t in texts:
                val = t.get("value", "")
                match = re.search(r"optimum.*?pH\s*(\d+(?:\.\d+)?)", val, re.IGNORECASE)
                if match:
                    return float(match.group(1))
    return None


def get_ec_number(uniprot_entry: dict) -> Optional[str]:
    """Extract EC number from UniProt entry cross-references."""
    desc = uniprot_entry.get("proteinDescription", {})
    for name_block in [
        desc.get("recommendedName", {}),
        *desc.get("alternativeNames", []),
    ]:
        for ec in name_block.get("ecNumbers", []):
            val = ec.get("value", "")
            if val:
                return val
    return None
