import pandas as pd

from configs.xviewer import XViewer_Connection
from configs.xsmiles import BB_XSMILES_FILEPATH, BB_XSMILES_PLACEHOLDERS_FILEPATH
from xviewer_queries.generic_query import generic_query

def query_xsmiles(
    library_name:str = None,
    cycle_name:str = None,
    xviewer_connection:XViewer_Connection = None,
) -> pd.DataFrame:
    if xviewer_connection is not None:
        df = generic_query(
            query=f"""
            SELECT T1.ID, T5.LIBRARY_NAME, T4.CYCLE_NAME, T1.SET_ID, T3.SET_NAME, T3.SET_TYPE, 
            T1.BBID, T2.XBBID, T1.XSMILES_ORIG, T1.XSMILES_FOR_ENUM, T1.INSTRUCTION, 
            T1.XSMILES_STATUS
            FROM XCRX_CHEM.CHEM_BB_XSMILES_V2 T1
            LEFT OUTER JOIN XCRX_LIBREG.LIBREG_BBS T2 ON
                T1.BBID = T2.ID 
            LEFT OUTER JOIN XCRX_CHEM.CHEM_BB_XSMILES_V2_SETS T3 ON
                T1.SET_ID = T3.ID 
            LEFT OUTER JOIN XCRX_LIBREG.LIBREG_LIBRARY_CYCLES T4 ON
                T3.LIBRARY_CYCLE_ID = T4.ID 
            LEFT OUTER JOIN XCRX_LIBREG.LIBREG_LIBRARIES T5 ON
                T4.LIB_ID = T5.ID 
            WHERE LIBRARY_NAME = '{library_name}' AND CYCLE_NAME = '{cycle_name}'
            ORDER BY LIBRARY_NAME, CYCLE_NAME, XBBID
            """,
            connection=xviewer_connection,
        )
    else:
        df = pd.read_csv(BB_XSMILES_FILEPATH).fillna("")
        df = df[df["LIBRARY_NAME"] == library_name]
        df = df[df["CYCLE_NAME"] == cycle_name]
        df = df.reset_index(drop=True)
    return df

def query_xsmiles_placeholders(
    library_name:str = None,
    cycle_name:str = None,
    xviewer_connection:XViewer_Connection = None,
) -> pd.DataFrame:
    if xviewer_connection is not None:
        df = generic_query(
            query=f"""
            SELECT * 
            FROM XCRX_CHEM.VW_LIBRARY_XSMILES_PLACEHOLDERS
            WHERE LIBRARY_NAME = '{library_name}' AND CYCLE_NAME = '{cycle_name}'
            """,
            connection=xviewer_connection,
        )
    else:
        df = pd.read_csv(BB_XSMILES_PLACEHOLDERS_FILEPATH).fillna("")
        df = df[df["LIBRARY_NAME"] == library_name]
        df = df[df["CYCLE_NAME"] == cycle_name]
        df = df.reset_index(drop=True)
    return df
